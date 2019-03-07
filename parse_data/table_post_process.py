import re
import random
import itertools
import functools
import pprint
from collections import Counter

from nltk.tokenize.punkt import PunktSentenceTokenizer, PunktParameters
import fitz
import cv2

from .utils import overlap_ratio, clean_text

CAPTION_PATTERN = (r'^((Supp(\.)?|Supplementa(l|ry))\s*)?((T|t)able|TABLE)\s*'
                   r'S?([0-9]+|I(?=[^I]|$)|II(?=[^I]|$)|III|IV|V|VI(?=[^I]|$)|VII(?=[^I]|$)|VIII|IX|X'
                   r'|A|B|C|D|E|F|G|H|I|J|K|L)')

punkt_param = PunktParameters()
punkt_param.abbrev_types = set(['fig'])
sent_tokenizer = PunktSentenceTokenizer(punkt_param)


def box_key(box, page_width):
    """key for sort table boxes
    """
    return (box[0] >= page_width / 2.5, box[1], box[0])


def get_caption(block):
    """check if the block is caption of a table
    """
    if block['type'] == 1:
        return None
    bbox = block['bbox']

    direction = None
    lines = []
    for line in block['lines']:
        direction = line['dir']
        lines.append(map(lambda x: x['text'].strip(), line['spans']))
    text = ' '.join(filter(bool, map(clean_text, itertools.chain(*lines))))

    match = re.match(CAPTION_PATTERN, text)
    if not match:
        return None

    #FIXME bug:wrong line
    flags = line['spans'][0]['flags']

    if re.match(CAPTION_PATTERN + r'\s*([:.]|$)', text):
        label = match.group(0)
    elif (flags & 2) or (flags & 16):
        label = match.group(0)
    elif len(sent_tokenizer.tokenize(text)) == 1:
        label = match.group(0)
    else:
        return None

    caption = {
        'text': text,
        'bbox': bbox,
        'label': label,
        'dir': direction,
    }
    return caption


def find_captions(page_dict):
    """find captions in the page
    """
    captions = []
    for block in page_dict['blocks']:
        caption = get_caption(block)
        if caption:
            captions.append(caption)
    return captions


def rotate(box, direct):
    """rotate box by direct
    Args:
        box: [x1, y1, x2, y2]
        direct: (cos(theta), sin(theta))

    Returns:
        rotated [x1, y1, x2, y2]
    """
    if direct == (0, -1):
        return [-box[3], box[0], -box[1], box[2]]
    return box


def caption_distance(table_box, caption_box, direct):
    """calcuate distance between given table and caption
    """
    if direct == (0, -1):
        table_box = rotate(table_box, direct)
        caption_box = rotate(caption_box, direct)

    c_center = (caption_box[0] + caption_box[2]) / 2
    c_bottom = caption_box[3]

    t_left = table_box[0]
    t_right = table_box[2]
    t_top = table_box[1]

    dx = max([0, t_left - c_center, c_center - t_right])
    dy = max([0, c_bottom - t_top])
    dist = dx * dx + dy * dy
    return dist


def merge_boxes(box1, box2):
    """merge two boxes
    """
    return (min(box1[0], box2[0]), min(box1[1], box2[1]),
            max(box1[2], box2[2]), max(box1[3], box2[3]))


def merge_tables(table_boxes):
    """merge tables if they overlap
    """
    classes = [i for i in range(len(table_boxes))]
    for i, table_box in enumerate(table_boxes):
        for j, table_box2 in enumerate(table_boxes):
            r = overlap_ratio(table_box, table_box2, extend=20)
            if r >= 0.01:
                classes = [classes[i] if k == classes[j] else k for k in classes]

    ret = []
    for cls in set(classes):
        table_boxes_cls = [table_boxes[i] for i, c in enumerate(classes) if c == cls]
        box = functools.reduce(merge_boxes, table_boxes_cls)
        ret.append(box)
    return ret


def merge_lines(lines):
    """merge lines in the pdf
    """
    def _merge(l1, l2):
        return {
            'dir': l1['dir'],
            'spans': l1['spans'] + l2['spans'],
            'bbox': merge_boxes(l1['bbox'], l2['bbox']),
        }

    groups = [i for i in range(len(lines))]
    for i, line1 in enumerate(lines):
        for j, line2 in enumerate(lines):
            if groups[i] == groups[j]:
                continue
            if (line1['dir'] == line2['dir'] and
                    overlap_ratio(line1['bbox'], line2['bbox'], extend=5)):
                groups = [groups[i] if g == groups[j] else g for g in groups]

    ret = []
    for cls in set(groups):
        group_cls = [lines[i] for i, c in enumerate(groups) if c == cls]
        line = functools.reduce(_merge, group_cls)
        ret.append(line)
    return ret


def adjust_table_box(table_box, lines):
    """adjust table_box by the blocks from pdf
    """
    bboxes = []
    for block in lines:
        block_box = block['bbox']
        if overlap_ratio(table_box, block_box) >= 0.3:
            bboxes.append(block_box)
    if not bboxes:
        return None
    box = functools.reduce(merge_boxes, bboxes)
    return box


def adjust_tables(page_dict, table_boxes):
    """adjust table_boxes by the blocks from pdf
    """
    lines = []
    for block in page_dict['blocks']:
        if block['type'] == 0:
            lines += block['lines']
    lines = list(filter(lambda line: any(span['text'].strip() for span in line['spans']), lines))

    table_boxes = [adjust_table_box(table_box, lines) for table_box in table_boxes]
    table_boxes = list(filter(bool, table_boxes))
    table_boxes = merge_tables(table_boxes)

    prev_len = 0
    while prev_len != len(lines):
        prev_len = len(lines)
        lines = merge_lines(lines)

    prev = None
    while prev != table_boxes:
        prev = table_boxes
        table_boxes = [adjust_table_box(table_box, lines) for table_box in table_boxes]
        table_boxes = list(filter(bool, table_boxes))
    table_boxes = merge_tables(table_boxes)

    page_width = page_dict['width']
    table_boxes.sort(key=functools.partial(box_key, page_width=page_width))
    return table_boxes, lines


def associate_caption_table(table_boxes, captions, last_page_caption):
    """assign tables to captions
    """
    dists = []
    for i, table_box in enumerate(table_boxes):
        for j, caption in enumerate(captions):
            direct = caption['dir']
            dists.append((caption_distance(table_box, caption['bbox'], direct), i, j))
    dists.sort()
    used_table, used_caption = set(), set()
    caption_assignment = [None for _ in table_boxes]
    for _, i, j in dists:
        if i in used_table or j in used_caption:
            continue
        caption_assignment[i] = j
        used_table.add(i)
        used_caption.add(j)

    is_continued = [(idx is None) for idx in caption_assignment]

    prev_caption_id = -1 if last_page_caption else None
    for i in range(len(table_boxes)):
        if caption_assignment[i] is None:
            caption_assignment[i] = prev_caption_id
        prev_caption_id = caption_assignment[i]

    ret = []
    for i, table_box in enumerate(table_boxes):
        if caption_assignment[i] is not None:
            j = caption_assignment[i]
            caption = last_page_caption if j == -1 else captions[j]
            ret.append({
                'bbox': table_boxes[i],
                'caption': caption,
                'continued': is_continued[i],
            })
    return ret


def table_post_process(page_dict, table_boxes, last_page_caption):
    """post process for table detection
    """
    table_boxes, lines = adjust_tables(page_dict, table_boxes)
    captions = find_captions(page_dict)
    tables = associate_caption_table(table_boxes, captions, last_page_caption)
    return tables
