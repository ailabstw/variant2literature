"""pdf utils
"""
import os
import itertools
import subprocess
import logging

from nltk.tokenize.punkt import PunktSentenceTokenizer, PunktParameters
import numpy as np
import fitz
import rpyc
import cv2

from .utils import clean_text, overlap_ratio, load_np
from .table_post_process import table_post_process

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CHAR_MAP = {
    'AdvPS586B': str.maketrans('.12', '>+-'),
    'AdvPSMP4': str.maketrans('[', '>'),
    'AdvPi1': str.maketrans('4', '>'),
    'AdvPi2': str.maketrans('?', '>'),
    'AdvP7DA6': str.maketrans('.', '>'),
    'AdvP4C4E74': str.maketrans('!' + b'\xc3\xbe'.decode('utf8'), '>+'),
    'AdvMPi-One': str.maketrans('1', '+'),
    'AdvMacMthSyN': str.maketrans('\x00', '-'),
    'GeneralSymbolsPB': str.maketrans(']', '>'),
    'Universal-GreekwithMathP': str.maketrans(b'\xef\xbf\xbd'.decode('utf8'), '>'),
}


def pdf_to_image(filename):
    """convert pdf to image
    """
    with open(os.devnull, 'w') as fnull:
        data = subprocess.check_output(['pdftoppm', '-jpeg', '-r', '150', filename], stderr=fnull)
    data = list(map(lambda x: x + b'\xff\xd9', data.split(b'\xff\xd9')[:-1]))
    images = []
    for page_data in data:
        np_arr = np.fromstring(page_data, np.uint8)
        image = cv2.imdecode(np_arr, cv2.IMREAD_COLOR)
        images.append(image)
    return images, data


def get_lines(block):
    """get text lines from the pdf block
    """
    if block['type'] == 1:
        return ''

    lines = []
    for line in block['lines']:
        line_texts, prev = [], dict()
        for span in line['spans']:
            text = span['text'].translate(CHAR_MAP.get(span['font'], str.maketrans('', '')))
            if prev and prev['size'] != span['size']:
                text = ' ' + text
            line_texts.append(text)
            prev = span
        lines.append(''.join(line_texts))
    return lines


def get_text_dicts(blocks):
    """collect all texts and their bounding boxes from blocks
    """

    text_dicts = []
    for block in blocks:
        if block['type'] == 1:
            continue
        for line in block['lines']:
            direction = line['dir']
            line_texts = []
            for span in line['spans']:
                text = span['text'].translate(CHAR_MAP.get(span['font'], str.maketrans('', '')))
                line_texts.append(text)
            line_text = ''.join(line_texts).strip()
            if line_text:
                text_dicts.append({
                    'orig_bbox': line['bbox'],
                    'bbox': rotate(line['bbox'], direction),
                    'text': ''.join(line_texts),
                })
    text_dicts.sort(key=lambda x: (x['bbox'][1], x['bbox'][0]))
    return text_dicts


class Rectangle:
    """text and bounding box
    """
    def __init__(self, bbox, text):
        self.bbox = bbox
        self.text = text

    def __repr__(self):
        return '[' + ', '.join(map('{:.3f}'.format, self.bbox)) + ']'

    def h_overlap(self, other):
        """horizontal overlap
        """
        if not isinstance(other, Rectangle):
            raise TypeError('should be Rectangle')
        min_r = min(self.bbox[2], other.bbox[2])
        max_l = max(self.bbox[0], other.bbox[0])
        return max(0, min_r - max_l)

    def v_overlap(self, other):
        """vertical overlap
        """
        if not isinstance(other, Rectangle):
            raise TypeError('should be Rectangle')
        min_b = min(self.bbox[3], other.bbox[3])
        max_t = max(self.bbox[1], other.bbox[1])
        return max(0, min_b - max_t)

    def v_overlap_ratio(self, other):
        """vertical overlap ratio
        """
        d = min(self.bbox[3] - self.bbox[1], other.bbox[3] - other.bbox[1])
        return self.v_overlap(other) / d

    def union(self, other):
        """union two rectangle
        """
        self.text = self.text + ' ' + other.text
        self.bbox = (
            min(self.bbox[0], other.bbox[0]),
            min(self.bbox[1], other.bbox[1]),
            max(self.bbox[2], other.bbox[2]),
            max(self.bbox[3], other.bbox[3]),
        )

    def intersect(self, other):
        """intersect two rectangle
        """
        self.text = self.text + ' ' + other.text
        self.bbox = (
            max(self.bbox[0], other.bbox[0]),
            max(self.bbox[1], other.bbox[1]),
            min(self.bbox[2], other.bbox[2]),
            min(self.bbox[3], other.bbox[3]),
        )


def group_by_lines(text_dicts):
    """group texts by lines
    """
    lines = []
    r_line = Rectangle([0, 0, 0, 0], '')
    for text_dict in text_dicts:
        r = Rectangle(text_dict['bbox'], text_dict['text'])
        if not lines or r_line.v_overlap_ratio(r) < 0.1:
            # if lines:
            #     print(r_line.text)
            lines.append([r])
            r_line = Rectangle(r.bbox, text_dict['text'])
        else:
            lines[-1].append(r)
            r_line.intersect(r)
    return lines


def find_column_positions(lines):
    """determine column positions
    """
    rights = []
    n_lines = len(lines)
    tcs = list(itertools.chain(*lines))
    tcs.sort(key=lambda x: x.bbox[0])

    i, cur, cur_tcs = 0, 0, []
    while i < len(tcs):
        if cur < n_lines * 0.8:
            cur_tcs.append(tcs[i])
            cur += 1
            i += 1

        if cur >= n_lines * 0.8 or i == len(tcs):
            cur_tcs.sort(key=lambda x: x.bbox[2])
            rights.append(cur_tcs[len(cur_tcs) // 4 * 3].bbox[2])
            cur, cur_tcs = 0, []
            while i > 0 and i < len(tcs) and tcs[i].bbox[2] - rights[-1] > rights[-1] - tcs[i].bbox[0]:
                i -= 1
            while i < len(tcs) and tcs[i].bbox[2] - rights[-1] < rights[-1] - tcs[i].bbox[0]:
                i += 1
    return rights


def find_row_positions(lines):
    """determine row positions
    """
    bottoms = []
    for line in lines:
        bottoms.append(max([r.bbox[3] for r in line]))
    bottoms = list(set(bottoms))
    bottoms.sort()
    return bottoms


def find_cell(bbox, rows, columns):
    """determine which cell the box belongs to
    """
    # choose which the bottom of the bounding box in
    r = len(rows) - 1
    for i, x in enumerate(rows):
        if bbox[3] <= x:
            r = i
            break

    # choose max horizontal overlap
    columns = [0] + columns
    c, max_overlap = 0, 0
    for i in range(1, len(columns)):
        if min(bbox[2], columns[i]) - max(bbox[0], columns[i - 1]) > max_overlap:
            max_overlap = min(bbox[2], columns[i]) - max(bbox[0], columns[i - 1])
            c = i - 1
    return r, c


def rotate(bbox, direction):
    """return rotated box coordinates
    """
    if direction == (0.0, -1.0):
        return [-bbox[3], bbox[0], -bbox[1], bbox[2]]
    return bbox


def aggregate_cell(dicts):
    """merge blocks in a cell
    """
    if not dicts:
        return '', (0.0, 0.0, 0.0, 0.0)

    dicts = sorted(dicts, key=lambda x: x['orig_bbox'][0])

    r = Rectangle(dicts[0]['orig_bbox'], '')
    for d in dicts:
        r.union(Rectangle(d['orig_bbox'], d['text']))

    text = clean_text(r.text)
    bbox = r.bbox
    return text, bbox


def construct_table(blocks):
    """construct a table from given pdf blocks
    """
    text_dicts = get_text_dicts(blocks)
    lines = group_by_lines(text_dicts)

    rows = find_row_positions(lines)
    n_rows = len(rows)

    columns = find_column_positions(lines)
    n_cols = len(columns)

    table_cells = [[[] for col in columns] for row in rows]
    for text_dict in text_dicts:
        r, c = find_cell(text_dict['bbox'], rows, columns)
        table_cells[r][c].append(text_dict)

    table = []
    for r in range(n_rows):
        row = []
        for c in range(n_cols):
            text, bbox = aggregate_cell(table_cells[r][c])
            row.append({
                'text': text,
                'bbox': bbox,
            })
        table.append(row)
    return table


def get_pdf_page_dict(page, ratio):
    """get dictionary of the page and adjust bounding boxes
    """
    page_dict = page.getText(output='dict')

    page_dict['width'] *= ratio
    page_dict['height'] *= ratio

    for block in page_dict['blocks']:
        if block['type'] == 0:
            block['bbox'] = [int(z * ratio) for z in block['bbox']]
            for line in block['lines']:
                line['bbox'] = [int(z * ratio) for z in line['bbox']]

        elif block['type'] == 1:
            block['width'] *= ratio
            block['height'] *= ratio
            block['bbox'] = [int(z * ratio) for z in block['bbox']]
    return page_dict


def split_sents(body, max_sent_len=2000):
    """split sentence if too long
    """
    ret = []
    for sent in body:
        while len(sent) > max_sent_len:
            ret.append(sent[:max_sent_len])
            sent = sent[max_sent_len:]
        ret.append(sent)
    return ret


def get_pdf_objects(filename, table_detect=True):  # pylint: disable=too-many-locals
    """extract body, table, table images from pdf
    """
    body, tables = [], []

    pages = fitz.open(filename)
    page_images, page_image_data = pdf_to_image(filename)

    prev_caption = None
    for i, page in enumerate(pages):
        ratio = page_images[i].shape[0] / page.rect[3]

        page_dict = get_pdf_page_dict(page, ratio)

        pred_table_boxes = find_tables(page_image_data[i]) if table_detect else []
        page_tables = table_post_process(page_dict, pred_table_boxes, prev_caption)
        prev_caption = page_tables[-1]['caption'] if page_tables else None

        # seperate body blocks and table blocks
        table_blocks = [[] for _ in page_tables]

        for block in page_dict['blocks']:
            if block['type'] == 1:
                continue
            for j, table in enumerate(page_tables):
                if (not table['continued'] and
                        overlap_ratio(block['bbox'], table['caption']['bbox']) > 0.5):
                    break
                elif overlap_ratio(block['bbox'], table['bbox']) > 0.5:
                    table_blocks[j].append(block)
                    break
            else:
                body += get_lines(block)

        # construct table
        for j, (blocks, table) in enumerate(zip(table_blocks, page_tables)):
            table['cells'] = construct_table(blocks)

        # crop table images
        for table in page_tables:
            x1, y1, x2, y2 = table['bbox']
            image = page_images[i][y1:y2, x1:x2, :]
            if image.size == 0:
                continue
            image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
            img_data = cv2.imencode('.jpg', image, [int(cv2.IMWRITE_JPEG_QUALITY), 75])[1].tostring()
            table['image'] = img_data

        tables += page_tables

    # sentence tokenize body text
    body = ' '.join(map(clean_text, body))
    punkt_param = PunktParameters()
    punkt_param.abbrev_types = set(['fig'])
    body = list(PunktSentenceTokenizer(punkt_param).tokenize(body))
    body = split_sents(body)

    return body, tables


def find_tables(img_data):
    """get table predictions
    """
    config = {'allow_all_attrs': True, 'sync_request_timeout': None}
    host = os.environ['LOAD_BALANCER_HOST']

    with rpyc.connect(host=host, port=18861, config=config) as conn:
        ret = conn.root.detect(img_data)
        tables = load_np(ret)
    return tables
