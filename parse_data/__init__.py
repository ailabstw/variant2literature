"""parse paper data
"""
import sys
import argparse
import time
import string
import tempfile
import os
import io
import multiprocessing
import logging
import json
import pickle
import gzip
import zlib
import hashlib
import base64
import shutil

import cv2

from .parse import parse_dir, PaperData

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_path(_id):
    """get path to put paper data
    """
    dirname = '/paper_data'
    m = hashlib.md5(_id.encode('ascii'))
    h = base64.b16encode(m.digest()).decode('ascii')
    h1, h2 = h[:2], h[2:4]
    path = os.path.join(dirname, h1, h2)
    return path


def dump_data(_id, parsed_data):
    """dump parsed data
    """
    output_dir = os.path.join(get_path(_id), _id)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    for idx, filename, data in parsed_data:
        os.makedirs(f'{output_dir}/{idx}/table_images')

        with open(f'{output_dir}/{idx}/filename.txt', 'w') as fout:
            fout.write(filename)

        with open(f'{output_dir}/{idx}/body.txt', 'w') as fout:
            fout.write(data.body)

        for i, table in enumerate(data.tables):
            if 'image' in table:
                with open(f'{output_dir}/{idx}/table_images/{i}.jpg', 'wb') as fout:
                    fout.write(table['image'])
                del table['image']

        with open(f'{output_dir}/{idx}/tables.json', 'w') as fout:
            fout.write(json.dumps(data.tables))


def truncate_data(parsed_data, max_body_len=100000):
    """truncate if too long (usually not real articles)
    """
    ret = []
    for idx, filename, data in parsed_data:
        if len(data.body) > max_body_len:
            data = PaperData(data.body[:max_body_len], data.tables)
        ret.append((idx, filename, data))
    return ret


def process(_id, dir_path, nxml_only=False, table_detect=True, save_data=False):
    """parse the files in the directory
    """
    parsed_data = parse_dir(dir_path,
                            nxml_only=nxml_only,
                            table_detect=table_detect)
    parsed_data = list(parsed_data)
    parsed_data = truncate_data(parsed_data)

    if save_data:
        dump_data(_id, parsed_data)

    return parsed_data
