"""main"""
import os
from os.path import basename
import gzip
import json
import argparse
import time
import logging
import multiprocessing
import pickle

from . import pytmvar

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def process(pmid, body, tables, var_extr):
    var_body_mentions = var_extr.extract(body, pmid)

    var_table_mentions = []
    for k, t in enumerate(tables):
        for i, row in enumerate(t['cells']):
            for j, cell in enumerate(row):
                mentions = var_extr.extract(cell['text'], pmid)
                for text, (start, end), var in mentions:
                    var_table_mentions.append((text, (k, i, j), (start, end), var))
    return var_body_mentions, var_table_mentions
