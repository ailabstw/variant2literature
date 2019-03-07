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

from . import pygnormplus

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def process(pmid, body, tables, gene_extr):
    gene_body_mentions = gene_extr.extract(body, pmid)

    gene_caption_mentions = []
    gene_table_mentions = []
    for k, t in enumerate(tables):
        if 'caption' in t:
            caption = t['caption']['text']
            for text, (start, end), gene in gene_extr.extract(caption, pmid):
                gene_caption_mentions.append((text, k, (start, end), gene))

        for i, row in enumerate(t['cells']):
            for j, cell in enumerate(row):
                mention = gene_extr.search_gene(cell['text'])
                if not mention:
                    continue
                text, (start, end), gene = mention
                gene_table_mentions.append((text, (k, i, j), (start, end), gene))

    return (gene_body_mentions,
            gene_caption_mentions,
            gene_table_mentions)
