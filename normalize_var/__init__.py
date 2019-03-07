# pylint: disable=protected-access
"""normalize hgvs/rsid to chromosome variants
"""
import sys
import os
import argparse
import multiprocessing
import json
import time
import logging
from typing import NamedTuple

from sqlalchemy import create_engine

sys.path.insert(0, '/app/mysqldb')
from models import (engine,  # pylint: disable=no-name-in-module, wrong-import-position
                    var_pmid)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class IndexKey(NamedTuple):
    """key of the paper indexing
    """
    chrom: str
    position: int
    ref: str
    alt: str
    paper_id: str
    file_idx: int
    in_table: bool
    table_idx: int
    row: int
    col: int
    start: int


def write_mysql(_id, values):
    """write data into mysql db
    """
    with engine.connect() as conn:
        conn.execute(var_pmid.delete().where(var_pmid.c._id == _id))
        for i in range(0, len(values), 1000):
            batch = values[i:i + 1000]
            conn.execute(var_pmid.insert(), batch)  # pylint: disable=no-value-for-parameter


def process(results, _id, var_normalizer):  # pylint: disable=too-many-locals
    """normalize HGVS variants to chromosome variants
    """
    vcf_dict = dict()

    with var_normalizer.connect():
        for r in results:
            for tx_name, vcf in var_normalizer.to_vcf([r.gene_id], r.var_json):
                chrom, position, ref, alt = vcf
                if len(ref) > 200 or len(alt) > 200:
                    continue
                key = IndexKey(
                    chrom=chrom, position=position, ref=ref, alt=alt,
                    paper_id=_id, file_idx=r.file_idx,
                    in_table=r.in_table, table_idx=r.table_idx,
                    row=r.row, col=r.col, start=r.start,
                )
                vcf_dict[key] = (tx_name, r.gene_id, r.end, r.variant, r.gene_start, r.gene_end)

    values = []
    for key, (tx_name, gene_id, end, mention, gene_start, gene_end) in vcf_dict.items():
        values.append({
            'chrom': key.chrom,
            'position': key.position,
            'ref': key.ref,
            'alt': key.alt,
            'gene_id': gene_id,
            'transcript': tx_name,
            '_id': key.paper_id,
            'file_idx': key.file_idx,
            'in_table': key.in_table,
            'table_idx': key.table_idx,
            'row': key.row,
            'col': key.col,
            'start': key.start,
            'end': end,
            'mention': mention,
            'gene_start': gene_start,
            'gene_end': gene_end,
        })
    write_mysql(_id, values)
