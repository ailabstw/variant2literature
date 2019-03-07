# pylint: disable=no-value-for-parameter
"""retrive record from mysql
"""
import sys
import logging
from itertools import starmap

from sqlalchemy.sql import select

sys.path.insert(0, '/app/mysqldb')
from models import rsid, transcript, gene, RSID, Transcript, Gene

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_gene_chrom(conn, gene_id):
    """return chromsome which the gene is on
    """
    query = select([gene]).where(gene.c.id == gene_id)
    results = conn.execute(query)
    result = next(results, None)
    if not result:
        return None
    return Gene(*result).chrom


def get_gene_tx(conn, gene_id):
    """return all transcripts belong to the gene
    """
    query = select([transcript]).where(transcript.c.gene_id == gene_id)
    results = conn.execute(query)
    return list(starmap(Transcript, results))


def get_refseq_tx(conn, refseq):
    """return all transcripts belong to the refseq
    """
    query = select([transcript]).where(transcript.c.name == refseq)
    results = conn.execute(query)
    return list(starmap(Transcript, results))


def get_rsid_chromvar(conn, name):
    """return all chromosome variant of the rsid
    """
    if not name[2:].isdigit():
        return []
    query = select([rsid]).where(rsid.c.name == int(name[2:]))
    results = conn.execute(query)
    return list(starmap(RSID, results))
