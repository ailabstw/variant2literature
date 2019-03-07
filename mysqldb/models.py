# pylint: disable=invalid-name
"""models for tables in mysql
"""
import gzip
import os
import csv
import logging
import time

from sqlalchemy.pool import NullPool
from sqlalchemy import create_engine, MetaData, Table, Column, Index
from sqlalchemy.types import Integer, BigInteger, String, Text
from sqlalchemy.schema import CreateTable
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.dialects.mysql import LONGBLOB, MEDIUMTEXT

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

metadata = MetaData()

paper_status_cnt = Table(
    'paper_status_cnt', metadata,
    Column('status', Integer(), primary_key=True),
    Column('cnt', Integer(), nullable=False, default=0),
    mysql_engine='MyISAM',
)

paper_status = Table(
    'paper_status', metadata,
    Column('_id', String(50), primary_key=True),
    Column('pmid', String(50), nullable=False),
    Column('pmcid', String(50), nullable=False),
    Column('path', String(100), nullable=False),
    Column('status', Integer(), nullable=False, default=0),
    mysql_engine='MyISAM',
)

Index('idx_pmid', paper_status.c.pmid)
Index('idx_pmcid', paper_status.c.pmcid)
Index('idx_status', paper_status.c.status)

var_pmid = Table(
    'var_pmid', metadata,
    Column('chrom', String(50), primary_key=True),
    Column('position', String(50), primary_key=True),
    Column('ref', String(200), primary_key=True),
    Column('alt', String(200), primary_key=True),
    Column('gene_id', Integer(), nullable=False),
    Column('transcript', String(50), nullable=False),
    Column('_id', String(50), primary_key=True),
    Column('file_idx', Integer(), primary_key=True),
    Column('in_table', Integer(), primary_key=True),
    Column('table_idx', Integer(), primary_key=True),
    Column('row', Integer(), primary_key=True),
    Column('col', Integer(), primary_key=True),
    Column('start', Integer(), primary_key=True),
    Column('end', Integer(), nullable=False),
    Column('mention', String(200), nullable=False),
    Column('gene_start', Integer(), nullable=False),
    Column('gene_end', Integer(), nullable=False),
    mysql_engine='MyISAM',
)

Index('idx_paper_id', var_pmid.c._id)
Index('idx_gene', var_pmid.c.gene_id)
Index('idx_var', var_pmid.c.chrom, var_pmid.c.position,
      var_pmid.c.ref, var_pmid.c.alt)
Index('idx_pos', var_pmid.c.chrom, var_pmid.c.position)

rsid = Table(
    'rsid', metadata,
    Column('_id', Integer(), primary_key=True, autoincrement=True),
    Column('name', BigInteger(), nullable=False),
    Column('chrom', String(50), nullable=False),
    Column('start', BigInteger(), nullable=False),
    Column('end', BigInteger(), nullable=False),
    Column('ref', Text(), nullable=False),
    Column('observed', Text(), nullable=False),
    mysql_engine='MyISAM',
)

transcript = Table(
    'transcript', metadata,
    Column('name', String(50), primary_key=True),
    Column('chrom', String(50), primary_key=True),
    Column('strand', String(2), nullable=False),
    Column('tx_start', BigInteger(), nullable=False),
    Column('tx_end', BigInteger(), nullable=False),
    Column('cds_start', BigInteger(), nullable=False),
    Column('cds_end', BigInteger(), nullable=False),
    Column('exon_starts', Text(), nullable=False),
    Column('exon_ends', Text(), nullable=False),
    Column('gene_id', Integer()),
    mysql_engine='MyISAM',
)

Index('idx_gene', transcript.c.gene_id)


gene = Table(
    'gene', metadata,
    Column('id', Integer(), primary_key=True),
    Column('symbol', String(50), nullable=False),
    Column('chrom', String(50), nullable=False),
    mysql_engine='MyISAM',
)

Index('idx_symbol', gene.c.symbol)


pubmed = Table(
    'pubmed', metadata,
    Column('id', Integer(), primary_key=True, autoincrement=True),
    Column('pmid', String(50), nullable=False),
    Column('year', Integer(), nullable=False),
    Column('authors', Text(), nullable=False),
    Column('journal', String(200), nullable=False),
    Column('title', Text(), nullable=False),
    mysql_charset='utf8mb4',
    mysql_engine='MyISAM',
)

index = Index('idx_pmid', pubmed.c.pmid)


class PaperStatus:
    """paper status
    """
    def __init__(self, _id, pmid, pmcid, path, status):
        self._id = _id
        self.pmid = pmid
        self.pmcid = pmcid
        self.path = path
        self.status = status


class PaperStatusCnt:
    """paper status cnt
    """
    def __init__(self, status, cnt):
        self.status = status
        self.cnt = cnt


class VarPmid:
    """variant to paper
    """
    def __init__(self, chrom, position, ref, alt,  # pylint: disable=too-many-arguments
                 gene_id, transcript, pmid, file_idx,
                 in_table, table_idx, row, col,
                 start, end, mention, gene_start, gene_end):

        self.chrom = chrom
        self.position = position
        self.ref = ref
        self.alt = alt
        self.gene_id = gene_id
        self.transcript = transcript
        self.pmid = pmid
        self.file_idx = file_idx
        self.in_table = in_table
        self.table_idx = table_idx
        self.row = row
        self.col = col
        self.start = start
        self.end = end
        self.mention = mention
        self.gene_start = gene_start
        self.gene_end = gene_end


class RSID:
    """rsid
    """
    def __init__(self, _id, name, chrom, start, end, ref, observed):  # pylint: disable=too-many-arguments
        self._id = _id
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ref = ref
        self.observed = observed


class Transcript:  # pylint: disable=too-many-instance-attributes
    """transcript
    """
    def __init__(self, name, chrom, strand, tx_start, tx_end, # pylint: disable=too-many-arguments
                 cds_start, cds_end, exon_starts, exon_ends, gene_id):
        self.name = name
        self.chrom = chrom
        self.strand = 1 if strand == '+' else -1
        self.tx_start = tx_start
        self.tx_end = tx_end
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.exon_starts = tuple(map(int, exon_starts.strip(',').split(',')))
        self.exon_ends = tuple(map(int, exon_ends.strip(',').split(',')))
        self.gene_id = gene_id
        self.is_coding = (cds_start != cds_end)


class Gene:
    """gene
    """
    def __init__(self, _id, symbol, chrom):
        self._id = _id
        self.symbol = symbol
        self.chrom = chrom


class Pubmed:
    """pubmed
    """
    def __init__(self, _id, pmid, year, authors, journal, title):
        self._id = _id
        self.pmid = pmid
        self.year = year
        self.authors = authors
        self.journal = journal
        self.title = title


def insert_rsid(engine):  # pylint: disable=too-many-locals
    """insert rsid into database
    """

    query = ('INSERT IGNORE INTO `rsid` '
             '(`name`, `chrom`, `start`, `end`, `ref`, `observed`) '
             'VALUES (%s, %s, %s, %s, %s, %s)')

    t = time.time()
    with gzip.open('/app/models/snp150.txt.gz') as f, engine.connect() as conn:
        conn.execute('truncate rsid')
        values, i = [], -1
        for i, line in enumerate(f):
            row = line.decode('ascii').split('\t', 10)
            name = row[4][2:]
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            ref = row[7]
            observed = row[9]

            values.append((name, chrom, start, end, ref, observed))

            if (i + 1) % 100000 == 0:
                conn.execute(query, values)
                print(i + 1, 'rsid inserted', time.time() - t)
                values = []

        conn.execute(query, values)
        print(i + 1, 'rsid inserted', time.time() - t)



def insert_transcripts(engine):
    """insert transcripts into database
    """
    print('loading transcripts ...')
    from join_gene_name_and_id import get_transcripts
    transcripts = get_transcripts()
    print('transcripts load OK')

    t = time.time()
    with engine.connect() as conn:
        conn.execute('truncate transcript')
        values, i = [], -1
        for i, row in enumerate(transcripts):
            values.append(row)

            if (i + 1) % 100000 == 0:
                conn.execute(transcript.insert(), values)  # pylint: disable=no-value-for-parameter
                print(i + 1, 'transcripts inserted', time.time() - t)
                values = []

        conn.execute(transcript.insert(), values)  # pylint: disable=no-value-for-parameter
        print(i + 1, 'transcripts inserted', time.time() - t)


def insert_genes(engine):
    """insert transcripts into database
    """

    t = time.time()
    with gzip.open('/app/models/Homo_sapiens.gene_info.gz') as f, engine.connect() as conn:
        conn.execute('truncate gene')
        f.readline()

        values, i = [], -1
        for i, line in enumerate(f):
            row = line.decode('ascii').split('\t')
            values.append({
                'id': row[1],
                'symbol': row[2],
                'chrom': row[6],
            })

            if (i + 1) % 10000 == 0:
                conn.execute(gene.insert(), values)  # pylint: disable=no-value-for-parameter
                print(i + 1, 'genes inserted', time.time() - t)
                values = []

        conn.execute(gene.insert(), values)  # pylint: disable=no-value-for-parameter
        print(i + 1, 'genes inserted', time.time() - t)


host = os.environ['MYSQL_HOST']
port = os.environ['MYSQL_PORT']
passwd = os.environ['MYSQL_ROOT_PASSWORD']
engine = create_engine(f'mysql+pymysql://root:{passwd}@{host}:{port}/gene',
                       poolclass=NullPool)


def create_db():
    _engine = create_engine(f'mysql+pymysql://root:{passwd}@{host}:{port}',
                            pool_pre_ping=True)
    with _engine.connect() as conn:
        conn.execute('create database IF NOT EXISTS gene')


def init():
    """init database
    """
    create_db()
    metadata.create_all(engine)

    insert_genes(engine)
    insert_transcripts(engine)
    insert_rsid(engine)

    index = Index('idx_name', rsid.c.name)
    index.create(bind=engine)


if __name__ == '__main__':
    init()
