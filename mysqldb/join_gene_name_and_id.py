"""collect transcripts from refseq, ensembl, ucsc
"""
import re
import gzip
from collections import defaultdict


def clean_text(name):
    name = re.sub(r'[^A-Za-z0-9]', '', name)
    return name.upper()


def get_transcripts():
    symbol_id = defaultdict(set)

    with gzip.open('/app/models/Homo_sapiens.gene_info.gz') as f:
        for line in f:
            line = line.decode('ascii')
            tax_id, gene_id, symbol, _, synonyms, _ = line.split('\t', 5)

            symbol_id[clean_text(symbol)].add(gene_id)
            symbol_id['LOC{}'.format(gene_id)].add(gene_id)

            if synonyms != '-':
                for synonym in synonyms.split('|'):
                    symbol_id[clean_text(synonym)].add(gene_id)

    tx_gene = dict()

    with gzip.open('/app/models/gene2refseq.gz') as f:
        for line in f:
            line = line.decode('ascii')
            if not line.startswith('9606\t'):
                continue
            tax_id, gene_id, _, rna_refseq, _ = line.split('\t', 4)
            rna_refseq = rna_refseq.split('.')[0]
            tx_gene[rna_refseq] = gene_id

    with gzip.open('/app/models/ensemblToGeneName.txt.gz') as f:
        for line in f:
            line = line.decode('ascii')
            ensid, gene_name = line.strip('\n').split(None, 1)
            for gene_id in symbol_id.get(clean_text(gene_name), []):
                tx_gene[ensid] = gene_id

    with gzip.open('/app/models/kgAlias.txt.gz') as f:
        for line in f:
            line = line.decode('ascii')
            kg_id, gene_name = line.strip('\n').split(None, 1)
            for gene_id in symbol_id.get(clean_text(gene_name), []):
                tx_gene[kg_id] = gene_id

    tx_files = [
        '/app/models/ncbiRefSeq.txt.gz',
        '/app/models/ensGene.txt.gz',
        '/app/models/knownGene.txt.gz'
    ]

    data = []
    for tx_file in tx_files:
        with gzip.open(tx_file) as f:
            for line in f:
                line = line.decode('ascii')
                (_, name, chrom, strand, tx_start, tx_end,
                 cds_start, cds_end, _, exon_starts, exon_ends, _) = line.split('\t', 11)

                name = name.split('.')[0]
                gene_id = tx_gene.get(name, None)
                if not gene_id:
                    continue

                data.append({
                    'name': name,
                    'chrom': chrom,
                    'strand': strand,
                    'tx_start': tx_start,
                    'tx_end': tx_end,
                    'cds_start': cds_start,
                    'cds_end': cds_end,
                    'exon_starts': exon_starts,
                    'exon_ends': exon_ends,
                    'gene_id': gene_id,
                })
    return data
