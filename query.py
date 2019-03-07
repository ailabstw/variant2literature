import sys
import re
import json
import itertools
from collections import defaultdict

from sqlalchemy.sql import and_

from var_utils import VarNormalizer
sys.path.insert(0, '/app/mysqldb')
from models import (engine,  # pylint: disable=no-name-in-module, wrong-import-position
                    gene, var_pmid, VarPmid)

ALL_TO_ONE = {
#{{{
    'ALANINE': 'A',
    'ALA': 'A',
    'ARGININE': 'R',
    'ARG': 'R',
    'ASPARAGINE': 'N',
    'ASN': 'N',
    'ASPARTATE': 'D',
    'ASPARTICACID': 'D',
    'ASP': 'D',
    'ASX': 'B',
    'CYSTEINE': 'C',
    'CYS': 'C',
    'GLUTAMATE': 'E',
    'GLUTAMICACID': 'E',
    'GLU': 'E',
    'GLUTAMINE': 'Q',
    'GLN': 'Q',
    'GLX': 'Z',
    'GLYCINE': 'G',
    'GLY': 'G',
    'HISTIDINE': 'H',
    'HIS': 'H',
    'ISOLEUCINE': 'I',
    'ILE': 'I',
    'LEUCINE': 'L',
    'LEU': 'L',
    'LYSINE': 'K',
    'LYS': 'K',
    'METHIONINE': 'M',
    'MET': 'M',
    'PHENYLALANINE': 'F',
    'PHE': 'F',
    'PROLINE': 'P',
    'PRO': 'P',
    'SERINE': 'S',
    'SER': 'S',
    'THREONINE': 'T',
    'THYMINE': 'T',
    'THR': 'T',
    'TRYPTOPHAN': 'W',
    'TRP': 'W',
    'TYROSINE': 'Y',
    'TYR': 'Y',
    'VALINE': 'V',
    'VAL': 'V',
    'STOP': 'X',
    'TER': 'X',
    '*': 'X',
}#}}}

STOP = r'Stop|STOP|stop'
P1CODE = r'A|R|N|D|B|C|E|Q|Z|G|H|I|L|K|M|F|P|S|T|W|Y|V|X|\*'
P3CODE = r'Ala|Arg|Asn|Asp|Asx|Cys|Glu|Gln|Glx|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter'
IVS = r'(IVS|ivs|Intron|intron)\d+([-+]\d+)?'
C_POS1 = r'([-*]?\d+)([+-]\d+)?'
C_POS2 = r'((' + IVS + ')|(' + C_POS1 + '))'
SUB = r'[>/]'
P0 = r'(' + P1CODE + r'|' + P3CODE + r'|' + STOP + r')'

P1 = r'(?P<wild_type>(' + P1CODE + r'|' + P3CODE + r'))'
P2 = r'(?P<mutant>(' + P1CODE + r'|' + P3CODE + r'))'
C_POS = r'(?P<start>(' + C_POS2 + r'))(_(?P<end>' + C_POS2 + r'))?'
P_RS = r'(?P<seq>p)'
C_RS = r'(?P<seq>[cmrg])'

patterns = [
    ('SUB', 'p', r'^(' + P_RS + r'\.)?' + P1 + r'(?P<start>\d+)' + P2 + r'$'),  # p.M34T
    ('SUB', 'c', r'^(' + C_RS + r'\.)?' + r'(?P<wild_type>[ACGT])' + r'(?P<start>\d+)' + r'(?P<mutant>[ACGT])' + r'$'),  # c.A34T
    ('SUB', 'c', r'^(' + C_RS + r'\.)?' + C_POS + r'(?P<wild_type>[ACGT]+)' + SUB + r'(?P<mutant>[ACGT]+)$'),  # c.123C>A
    ('SUB', 'p', r'^(' + P_RS + r'\.)?' + r'(?P<start>\d+)' + P1 + SUB + r'' + P2 + '[\])]?;?$'),  # p.123C>A
    ('DEL', 'c', r'^(' + C_RS + r'\.)?' + C_POS + r'del(?P<wild_type>[ACGT\d]*)(bp)?[\])]?;?$'),  # c.35delG
    ('DUP', 'c', r'^(' + C_RS + r'\.)?' + C_POS + r'dup(?P<wild_type>[ACGT\d]*)(bp)?[\])]?;?$'),  # c.35dupG
    ('INS', 'c', r'^(' + C_RS + r'\.)?' + C_POS + r'ins(?P<mutant>[ACGT]*)(bp)?[\])]?;?$'),  # c.35insG
    ('IND', 'c', r'^(' + C_RS + r'\.)?' + C_POS + r'del(?P<wild_type>[ACGT\d]*)ins(?P<mutant>[ACGT]*)$'),  # c.35delinsG
    ('RSID', '', r'(?P<rsid>(^rs\d+$))'),
]


def error(msg):
    print(msg)
    exit(0)

fout = sys.stdout
if len(sys.argv) > 1:
    fout = open(sys.argv[1], 'w')

gene_name = input('Gene: ')
variant_name = input('Variant: ')
print()

if gene_name.isdigit():
    gene_id = int(gene_name)

else:
    with engine.connect() as conn:
        query = gene.select().where(gene.c.symbol == gene_name)
        rows = conn.execute(query)
        row = next(rows, None)
        if not row:
            error('gene symbol not found')
    gene_id = row[0]

print('GeneID:', gene_id, file=fout)


for j, (mutation_type, rs, pattern) in enumerate(patterns):
    m = re.match(pattern, variant_name)
    if m:
        d = m.groupdict()
        break
else:
    error('cannot reconigze the variant or variant type not supported')

d['mut_type'] = mutation_type
if 'rsid' not in d:
    if not d.get('end', None):
        d['end'] = d['start']

    if 'seq' in d and d['seq']:
        seq = d['seq']
    else:
        seq = rs
    del d['seq']
    d['seq_types'] = [seq]

print('VariantJson:', json.dumps(d), file=fout)
print(file=fout)


var_normalizer = VarNormalizer()
vcfs = []
with var_normalizer.connect():
    for tx, vcf in var_normalizer.to_vcf([gene_id], json.dumps(d)):
        vcfs.append((tx, vcf))

results = []
with engine.connect() as conn:
    for tx, (chrom, pos, ref, alt) in vcfs:
        cond = and_(
            var_pmid.c.chrom == chrom,
            var_pmid.c.position == pos,
            var_pmid.c.ref == ref,
            var_pmid.c.alt == alt,
        )
        query = var_pmid.select().where(cond)
        rows = conn.execute(query)
        results += list(itertools.starmap(VarPmid, rows))

d = defaultdict(dict)
for r in results:
    d[r.pmid][(r.file_idx, r.table_idx, r.row, r.col, r.start)] = r.mention

print(f'{len(d)} results found:', file=fout)
for pmid, pmid_dict in d.items():
    print('pmid:', pmid, file=fout)
    for mention in pmid_dict.values():
        print('\tvariant:', mention, file=fout)

if len(sys.argv) > 1:
    fout.close()
