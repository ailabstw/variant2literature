"""normalize variant to vcfs
"""
import os
import json
import re

from .hgvs import HGVS
from .utils import Position, RNAVariant, ProteinVariant, MTVariant


class VarNormalizer:
    """convert variants to vcf
    """

    def __init__(self, load_all_genome=False):
        host = os.environ['MYSQL_HOST']
        port = os.environ['MYSQL_PORT']
        passwd = os.environ['MYSQL_ROOT_PASSWORD']
        db_uri = f'mysql+pymysql://root:{passwd}@{host}:{port}/gene'

        self.hgvs = HGVS(db_uri=db_uri,
                         ref_genome_path='/app/models/ucsc.hg19.fasta',
                         load_all_genome=load_all_genome)

    def connect(self):
        """return a Connection object
        """
        return self.hgvs.connect()

    def parse_position(self, pos_str, tx):
        """parse position string
        """
        if '?' in pos_str:
            return None

        try:
            offset, utr3p = 0, False
            pattern = '^(?P<pos>([+*-]?[0-9]+)|((IVS|Intron|ivs|intron)(-?)[0-9]+))(?P<offset>([+-][0-9]+)?)'
            m = re.search(pattern, pos_str)
            d = m.groupdict()
            pos = d['pos']
            if d['offset']:
                offset = int(d['offset'])

            if pos.upper().startswith('I'):
                if pos.upper().startswith('IVS'):
                    k = int(pos[3:])
                elif pos.upper().startswith('INTRON'):
                    k = int(pos[6:])

                if not tx:
                    return None
                position = self.hgvs.get_intron_pos(k, offset, tx)
                if not position:
                    return None
                pos = position.pos
                utr3p = position.utr3p

            elif pos[0] == '*':
                utr3p = True
                pos = pos[1:]
            pos = int(pos)

        except Exception:
            # print('parse fail: {}'.format(pos_str))
            # raise
            return None

        return Position(pos, offset, utr3p)

    def fix_variant(self, tx, mut_type, start, end, ref, alt):  # pylint: disable=too-many-arguments
        """fix the wrong representation of the variant if possible
        """
        if tx and start == end and start.offset < 0:
            if not self.hgvs.is_neg_offset_intron(tx, start):
                start = Position(start.pos, 0, False)
                end = Position(-start.offset, 0, False)

        if mut_type in {'SUB', 'DEL'} and ref:
            if ref.isdigit():
                ref_len = int(ref)
                ref = ''
            else:
                ref_len = len(ref)

            if tx:
                end = self.hgvs.increase_pos(start, ref_len - 1, tx)
            else:
                end = start + ref_len - 1

        if mut_type == 'INS':
            if tx:
                end = self.hgvs.increase_pos(start, 1, tx)
            else:
                end = start + 1

        if end is None:
            return None

        return start, end, ref, alt

    def parse_rna_var(self, tx, var, fix=True):
        """parse rna variant
        """
        start = self.parse_position(var['start'], tx)
        end = self.parse_position(var['end'], tx)
        mut_type = var['mut_type']
        ref, alt = var.get('wild_type', ''), var.get('mutant', '')

        if not start:
            return None

        if fix:
            fixed = self.fix_variant(tx, mut_type, start, end, ref, alt)
            if not fixed:
                return None
            start, end, ref, alt = fixed

        if re.search(r'[^ACGT]', ref):
            ref = ''
        if re.search(r'[^ACGT]', alt):
            return None
        if mut_type == 'INS' and not alt:
            return None

        dup = var.get('dup', '')
        dup = int(dup) if dup.isdigit() else 2

        return RNAVariant(mut_type, start, end, ref, alt, dup)

    def to_vcf_protein(self, gene_ids, var_dict):
        """convert protein var to vcf
        """
        if not var_dict['start'].isdigit():
            return

        start = int(var_dict['start'])
        ref, alt = var_dict['wild_type'], var_dict['mutant']
        protein_var = ProteinVariant(start, ref, alt)
        for gene_id in gene_ids:
            for tx in self.hgvs.get_gene_tx(gene_id):
                for vcf in self.hgvs.protein_to_chrom(tx, protein_var):
                    yield (tx.name, vcf)

    def to_vcf_mt(self, var_dict, fix=True):
        """convert mt var to vcf
        """
        start = self.parse_position(var_dict['start'], tx=False).pos
        end = self.parse_position(var_dict['end'], tx=False).pos
        mut_type = var_dict['mut_type']
        ref, alt = var_dict.get('wild_type', ''), var_dict.get('mutant', '')

        if fix:
            fixed = self.fix_variant(None, mut_type, start, end, ref, alt)
            if not fixed:
                return
            start, end, ref, alt = fixed

        if re.search(r'[^ACGT]', ref):
            ref = ''
        if re.search(r'[^ACGT]', alt):
            return
        if mut_type == 'INS' and not alt:
            return

        dup = var_dict.get('dup', '')
        dup = int(dup) if dup.isdigit() else 2

        mt_var = MTVariant(mut_type, start, end, ref, alt, dup)

        for vcf in self.hgvs.mt_to_chrom(mt_var):
            yield ('Mitochondrial', vcf)


    def to_vcf_rna(self, gene_ids, var_dict, fix=True):
        """convert rna var to vcf
        """
        for gene_id in gene_ids:
            chrom = self.hgvs.get_gene_chrom(gene_id)
            if chrom == 'MT':
                yield from self.to_vcf_mt(var_dict, fix)
                return

            for tx in self.hgvs.get_gene_tx(gene_id):
                var = self.parse_rna_var(tx, var_dict, fix)
                if not var:
                    continue
                for vcf in self.hgvs.rna_to_chrom(tx, var):
                    yield (tx.name, vcf)

    def to_vcf_rsid(self, rsid):
        """convert rsid to vcf
        """
        for vcf in self.hgvs.rsid_to_chrom(rsid):
            yield ('', vcf)

    def to_vcf(self, gene_ids, var_json, fix=True):
        """convert variant json to vcf
        """
        var = json.loads(var_json)


        if var['mut_type'] == 'VCF':
            chrom = var['chrom']
            position = int(var['position'])
            ref = var['ref']
            alt = var['alt']
            yield ('', (chrom, position, ref, alt))

        elif var['mut_type'].startswith('RSID+'):
            var['mut_type'] = var['mut_type'][5:]
            rsid_vcfs = set(self.to_vcf_rsid(var['rsid']))
            if 'seq_types' in var and 'p' in var['seq_types'] and var['mut_type'] == 'SUB':
                hgvs_vcfs = set(self.to_vcf_protein(gene_ids, var))

            elif 'seq_types' in var and {'c', 'n', 'r'} & set(var['seq_types']):
                hgvs_vcfs = set(self.to_vcf_rna(gene_ids, var, fix))

            elif 'seq_types' in var and {'m'} & set(var['seq_types']):
                hgvs_vcfs = set(self.to_vcf_mt(var))
            else:
                return

            vcfs = rsid_vcfs & hgvs_vcfs
            yield from vcfs

        elif var['mut_type'] == 'RSID':
            yield from self.to_vcf_rsid(var['rsid'])

        # for protein, only supporting substitution
        elif 'seq_types' in var and 'p' in var['seq_types'] and var['mut_type'] == 'SUB':
            yield from self.to_vcf_protein(gene_ids, var)

        elif 'seq_types' in var and {'c', 'n', 'r'} & set(var['seq_types']):
            yield from self.to_vcf_rna(gene_ids, var, fix)

        elif 'seq_types' in var and {'m'} & set(var['seq_types']):
            yield from self.to_vcf_mt(var)
