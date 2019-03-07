"""hgvs
"""
import logging
from sqlalchemy import create_engine
from .seqdb import SequenceFileDB
from .utils import (rna_to_protein, protein_to_rna, three_to_one,
                    rev_p1, revcomp, Position, ChromVariant)
from .gene_db import get_gene_tx, get_refseq_tx, get_rsid_chromvar, get_gene_chrom

logger = logging.getLogger(__name__)


class VCF:
    """normalize chromsome variant to vcf format
    """
    def __init__(self, genome, chrom_var):
        self.genome = genome
        self.chrom = chrom_var.chrom
        self.pos = chrom_var.pos
        self.ref = list(chrom_var.ref)
        self.alt = list(chrom_var.alt)
        self.done = False
        self._normalize()

    def __iter__(self):
        ref = ''.join(self.ref)
        alt = ''.join(self.alt)
        yield from [self.chrom, self.pos, ref, alt]

    def as_tuple(self):
        """return vcf as tuple
        """
        return tuple(self.__iter__())

    def _ensure_non_empty(self):
        if not self.ref or not self.alt:
            if self.pos > 1:
                self.pos -= 1
                res = self.genome[self.chrom][self.pos - 1].upper()
                self.ref = [res] + self.ref
                self.alt = [res] + self.alt
            else:
                res = self.genome[self.chrom][self.pos + len(self.ref) - 1].upper()
                self.ref = self.ref + [res]
                self.alt = self.alt + [res]
            self.done = False

    def _right_trim(self):
        if len(self.ref) == 1 and len(self.alt) == 1:
            return
        if all(x == 'N' for x in self.ref + self.alt):
            return

        while (self.ref and self.alt and self.ref[-1] == self.alt[-1] and
               ((len(self.ref) > 1 and len(self.alt) > 1) or self.pos > 1)):
            self.ref.pop()
            self.alt.pop()
            self.done = False

    def _left_trim(self):
        while self.ref[0] == self.alt[0] and min(len(self.ref), len(self.alt)) > 1:
            self.pos += 1
            self.ref.pop(0)
            self.alt.pop(0)

    def _normalize(self):
        if self.pos <= 0:
            return
        if self.pos + len(self.ref) - 1 >= len(self.genome[self.chrom]):
            return
        self._ensure_non_empty()
        while not self.done:
            self.done = True
            self._right_trim()
            self._ensure_non_empty()
        self._left_trim()


class HGVS:
    """convert hgvs names and some utils
    """
    def __init__(self, db_uri, ref_genome_path,
                 load_all_genome=False):
        self.engine = create_engine(db_uri, pool_pre_ping=True)
        self.load_all_genome = load_all_genome
        self.genome = SequenceFileDB(ref_genome_path, load_all=load_all_genome)
        self.conn = self.engine.connect()
        self.conn.close()

    def connect(self):
        """return a Connection object
        """
        self.conn = self.engine.connect()
        return self.conn

    def get_utr_len(self, tx, pos, exon_starts, exon_ends):
        """return length of the untranslated region
        """
        utr = 0
        for exon_start, exon_end in zip(exon_starts, exon_ends):
            if tx.strand * exon_start <= tx.strand * pos < tx.strand * exon_end:
                utr += tx.strand * (pos - exon_start)
                break
            else:
                utr += tx.strand * (exon_end - exon_start)
        return utr

    def is_neg_offset_intron(self, tx, position):
        """return if the position is a valid intron position in the transcript
        """
        chrom_pos = self.rna_to_chrom_pos(tx, Position(position.pos, -1, position.utr3p))
        if not chrom_pos:
            return False
        exon_starts = rev_p1(tx.exon_ends) if tx.strand == -1 else tx.exon_starts

        if chrom_pos not in exon_starts:
            return False
        idx = exon_starts.index(chrom_pos)

        if idx > 0:
            half_intron_len = abs(exon_starts[idx] - exon_starts[idx - 1] + tx.strand) // 2
            if -position.offset > half_intron_len:
                return False
        return True

    def rna_to_chrom_pos(self, tx, position):  # pylint: disable=too-many-locals
        """convert the RNA position to the chromosome position
        """
        strand = tx.strand
        if tx.strand == -1:
            exon_starts, exon_ends = rev_p1(tx.exon_ends), rev_p1(tx.exon_starts)
            cds_start, cds_end = rev_p1([tx.cds_start, tx.cds_end])
            tx_start, tx_end = rev_p1([tx.tx_start, tx.tx_end])
        else:
            exon_starts, exon_ends = tx.exon_starts, tx.exon_ends
            cds_start, cds_end = tx.cds_start, tx.cds_end
            tx_start, tx_end = tx.tx_start, tx.tx_end

        if tx.is_coding:
            utr5p_len = self.get_utr_len(tx, cds_start, exon_starts, exon_ends)
            utr3p_len = self.get_utr_len(tx, cds_end, exon_starts, exon_ends)
        else:
            utr5p_len = 0
            utr3p_len = strand * (tx_end - tx_start)

        pos = position.pos
        chrom_pos = tx_start
        if position.utr3p:
            pos += utr3p_len
        elif pos < 0:
            pos += utr5p_len + 1
        else:
            pos += utr5p_len

        if pos < 0:
            return chrom_pos + (pos + position.offset) * strand

        j = 0
        while pos > 0:
            step = min(pos, max(0, strand * (exon_ends[j] + strand - chrom_pos)))
            chrom_pos += step * strand
            pos -= step
            while strand * chrom_pos >= strand * exon_ends[j] + 1:
                j += 1
                if j == len(exon_ends):
                    if not position.utr3p:
                        return None
                    return chrom_pos + (pos + position.offset) * strand
            chrom_pos = strand * max(strand * chrom_pos, strand * exon_starts[j] + 1)

        return chrom_pos + position.offset * strand

    def mt_to_chrom(self, mt_var):
        """convert a mt variant to a chromosome variant
        """
        start = mt_var.start
        end = mt_var.end if mt_var.end else mt_var.start
        ref, alt = mt_var.ref, mt_var.alt

        if mt_var.mut_type == 'INS':
            start = end

        start_chrom_pos = start
        end_chrom_pos = end

        ref_seq = self.genome['chrM'][start_chrom_pos - 1:end_chrom_pos].upper()

        if ref:
            if ref_seq != ref:
                logger.debug('%s != %s at %s', ref_seq, ref, start_chrom_pos)
                return

        if mt_var.mut_type in {'DEL', 'INDEL'}:
            ref = ref_seq

        if mt_var.mut_type in {'DUP'}:
            ref = ''
            alt = ref_seq * (mt_var.dup - 1)

        yield VCF(self.genome, ChromVariant('chrM', start_chrom_pos, ref, alt)).as_tuple()

    def rna_to_chrom(self, tx, rna_var):
        """convert a rna variant to a chromosome variant
        """
        start = rna_var.start
        end = rna_var.end if rna_var.end else rna_var.start
        ref, alt = rna_var.ref, rna_var.alt
        if tx.strand == -1:
            ref, alt = revcomp(ref), revcomp(alt)
            start, end = end, start

        if rna_var.mut_type == 'INS':
            start = end

        start_chrom_pos = self.rna_to_chrom_pos(tx, start)
        end_chrom_pos = self.rna_to_chrom_pos(tx, end)
        if not start_chrom_pos or not end_chrom_pos:
            return

        ref_seq = self.genome[tx.chrom][start_chrom_pos - 1:end_chrom_pos].upper()

        if ref:
            if ref_seq != ref:
                logger.debug('%s != %s at %s', ref_seq, ref, start_chrom_pos)
                return

        if rna_var.mut_type in {'DEL', 'INDEL'}:
            ref = ref_seq

        if rna_var.mut_type in {'DUP'}:
            ref = ''
            alt = ref_seq * (rna_var.dup - 1)

        yield VCF(self.genome, ChromVariant(tx.chrom, start_chrom_pos, ref, alt)).as_tuple()

    def protein_to_chrom(self, tx, protein_var):
        """convert a protein variant to a chromosome variant
        """
        ref, alt = protein_var.ref.upper(), protein_var.alt.upper()
        alt = ref if alt == '=' else alt
        ref, alt = three_to_one(ref), three_to_one(alt)
        pos = protein_var.pos

        if tx.strand == -1:
            chrom_pos = self.rna_to_chrom_pos(tx, Position(pos * 3))
        else:
            chrom_pos = self.rna_to_chrom_pos(tx, Position(pos * 3 - 2))

        if not chrom_pos:
            return

        ref_seq = self.genome[tx.chrom][chrom_pos - 1:chrom_pos + 2].upper()

        if len(ref_seq) < 3:
            logger.debug('out of range')
            return

        if tx.strand == -1:
            ref_seq = revcomp(ref_seq)

        if rna_to_protein(ref_seq) != ref:
            logger.debug('%s != %s', rna_to_protein(ref_seq), ref)
            return

        if tx.strand == -1:
            ref_seq = revcomp(ref_seq)

        for alt_seq in protein_to_rna(alt):
            if tx.strand == -1:
                alt_seq = revcomp(alt_seq)
            yield VCF(self.genome, ChromVariant(tx.chrom, chrom_pos, ref_seq, alt_seq)).as_tuple()

    def refseq_rna_to_chrom(self, refseq, rna_var):
        """convert a rna variant to chromosome variant according to the refseq
        """
        for tx in get_refseq_tx(self.conn, refseq):
            logger.debug(tx.name)
            yield from self.rna_to_chrom(tx, rna_var)

    def gene_rna_to_chrom(self, gene_id, rna_var):
        """convert a rna variant to chromosome variant according to the gene
        """
        for tx in get_gene_tx(self.conn, gene_id):
            logger.debug(tx.name)
            yield from self.rna_to_chrom(tx, rna_var)

    def refseq_protein_to_chrom(self, refseq, protein_var):
        """convert a protein variant to chromosome variant according to the refseq
        """
        for tx in get_refseq_tx(self.conn, refseq):
            logger.debug(tx.name)
            yield from self.protein_to_chrom(tx, protein_var)

    def gene_protein_to_chrom(self, gene_id, protein_var):
        """convert a protein variant to chromosome variant according to the gene
        """
        for tx in get_gene_tx(self.conn, gene_id):
            logger.debug(tx.name)
            yield from self.protein_to_chrom(tx, protein_var)

    def exon_distance(self, strand, pos, start, end):
        """return the distance to the nearest exon
        """
        if strand * start < strand * pos <= strand * end:
            return 0
        return min(abs(pos - (start + strand)), abs(pos - end))

    def chrom_to_rna_pos(self, tx, chrom_pos):  # pylint: disable=too-many-locals
        """convert the chromosome position to the rna position given the transcript
        """
        strand = tx.strand
        if tx.strand == -1:
            exon_starts, exon_ends = rev_p1(tx.exon_ends), rev_p1(tx.exon_starts)
            cds_start, cds_end = rev_p1([tx.cds_start, tx.cds_end])
            tx_start, tx_end = rev_p1([tx.tx_start, tx.tx_end])
        else:
            exon_starts, exon_ends = tx.exon_starts, tx.exon_ends
            cds_start, cds_end = tx.cds_start, tx.cds_end
            tx_start, tx_end = tx.tx_start, tx.tx_end

        if tx.is_coding:
            utr5p_len = self.get_utr_len(tx, cds_start, exon_starts, exon_ends)
            utr3p_len = self.get_utr_len(tx, cds_end, exon_starts, exon_ends)
        else:
            utr5p_len = 0
            utr3p_len = strand * (tx_end - tx_start)

        min_exon_distance = min(map(lambda exon: self.exon_distance(strand, chrom_pos, *exon),
                                    zip(exon_starts, exon_ends)))

        pos, offset = 0, 0
        for exon_start, exon_end in zip(exon_starts, exon_ends):
            if self.exon_distance(strand, chrom_pos, exon_start, exon_end) == min_exon_distance:
                if strand * chrom_pos > strand * exon_end:
                    pos += strand * exon_end - strand * exon_start
                    offset = min_exon_distance
                elif strand * chrom_pos <= strand * exon_start:
                    pos += 1
                    offset = -min_exon_distance
                else:
                    pos += strand * chrom_pos - strand * exon_start
                break
            else:
                pos += strand * exon_end - strand * exon_start

        if pos > utr3p_len:
            return pos - utr3p_len, offset, True
        if pos <= utr5p_len:
            return pos - utr5p_len - 1, offset, False

        return pos - utr5p_len, offset, False

    def rsid_to_chrom(self, rsid):
        """return the chromosome variants of the rsid
        """
        for var in get_rsid_chromvar(self.conn, rsid):
            ref = var.ref.replace('-', '')
            for alt in var.observed.split('/'):
                alt = alt.replace('-', '')
                yield VCF(self.genome, ChromVariant(var.chrom, var.start + 1, ref, alt)).as_tuple()

    def get_gene_tx(self, gene_id):
        """return all transcripts in the gene
        """
        return get_gene_tx(self.conn, gene_id)

    def get_gene_chrom(self, gene_id):
        """return the chomosome which the gene belongs to
        """
        return get_gene_chrom(self.conn, gene_id)

    def get_refseq_tx(self, refseq):
        """return all transcripts named refseq
        """
        return get_refseq_tx(self.conn, refseq)

    def increase_pos(self, position, length, tx):
        """increase rna position with `length`
        """
        chrom_pos = self.rna_to_chrom_pos(tx, Position(position.pos, position.offset, position.utr3p))
        if not chrom_pos:
            return None
        chrom_pos += tx.strand * length
        pos, offset, utr3p = self.chrom_to_rna_pos(tx, chrom_pos)
        return Position(pos, offset, utr3p)

    def get_intron_pos(self, k, offset, tx):
        """get the rna position of (IVS k)
        """
        if tx.strand == -1:
            exon_starts, exon_ends = rev_p1(tx.exon_ends), rev_p1(tx.exon_starts)
        else:
            exon_starts, exon_ends = tx.exon_starts, tx.exon_ends

        if offset < 0 and k < len(exon_starts):
            chrom_pos = exon_starts[k] + tx.strand
        elif offset > 0 and k - 1 < len(exon_ends):
            chrom_pos = exon_ends[k - 1]
        else:
            return None

        pos, offset, utr3p = self.chrom_to_rna_pos(tx, chrom_pos)
        return Position(pos, offset, utr3p)
