"""util functions
"""
from typing import DefaultDict, List, Dict, NamedTuple, Optional
from collections import defaultdict


RNA_PROTEIN = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "X", "TAG": "X",
    "TGT": "C", "TGC": "C", "TGA": "X", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


THREE_TO_ONE = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    'TER': 'X', 'STOP': 'X', '*': 'X',
}

PROTEIN_RNA: DefaultDict[str, List[str]] = defaultdict(list)
for key, value in RNA_PROTEIN.items():
    PROTEIN_RNA[value].append(key)

ONE_TO_THREE: Dict[str, str] = dict()
for key, value in THREE_TO_ONE.items():
    ONE_TO_THREE[value] = key


def rna_to_protein(seq):
    """convert 3 nucleotide to amino acid
    """
    return RNA_PROTEIN.get(seq)


def protein_to_rna(res):
    """convert amino acid to 3 nucleotide
    """
    return PROTEIN_RNA.get(res, list())


def three_to_one(res):
    """convert 3 code amino acid to 1 code
    """
    return THREE_TO_ONE.get(res, res)


def one_to_three(res):
    """convert 1 code amino acid to 3 code
    """
    return ONE_TO_THREE.get(res, res)


class Position(NamedTuple):
    """a rna position
    """
    pos: int
    offset: int = 0
    utr3p: bool = False

    def __str__(self):
        s = f'{self.pos}'
        if self.offset:
            s += f'{self.offset:+d}'
        if self.utr3p:
            s = '*' + s
        return s


class MTVariant(NamedTuple):
    """a rna variant
    """
    mut_type: str
    start: int
    end: Optional[int] = None
    ref: str = ''
    alt: str = ''
    dup: int = 2


class RNAVariant(NamedTuple):
    """a rna variant
    """
    mut_type: str
    start: Position
    end: Optional[Position] = None
    ref: str = ''
    alt: str = ''
    dup: int = 2


class ProteinVariant(NamedTuple):
    """a protein variant (currently only substitution)
    """
    pos: int
    ref: str
    alt: str


class ChromVariant(NamedTuple):
    """a chromosome variant
    """
    chrom: str
    pos: int
    ref: str
    alt: str


def revcomp(seq):
    """return reversed and complementary sequence
    """
    comp = {
        'A': 'T', 'C': 'G', 'G': 'C',
        'T': 'A', 'N': 'N', '-': '-',
    }
    return ''.join(map(comp.get, reversed(seq)))


def rev_p1(lst):
    """reverse and plus 1
    """
    return list(map(lambda x: x + 1, reversed(lst)))
