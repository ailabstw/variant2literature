"""CRF features"""
from typing import List, Dict, Tuple, Pattern
import re
import bisect

from nltk.tag.api import TaggerI


PATTERN_SPECIAL_CHAR = [
    (re.compile(r'^[;:,.->+_]$'), '-SpecificC1-'),
    (re.compile(r'^[()]$'), '-SpecificC2-'),
    (re.compile(r'^[{}]$'), '-SpecificC3-'),
    (re.compile(r'^[[]]$'), '-SpecificC4-'),
    (re.compile(r'^[\\/]$'), '-SpecificC5-'),
]

PATTERN_CHROM_KEY = [
    (re.compile(r'^(q|p|qter|pter|XY|t)$'), '-ChroKey-'),
]

PATTERN_MUTTYPE = [
    (re.compile(r'^(del|ins|dup|tri|qua|con|delins|indel)$'), '-MutatType-'),
    (re.compile(r'(fs|fsX|fsx)'), '-FrameShiftType-')
]

PATTERN_MUTWORD = [
    (re.compile(r'^(deletion|delta|elta|insertion|repeat|inversion|deletions|insertions|'
                'repeats|inversions)$'),
     '-MutatWord-'),
]

PATTERN_BASE = [
    (re.compile(r'^(single|a|one|two|three|four|five|'
                'six|seven|eight|nine|ten|[0-9]+)$'), '-Base-'),
    (re.compile(r'^(kb|mb)$'), '-Byte-'),
    (re.compile(r'(base|bases|pair|amino|acid|acids|codon|'
                'postion|postions|bp|nucleotide|nucleotides)'), '-bp-'),
]

PATTERN_TYPE1 = [
    (re.compile(r'^[cgrm]$'), '-Type1-'),
    (re.compile(r'^(ivs|ex|orf)$'), '-Type1_2-'),
]

PATTERN_TYPE2 = [
    (re.compile(r'^p$'), '-Type2-'),
]

PATTERN_DNASYM = [(re.compile(r'^[ATCGUatcgu]$'), '-DNASym-')]

PATTERN_RSCODE = [(re.compile(r'^(rs|RS|Rs)$'), '-RScode-')]

PATTERN_PROTEIN_FULL = (
    re.compile(r'(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|'
               'glycine|aspartate|methionine|threonine|histidine|aspartic|'
               'asparticacid|arginine|asparagine|tryptophan|proline|phenylalanine|'
               'cysteine|serine|glutamate|tyrosine|stop|frameshift)'), '-ProteinSymFull-',
)
PATTERN_PROTEIN_TRI = (
    re.compile(r'^(cys|ile|ser|gln|met|asn|pro|lys|asp'
               '|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr|fs|fsx)$'), '-ProteinSymTri-',
)
PATTERN_PROTEIN_TRI_SUB = (
    re.compile(r'^(ys|le|er|ln|et|sn|ro|ys|sp|hr|he|la|ly|is|eu|rg|rp|al|lu|yr)$'),
    '-ProteinSymTriSub-',
)
PATTERN_PROTEIN_CHAR = (re.compile(r'^[CISQMNPKDTFAGHLRWVEYX]$'), '-ProteinSymChar-')


def get_hgvs_offsets(document: str, patterns: List[Pattern]) -> List[Tuple[int, int]]:
    """find offsets in the text that match HGVS patterns
    """
    offsets = []
    for pattern in patterns:
        for m in re.finditer(pattern, document):
            offsets.append((m.start(2), m.end(2)))
    offsets.sort()
    return offsets


def get_hgvs(offset: Tuple[int, int],
             offsets: List[Tuple[int, int]],
             name: str) -> str:
    """check the offset of some token in any HGVS mention
    """
    idx = bisect.bisect_right(offsets, (offset[0], float('inf')))
    if idx > 0 and offset[0] < offsets[idx - 1][1]:
        return name
    return 'O'


def get_count_features(token: str) -> str:
    """count digits and letters
    """
    num_d = sum(c.isdigit() for c in token)
    f_num_d = 'N:{}'.format(num_d if num_d < 4 else '4+')

    num_uc = sum(c.isupper() for c in token)
    f_num_uc = 'U:{}'.format(num_uc if num_uc < 4 else '4+')

    num_lc = sum(c.islower() for c in token)
    f_num_lc = 'L:{}'.format(num_lc if num_lc < 4 else '4+')

    num_c = len(token)
    f_num_c = 'A:{}'.format(num_c if num_c < 4 else '4+')
    return ' '.join([f_num_d, f_num_uc, f_num_lc, f_num_c])


def get_hgvs_feature(offset: Tuple[int, int],
                     offsets_protein: List[Tuple[int, int]],
                     offsets_dna: List[Tuple[int, int]],
                     offsets_snp: List[Tuple[int, int]]) -> str:
    """whether the token matches any HGVS pattern
    """
    ret = get_hgvs(offset, offsets_protein, 'ProteinMutation')
    ret = get_hgvs(offset, offsets_dna, 'DNAMutation') if ret == 'O' else ret
    ret = get_hgvs(offset, offsets_snp, 'SNP') if ret == 'O' else ret
    return ret


def get_regex_fea(token: str, patterns: List[Tuple[Pattern, str]]) -> str:
    """mutation regex feature

    Args:
        token: a token from the tokenized text
        patterns: patterns to match
    Returns:
        the pattern name that the token matches
    """
    for pattern, name in patterns:
        if pattern.search(token):
            return name
    return '__nil__'


def get_protein_symbol(token: str, prev_token: str, prev_char: str) -> str:
    """protein symbol feature

    Args:
        token: a token from the tokenized text
        prev_token: the previous token before `token`
        prev_char: the previous character before `token`
    Returns:
        whether the token is a protein symbol
    """
    if PATTERN_PROTEIN_FULL[0].search(token):
        return PATTERN_PROTEIN_FULL[1]

    elif PATTERN_PROTEIN_TRI[0].search(token):
        return PATTERN_PROTEIN_TRI[1]

    elif (PATTERN_PROTEIN_TRI_SUB[0].search(token) and
          PATTERN_PROTEIN_CHAR[0].search(prev_token) and
          prev_char != ' '):
        return PATTERN_PROTEIN_TRI_SUB[1]

    elif PATTERN_PROTEIN_CHAR[0].search(token):
        return PATTERN_PROTEIN_CHAR[1]

    return '__nil__'


def get_pos_tags(tokens: List[str], tagger: TaggerI) -> Dict[str, str]:
    """POS tags features

    Args:
        tokens: the tokens from the tokenized text
    Returns:
        space concatenated pattern features
    """
    tags = tagger.tag(tokens)

    pos_dict = {}
    for token, tag in tags:
        if re.search(r'[^0-9A-Za-z]', token):
            pos_dict[token] = token
        else:
            pos_dict[token] = tag
    return pos_dict


def get_pattern(token: str) -> str:
    """digits, letters pattern features

    Args:
        token: a token from the tokenized text
    Returns:
        space concatenated pattern features
    """
    f_p1 = re.sub(r'[A-Z]', 'A', token)
    f_p1 = re.sub(r'[a-z]', 'a', f_p1)
    f_p1 = re.sub(r'[0-9]', '0', f_p1)
    f_p1 = 'P1:' + f_p1

    f_p2 = re.sub(r'[A-Za-z]', 'a', token)
    f_p2 = re.sub(r'[0-9]', '0', f_p2)
    f_p2 = 'P2:' + f_p2

    f_p3 = re.sub(r'[A-Z]+', 'A', token)
    f_p3 = re.sub(r'[a-z]+', 'a', f_p3)
    f_p3 = re.sub(r'[0-9]+', '0', f_p3)
    f_p3 = 'P3:' + f_p3

    f_p4 = re.sub(r'[A-Za-z]+', 'a', token)
    f_p4 = re.sub(r'[0-9]+', '0', f_p4)
    f_p4 = 'P4:' + f_p4
    return ' '.join([f_p1, f_p2, f_p3, f_p4])


def get_prefix(token: str) -> str:
    """prefix features

    Args:
        token: a token from the tokenized text
    Returns:
        space concatenated prefix features
    """
    return ' '.join(token[:i] if len(token) >= i else '__nil__' for i in range(1, 6))


def get_suffix(token: str) -> str:
    """suffix features

    Args:
        token: a token from the tokenized text
    Returns:
        space concatenated suffix features
    """
    return ' '.join(token[-i:] if len(token) >= i else '__nil__' for i in range(1, 6))
