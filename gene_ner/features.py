"""CRF features"""
from typing import List, Dict, Tuple, Pattern
import re
import bisect


PATTERN_SPECIAL_CHAR = [
    (re.compile(r'^[;:,.>+_-]$'), '-SpecificC1-'),
    (re.compile(r'^[()]$'), '-SpecificC2-'),
    (re.compile(r'^[{}]$'), '-SpecificC3-'),
    (re.compile(r'^[[]]$'), '-SpecificC4-'),
    (re.compile(r'^[\\/]$'), '-SpecificC5-'),
]


PATTERN_CHEM_PRE_SUF = [
    (re.compile(r'^.*(yl|ylidyne|oyl|sulfonyl)$'), '-CHEMinlineSuffix-'),
    (re.compile(r'^(meth|eth|prop|tetracos).*$'), '-CHEMalkaneStem-'),
    (re.compile(r'^(di|tri|tetra).*$'), '-CHEMsimpleMultiplier-'),
    (re.compile(r'^(benzen|pyridin|toluen).*$'), '-CHEMtrivialRing-'),
    (re.compile(r'^.*(one|ol|carboxylic|amide|ate|acid|ium|ylium|ide|uide|iran|olan|inan|pyrid|acrid|'
                 'amid|keten|formazan|fydrazin)(s|)$'),
     '-CHEMsuffix-'),
]

PATTERN_PROTEIN_SYM = [
    (re.compile(r'^.*(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|glycine|aspartate|'
                 'methionine|threonine|histidine|aspartic|asparticacid|arginine|asparagine|tryptophan|'
                 'proline|phenylalanine|cysteine|serine|glutamate|tyrosine|stop|frameshift).*$'),
     '-ProteinSymFull-'),
    (re.compile(r'^(cys|ile|ser|gln|met|asn|pro|lys|asp|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr|fs|fsx)$'), '-ProteinSymTri-'),
    (re.compile(r'^[CISQMNPKDTFAGHLRWVEYX]$'), '-ProteinSymChar-'),
]


def get_ctd_gene_state(offset, ctd_offsets):
    """check the offset of some token in any HGVS mention
    """
    idx = bisect.bisect_right(ctd_offsets, (offset[0], float('inf')))
    if idx > 0 and offset[0] < ctd_offsets[idx - 1][1]:
        if offset[0] == ctd_offsets[idx - 1][0]:
            return 'CTDGene_B'
        elif offset[1] == ctd_offsets[idx - 1][1]:
            return 'CTDGene_E'
        return 'CTDGene_I'
    return '__nil__'


def get_abb_state(token):
    return '__nil__'

def get_mention_type(idx, tokens):
    token = tokens[idx]
    if re.match(r'^(ytochrome|cytochrome)$', token):
        return '-Type_cytochrome-'

    if re.match(r'^.*target$', token):
        return '-Type_target-'

    if re.match(r'^.*(irradiation|hybrid|fusion|experiment|gst|est|gap|antigen)$', token):
        return '-Type_ExperimentNoun-'

    if re.match(r'^.*(irradiation|hybrid|fusion|experiment|gst|est|gap|antigen)$', token):
        return '-Type_ExperimentNoun-'

    if re.match(r'^.*(irradiation|hybrid|fusion|experiment|gst|est|gap|antigen)$', token):
        return '-Type_ExperimentNoun-'

    if re.match(r'^.*(disease|disorder|dystrophy|deficiency|syndrome|dysgenesis|cancer|injury|neoplasm|diabetes|diabete)$', token):
        return '-Type_Disease-'

    if re.match(r'^.*(motif|domain|omain|binding|site|region|sequence|frameshift|finger|box).*$', token):
        return '-Type_DomainMotif-'

    if (token == '-' and idx < len(tokens) - 1 and
            re.match(r'^.*(motif|domain|omain|binding|site|region|sequence|frameshift|finger|box).*$', tokens[idx + 1])):
        return '-Type_DomainMotif-'

    if re.match(r'^[rmc]$', token) and idx < len(tokens) - 1 and tokens[idx + 1] in {'DNA', 'RNA'}:
        return '-Type_DomainMotif-'

    if re.match(r'^.*(famil|complex|cluster|proteins|genes|factors|transporter|proteinase|membrane|ligand|enzyme|channels|tors|ase|ases)$', token):
        return '-Type_Family-'

    if re.match(r'^marker', token.lower()):
        return '-Type_Marker-'

    if (token == '.*cell.*' or (idx < len(tokens) - 1 and tokens[idx + 1] == 'cell' and
                                           re.match(r'^(T|B|monocytic|cancer|tumor|myeloma|epithelial|crypt)$', token))):
        return '-Type_Cell-'

    if token == '.*chromosome.*':
        return '-Type_Chromosome-'

    if (re.match(r'^[pq]$', token) and ((idx < len(tokens) - 1 and re.match('^[0-9]+$', tokens[idx + 1])) or
                                        (idx > 0 and re.match('^[0-9]+$', tokens[idx - 1])))):
        return '-Type_ChromosomeStrain-'

    if re.match(r'^.*(related|regulated|associated|correlated|reactive).*$', token):
        return '-Type_relation-'

    if re.match(r'^.*(polymorphism|mutation|deletion|insertion|duplication|genotype|genotypes).*$', token.lower()):
        return '-Type_VariationTerms-'

    if re.match(r'.*(oxidase|transferase|transferases|kinase|kinese|subunit|unit|receptor|adrenoceptor|transporter|'
                 'regulator|transcription|antigen|protein|gene|factor|member|molecule|channel|deaminase|spectrin).*', token):
        return '-Type_suffix-'

    if (re.match(r'^[-_(]$', token) and idx < len(tokens) - 1 and
            re.match(r'^.*(alpha|beta|gamma|delta|theta|kappa|zeta|sigma|omega|i|ii|iii|iv|v|vi|[abcdefgyr])$', tokens[idx + 1].lower())):
        return '-Type_strain-'

    if re.match(r'^(alpha|beta|gamma|delta|theta|kappa|zeta|sigma|omega|i|ii|iii|iv|v|vi|[abcdefgyr])$', token):
        return '-Type_strain-'

    return '__nil__'


def get_ws_fea(idx, offsets):
    if idx == 0:
        f_wsb = 'WSB:1st'
    elif offsets[idx][0] == offsets[idx - 1][1]:
        f_wsb = 'WSB:NoGap'
    else:
        f_wsb = 'WSB:Gap'

    if idx == len(offsets) - 1:
        f_wsf = 'WSF:last'
    elif offsets[idx][1] == offsets[idx + 1][0]:
        f_wsf = 'WSF:NoGap'
    else:
        f_wsf = 'WSF:Gap'
    return ' '.join([f_wsb, f_wsf])


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
