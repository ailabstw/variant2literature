"""utils
"""
import re
import itertools
import logging

logger = logging.getLogger(__name__)

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

MAX_TOKEN_LEN = 50


def get_offsets(document, tokens):
    idx, offsets = 0, []
    for token in tokens:
        idx = document.find(token, idx)
        offsets.append((idx, idx + len(token)))
        idx += len(token)
    return offsets


def cut_token(token):
    """split the token if too long
    """
    return [token[i:i + MAX_TOKEN_LEN] for i in range(0, len(token), MAX_TOKEN_LEN)]


def tokenize(document):
    doc_copy = document[:]
    doc_copy = re.sub(r'(?<=[A-Za-z])(?=[0-9])', ' ', doc_copy)
    doc_copy = re.sub(r'(?<=[0-9])(?=[A-Za-z])', ' ', doc_copy)
    doc_copy = re.sub(r'(?<=[A-Z])(?=[a-z])', ' ', doc_copy)
    doc_copy = re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', doc_copy)
    doc_copy = re.sub(r'(?=fs)', ' ', doc_copy)
    doc_copy = re.sub(r'(?=[\W\-_])|(?<=[\W\-_])', ' ', doc_copy)

    tokens = doc_copy.split()
    tokens = list(itertools.chain(*map(cut_token, tokens)))
    offsets = get_offsets(document, tokens)
    return tokens, offsets


def readlines(filename):
    ret = []
    with open(filename) as f:
        for line in f:
            ret.append(line.strip())
    return ret
