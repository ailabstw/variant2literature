"""utils
"""
import re
import logging
import threading
import functools
import itertools

logger = logging.getLogger(__name__)


MAX_TOKEN_LEN = 50


def get_offsets(document, tokens):
    idx, offsets = 0, []
    for token in tokens:
        idx = document.find(token, idx)
        offsets.append((idx, idx + len(token)))
        idx += len(token)
    return offsets


def cut_token(token):
    return [token[i:i + MAX_TOKEN_LEN] for i in range(0, len(token), MAX_TOKEN_LEN)]


def tokenize(document):
    doc_copy = document[:]
    doc_copy = re.sub(r'(?<=[A-Za-z])(?=[0-9])', ' ', doc_copy)
    doc_copy = re.sub(r'(?<=[0-9])(?=[A-Za-z])', ' ', doc_copy)
    doc_copy = re.sub(r'(?=[\W\-_])|(?<=[\W\-_])', ' ', doc_copy)

    tokens = doc_copy.split()
    tokens = list(itertools.chain(*map(cut_token, tokens)))
    offsets = get_offsets(document, tokens)
    return tokens, offsets
