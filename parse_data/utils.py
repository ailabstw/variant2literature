"""utils
"""
import html
import re
import io
import logging
import string
import threading
import functools

import numpy as np
import unidecode

from .greek_alphabet import GREEK_ALPHABETS

GREEK_ALPHABETS_TRANS = str.maketrans({k: v.lower() + ' ' for k, v in GREEK_ALPHABETS.items()})
SEP_PATTERN = re.compile('(?<=[{p}])(?=[^{p}])|(?<=[^{p}])(?=[{p}])'.format(p=string.printable))

logger = logging.getLogger(__name__)


def dump_np(data):
    """dump numpy array
    """
    f = io.BytesIO()
    np.savez(f, data=data)
    f.seek(0)
    return f.read()


def load_np(data):
    """load dumped numpy array
    """
    return np.load(io.BytesIO(data))['data']


def clean_text(x):
    """clean text
    """
    x = html.unescape(x)
    x = SEP_PATTERN.sub(' ', x)
    x = x.replace('\n', ' ')
    x = x.translate(GREEK_ALPHABETS_TRANS)
    x = unidecode.unidecode(x)
    x = ' '.join(x.strip().split())
    x = ''.join(filter(lambda c: c in string.printable, x))
    return x


def overlap_ratio(box1, box2, extend=0):
    """calcuate overlap ratio between two boxes
    """
    h = max(0, min(box1[2], box2[2]) - max(box1[0], box2[0]) + extend)
    w = max(0, min(box1[3], box2[3]) - max(box1[1], box2[1]) + extend)
    area = h * w
    area1 = (box1[2] - box1[0]) * (box1[3] - box1[1])
    area2 = (box2[2] - box2[0]) * (box2[3] - box2[1])
    if area1 <= 0 or area2 <= 0:
        return 0
    return area / min(area1, area2)


class FuncThread(threading.Thread):
    """function thread
    """

    def __init__(self, target, args, kwargs):
        threading.Thread.__init__(self)
        self.target = target
        self.args = args
        self.kwargs = kwargs
        self.result = None

    def run(self):
        self.result = self.target(*self.args, **self.kwargs)


def timeout(time):
    """timeout a function
    """

    def _timeout(fun):
        @functools.wraps(fun)
        def wrapper(*args, **kwargs):
            t = FuncThread(target=fun, args=args, kwargs=kwargs)
            t.start()
            t.join(time)
            if t.isAlive():
                raise TimeoutError
            return t.result
        return wrapper
    return _timeout
