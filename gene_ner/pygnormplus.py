"""tmVar in python"""
import json
import re
import logging
from collections import defaultdict

from nltk.stem.snowball import SnowballStemmer
import CRFPP

from . import features
from . import utils
from .gene_normalization import GeneNormalizer
from .gnormplus_pt import GenePrefixTree

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


GENE_PT = GenePrefixTree()


class Extractor():
    """extract variant mentions
    """

    def __init__(self):
        self.tagger = CRFPP.Tagger("-m /app/models/GNR.Model")
        self.normalizer = GeneNormalizer()
        self.stemmer = SnowballStemmer('english')
        self.gene_dict = self.load_gene_symbols()

    def load_gene_symbols(self):
        genes = set()
        with open('/app/models/gene_symbols.txt') as f:
            for line in f:
                genes.add(line.strip('\n'))
        return genes

    def search_gene(self, text):
        text_ = re.sub(r'[^0-9A-Za-z_.\'@+-]', ' ', text)
        for token in text_.split():
            if token in self.gene_dict:
                start = text.find(token)
                end = start + len(token)
                mention = text[start:end]
                gene_id = self.normalizer.normalize_one(mention)
                if gene_id:
                    return (mention, (start, end), gene_id)
        return None

    def filter_gene(self, mention):
        tokens, _ = utils.tokenize(mention.lower())
        if re.match(r'^fig\d+$', ''.join(tokens[:2])):
            return False
        if re.match(r'^chr\d+$', ''.join(tokens[:2])):
            return False
        return bool(GENE_PT.search(tokens))

    def extract(self, text, filename):
        """extract variant mention from text lines
        """
        cnt = 0
        offset, mention_offsets = 0, []
        for line in text.split('\n'):
            for start, end in self.extract_sent(line):
                mention_offsets.append((start + offset, end + offset))
            offset += len(line) + 1
        results = self.postprocess(text, mention_offsets)
        return results

    def _normalize(self, text, offsets, max_text_block=5000):
        ret = []
        j, n = 0, len(text)
        for i in range(0, len(text), max_text_block):
            block_offsets = []
            max_end = i + max_text_block
            while j < len(offsets) and offsets[j][0] < i + max_text_block:
                block_offsets.append((offsets[j][0] - i, offsets[j][1] - i))
                max_end = min(len(text), max(max_end, offsets[j][1]))
                j += 1
            if not block_offsets:
                continue
            ret += self.normalizer.normalize(text[i:max_end], block_offsets)
        return ret

    def postprocess(self, text, offsets):
        # gene_ids = self.normalizer.normalize(text, offsets)
        gene_ids = self._normalize(text, offsets)
        ret = []
        for gene_id, (start, end) in zip(gene_ids, offsets):
            mention = text[start:end]
            if gene_id and self.filter_gene(mention):
                ret.append((mention, (start, end), gene_id))
        return ret

    def extract_sent(self, document):  # pylint: disable=too-many-locals
        """extract variant mention from one sentence
        """
        tokens, offsets = utils.tokenize(document)
        ctd_offsets = GENE_PT.search_tokens(tokens, offsets)

        self.tagger.clear()

        feas_all = []
        for idx, (token, offset) in enumerate(zip(tokens, offsets)):
            f_stem = self.stemmer.stem(token.lower())
            f_ws = features.get_ws_fea(idx, offsets)
            f_num = features.get_count_features(token)
            f_spchar = features.get_regex_fea(token, features.PATTERN_SPECIAL_CHAR)
            f_chem_pre_suf = features.get_regex_fea(token, features.PATTERN_CHEM_PRE_SUF)
            f_mention_type = features.get_mention_type(idx, tokens)

            f_prefix = features.get_prefix(token)
            f_suffix = features.get_suffix(token)
            f_ctd_gene_state = features.get_ctd_gene_state(offset, ctd_offsets)
            f_abb_state = features.get_abb_state(token)

            f_protein_sym = '__nil__'
            _f_protein_sym = features.get_regex_fea(token, features.PATTERN_PROTEIN_SYM)
            if _f_protein_sym != '__nil__':
                f_chem_pre_suf = _f_protein_sym

            feas = ' '.join([token, f_stem, f_ws, f_num, f_spchar, f_chem_pre_suf,
                             f_mention_type, f_protein_sym, f_prefix, f_suffix, f_ctd_gene_state, f_abb_state])
            self.tagger.add(feas)
            feas_all.append(feas)

        self.tagger.parse()

        tags = []
        for i in range(self.tagger.size()):
            tags.append(self.tagger.y2(i))

        mention_offsets = self.process_crf_output(offsets, tags)
        return mention_offsets

    def process_crf_output(self, offsets, tags):
        i, mention_offsets = 0, []
        while i < len(tags):
            if tags[i] == 'Gene_B':
                start = offsets[i][0]
                i += 1
                while i < len(tags) and tags[i] not in {'Gene_B', 'O'}:
                    i += 1
                end = offsets[i - 1][1]
                mention_offsets.append((start, end))
            else:
                i += 1
        return mention_offsets


if __name__ == '__main__':
    e = Extractor()
    text = 'blabla. MYO15A. blabla'
    for start, end in e.extract_sent(text):
        print(text[start:end])
