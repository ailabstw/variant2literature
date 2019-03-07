"""tmVar in python"""
import json
import re
import logging
from collections import defaultdict

from nltk.stem.snowball import SnowballStemmer
import CRFPP

from . import features
from . import utils

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


PAD_TEXT = 'blablabla, '
MAX_GENE_LEGNTH = 2305000
MAX_RNA_LENGTH = 109224


class Extractor():
    """extract variant mentions
    """

    def __init__(self):
        self.tagger = CRFPP.Tagger("-m /app/models/MentionExtractionUB.Model")
        self.stemmer = SnowballStemmer('english')
        # self.pos_tagger = PerceptronTagger()
        self.regex_dna_mutation_str = utils.readlines('/app/models/tmvar_regexes/DNAMutation.RegEx.txt')
        self.regex_protein_mutation_str = utils.readlines('/app/models/tmvar_regexes/ProteinMutation.RegEx.txt')
        self.regex_snp_mutation_str = utils.readlines('/app/models/tmvar_regexes/SNP.RegEx.txt')

    def extract(self, text, filename):
        """extract variant mention from text lines
        """
        results = []
        offset = 0
        for line in text.split('\n'):
            mentions = self.extract_sent(PAD_TEXT + line + ', ' + PAD_TEXT)
            results += self.postprocess(mentions, filename, offset)
            offset += len(line) + 1
        return results

    def extract_sent(self, document):  # pylint: disable=too-many-locals
        """extract variant mention fomr one sentence
        """
        tokens, offsets = utils.tokenize(document)

        offsets_protein = features.get_hgvs_offsets(document, self.regex_protein_mutation_str)
        offsets_dna = features.get_hgvs_offsets(document, self.regex_dna_mutation_str)
        offsets_snp = features.get_hgvs_offsets(document, self.regex_snp_mutation_str)
        # pos_dict = features.get_pos_tags(tokens, self.pos_tagger)

        self.tagger.clear()

        for idx, (token, offset) in enumerate(zip(tokens, offsets)):

            # f_pos = pos_dict.get(token, '_NULL_')
            f_pos = 'NN'
            f_stem = self.stemmer.stem(token.lower())
            f_num = features.get_count_features(token)

            f_spchar = features.get_regex_fea(token, features.PATTERN_SPECIAL_CHAR)
            f_chrom_key = features.get_regex_fea(token, features.PATTERN_CHROM_KEY)
            f_mut_type = features.get_regex_fea(token.lower(), features.PATTERN_MUTTYPE)
            f_mut_word = features.get_regex_fea(token.lower(), features.PATTERN_MUTWORD)
            f_mut_article = features.get_regex_fea(token.lower(), features.PATTERN_BASE)
            f_type1 = features.get_regex_fea(token.lower(), features.PATTERN_TYPE1)
            f_type2 = features.get_regex_fea(token, features.PATTERN_TYPE2)
            f_dna_sym = features.get_regex_fea(token, features.PATTERN_DNASYM)
            f_rs_code = features.get_regex_fea(token, features.PATTERN_RSCODE)

            prev_token = tokens[idx - 1] if idx > 0 else ''
            prev_char = document[offset[0] - 1] if offset[0] > 0 else ''
            f_protein_sym = features.get_protein_symbol(token, prev_token, prev_char)

            f_pattern = features.get_pattern(token)
            f_prefix = features.get_prefix(token)
            f_suffix = features.get_suffix(token)
            f_hgvs = features.get_hgvs_feature(offset, offsets_protein,
                                               offsets_dna, offsets_snp)

            feas = ' '.join([token, f_stem, f_pos, f_num, f_spchar, f_chrom_key, f_mut_type,
                             f_mut_word, f_mut_article, f_type1, f_type2, f_dna_sym, f_protein_sym,
                             f_rs_code, f_pattern, f_prefix, f_suffix, f_hgvs])
            self.tagger.add(feas)

        self.tagger.parse()

        tags = []
        for i in range(self.tagger.size()):
            tags.append(self.tagger.y2(i))

        if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
            for tag, token in zip(tags, tokens):
                print(token, tag)

        mentions = self.get_mention_from_preds(document, tokens, offsets, tags)
        return mentions

    def mention_end_cond(self, mention, tag, prev_token):
        """conditions to determine an end of a variant mention
        """
        first_tag = mention[0][1]

        if first_tag == 'R' and tag == 'R' and not prev_token.isdigit():
            return False
        if first_tag != 'R' and tag not in ('O', 'A', 'R'):
            return False
        return True

    def merge_tokens(self, mention):
        """merge tokens with the same labels

        ex: A(W) rg(W) 123(P) I(M) le(M) --> Arg(W) 123(P) Ile(M)
        """
        merged_tokens, merged_tags = [], []
        idx = 0
        while idx < len(mention):
            _, cur_tag = mention[idx]
            cur_token = ''
            while idx < len(mention):
                token, tag = mention[idx]
                if tag != cur_tag:
                    break
                cur_token += token
                idx += 1
            merged_tokens.append(cur_token)
            merged_tags.append(cur_tag)
        return list(zip(merged_tokens, merged_tags))

    def get_mention_from_preds(self, document, tokens, offsets, tags):
        """get variant mentions from CRF preditions
        """
        idx, ret = 0, []
        while idx < len(tokens):
            if tags[idx] != 'O':
                start = offsets[idx][0]
                mention = [(tokens[idx], tags[idx])]
                idx += 1
                while idx < len(tokens) and not self.mention_end_cond(mention, tags[idx], tokens[idx - 1]):
                    mention.append((tokens[idx], tags[idx]))
                    idx += 1
                end = offsets[idx - 1][1]
                ret.append((document[start:end], self.merge_tokens(mention), (start, end)))
            else:
                idx += 1
        return ret

    def get_element(self, tokens, k, default=''):
        """get k-th elements
        """
        if len(tokens) <= k:
            return default
        return tokens[k]

    def clean_residue(self, residues):
        """clean wild type and mutant string
        """
        ret = []
        for residue in residues:
            residue = re.sub(r'[^A-Za-z0-9*]', '', residue)
            if residue:
                ret.append(residue.upper())
        return ret

    def clean_position(self, positions):
        """clean position string
        """
        ret = []
        for pos in positions:
            for part in pos.split('_'):
                part = re.sub(r'[()\[\]]', '', part)
                if part and re.match('[A-Za-z0-9*+-]', part):
                    ret.append(part)
        return ret

    def clean_seq_type(self, seq_types):
        """clean seq type string
        """
        ret = []
        for seq_type in seq_types:
            seq_type = seq_type.lower()
            if seq_type in ['c', 'n', 'p', 'g', 'chr']:
                ret.append(seq_type)
        return ret

    def check_mutation_type(self, tag_dict):  # pylint: disable=too-many-statements
        """check mutation type
        """
        tag_dict['T'] = list(map(str.upper, tag_dict['T']))
        tag_dict['A'] = self.clean_seq_type(tag_dict['A'])
        tag_dict['P'] = self.clean_position(tag_dict['P'])
        tag_dict['W'] = self.clean_residue(tag_dict['W'])
        tag_dict['M'] = self.clean_residue(tag_dict['M'])

        has_del = any('DEL' in t for t in tag_dict['T'])
        has_ins = any('INS' in t for t in tag_dict['T'])
        has_dup = any('DUP' in t for t in tag_dict['T'])

        var_dict = {}

        if 'R' in tag_dict:
            var_dict['mut_type'] = 'RSID'
            var_dict['rsid'] = self.get_element(tag_dict['R'], 0)
            return var_dict

        var_dict['seq_types'] = tag_dict['A'][:1]

        if 'F' in tag_dict:
            var_dict['mut_type'] = 'FS'
            var_dict['wild_type'] = self.get_element(tag_dict['W'], 0)
            var_dict['start'] = self.get_element(tag_dict['P'], 0)
            var_dict['end'] = var_dict['start']
            var_dict['mutant'] = self.get_element(tag_dict['M'], 0)
            var_dict['frameshift'] = self.get_element(tag_dict['S'], 0)

        elif has_del and has_ins:
            var_dict['mut_type'] = 'INDEL'
            var_dict['start'] = self.get_element(tag_dict['P'], 0)
            var_dict['end'] = self.get_element(tag_dict['P'], 1, default=var_dict['start'])
            var_dict['mutant'] = self.get_element(tag_dict['M'], 0)

        elif has_del:
            var_dict['mut_type'] = 'DEL'
            var_dict['start'] = self.get_element(tag_dict['P'], 0)
            var_dict['end'] = self.get_element(tag_dict['P'], 1, default=var_dict['start'])
            var_dict['wild_type'] = self.get_element(tag_dict['M'], 0)

        elif has_ins:
            var_dict['mut_type'] = 'INS'
            var_dict['start'] = self.get_element(tag_dict['P'], 0)
            var_dict['end'] = self.get_element(tag_dict['P'], 1, default=var_dict['start'])
            var_dict['mutant'] = self.get_element(tag_dict['M'], 0)

        elif has_dup:
            var_dict['mut_type'] = 'DUP'
            var_dict['start'] = self.get_element(tag_dict['P'], 0)
            var_dict['end'] = self.get_element(tag_dict['P'], 1, default=var_dict['start'])
            var_dict['dup'] = self.get_element(tag_dict['D'], 0)
            var_dict['wild_type'] = self.get_element(tag_dict['M'], 0)

        else:
            var_dict['mut_type'] = 'SUB'
            var_dict['wild_type'] = self.get_element(tag_dict['W'], 0)
            var_dict['start'] = self.get_element(tag_dict['P'], 0)
            var_dict['end'] = var_dict['start']
            var_dict['mutant'] = self.get_element(tag_dict['M'], 0)

        if not var_dict['seq_types']:
            self.get_possible_seq_type(var_dict, tag_dict)

        if 'mutant' in var_dict:
            var_dict['mutant'] = utils.ALL_TO_ONE.get(var_dict['mutant'], var_dict['mutant'])
        if 'wild_type' in var_dict:
            var_dict['wild_type'] = utils.ALL_TO_ONE.get(var_dict['wild_type'], var_dict['wild_type'])

        return var_dict

    def get_possible_seq_type(self, var_dict, tag_dict):
        """return possible seq types
        """
        seq_types = {'p', 'c', 'n', 'g', 'chr'}

        if re.search('[^ACGT0-9]', var_dict.get('wild_type', '')):
            seq_types = {'p'}

        if re.search('[^ACGT0-9]', var_dict.get('mutant', '')):
            seq_types = {'p'}

        if var_dict['mut_type'] == 'FS':
            seq_types = {'p'}

        for pos in [var_dict['start'], var_dict['end']]:
            seq_types &= self.check_pos_length(pos)

        if 'p' in seq_types and tag_dict['mention_type'] == 'protein_change':
            seq_types = {'p'}

        var_dict['seq_types'] = []
        for sts in [['c', 'n'], ['c'], ['p'], ['n'], ['g'], ['chr']]:
            if all(st in seq_types for st in sts):
                var_dict['seq_types'] = sts
                break

    def check_pos_length(self, pos):
        """check if the position is too large
        """
        if pos.upper().startswith('IVS') or pos.upper().startswith('INTRON'):
            return {'c', 'n'}

        mch = re.match(r'[*+-]{0,1}[0-9]+', pos)
        if not mch:
            return set()

        main_pos = mch.group(0)
        if main_pos[0] == '*':
            return {'c', 'n'}
        main_pos = int(mch.group(0))

        if main_pos > MAX_GENE_LEGNTH:
            return {'chr'}
        if main_pos > MAX_RNA_LENGTH:
            return {'g', 'chr'}
        if main_pos > MAX_RNA_LENGTH // 3:
            return {'c', 'n', 'g', 'chr'}

        return {'p', 'c', 'n', 'g', 'chr'}

    def check_variant(self, tag_dict, text, filename):
        """check whether the mention a valid variant
        """
        if tag_dict['mut_type'] == 'SUB':
            if not tag_dict['wild_type'] or not tag_dict['mutant']:
                logger.debug('sub no W / M: (id:%s) %s', filename, text)
                return False

        if tag_dict['mut_type'] != 'RSID':
            if 'seq_types' not in tag_dict or not tag_dict['seq_types']:
                logger.debug('no sequence types: (id:%s) %s', filename, text)
                return False

            if 'start' not in tag_dict or not tag_dict['start']:
                logger.debug('no start position: (id:%s) %s', filename, text)
                return False

        if tag_dict['mut_type'] in ['INS', 'INDEL']:
            if re.search('[^A-Z]', tag_dict['mutant']):
                logger.debug('wrong mutant for insertion: (%s) %s', filename, text)
                return False

        return True

    def get_mention_type(self, mentions) -> str:
        """return mention type

        Returns:
            `protein_change` (ex: Arg123Gln) or `others`
        """
        tag_seq = ''.join(tag for _, tag in mentions)
        if 'WPM' in tag_seq:
            return 'protein_change'
        return 'others'

    def postprocess(self, mentions, filename, offset):
        """post processing
        """
        ret = []
        for text, mention, (start, end) in mentions:
            tag_dict = defaultdict(list)
            for token, tag in mention:
                tag_dict[tag].append(token)
            tag_dict['mention_type'] = self.get_mention_type(mention)
            tag_dict = self.check_mutation_type(tag_dict)
            if self.check_variant(tag_dict, text, filename):
                start += offset - len(PAD_TEXT)
                end += offset - len(PAD_TEXT)
                ret.append((text, (start, end), json.dumps(tag_dict)))
        return ret
