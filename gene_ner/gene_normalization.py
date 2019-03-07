import re
import csv
import time

# symbol_dict: symbol -> {name: gene_id}，gene_id為string，symbol都是一對一
# synonym_dict: synonyms -> {name: [gene_ids]}
# multiple_name_dict: full name + other names -> {name: [gene_ids]}
# single, multiple的name只留字母以及數字，沒有' '
# gene_ids為list，因為一個名字可能對應多個id
# gene_id_set_dict: all names -> {gene id: set(name_split)}
# name都經過normalize處理，只留下字母、數字，且字母數字之間加入' '
# gene_id_multiple_dict: {gene id:[names]}，有' '

def ref_dict(file):
    symbol_dict = {}
    synonym_dict = {}
    multiple_name_dict = {}
    gene_id_set_dict = {}
    gene_id_multiple_dict = {}

    with open(file, 'r') as f:
        rows = csv.DictReader(f)

        for row in rows:
            gene_id = row['gene_id']
            single_name = []
            multiple_name = []
            gene_id_multiple = []
            gene_id_set = set()

            # create single_name_dict, symbol and synonyms
            symbol_nor = name_normalize(row['symbol'])
            symbol_ns = re.sub('\s','',symbol_nor)
            symbol_split = symbol_nor.split(' ')
            symbol_dict[symbol_ns] = gene_id
            gene_id_set = gene_id_set | set(symbol_split)

            synonyms = row['synonyms']
            if synonyms != '-':
                synonyms_list = synonyms.split('|')
                for i in synonyms_list:
                    i_nor = name_normalize(i)
                    i_ns = re.sub('\s','',i_nor)
                    i_split = i_nor.split(' ')
                    single_name.append(i_ns)
                    gene_id_set = gene_id_set | set(i_split)
            synonym_dict = create_dict(synonym_dict, single_name, gene_id)

            # create multiple_name_dict, full name and other names
            full_name = row['full_name']
            if full_name != '-':
                full_name_nor = name_normalize(full_name)
                full_name_ns = re.sub('\s','',full_name_nor)
                full_name_split = full_name_nor.split(' ')
                multiple_name.append(full_name_ns)
                gene_id_set = gene_id_set | set(full_name_split)
                gene_id_multiple.append(full_name_split)

            other_names = row['other_names']
            if other_names != '-':
                other_names_list = other_names.split('|')
                for name in other_names_list:
                    name_nor = name_normalize(name)
                    name_ns = re.sub('\s','',name_nor)
                    name_split = name_nor.split(' ')
                    multiple_name.append(name_ns)
                    gene_id_set = gene_id_set | set(name_split)
                    gene_id_multiple.append(name_split)

            multiple_name_dict = create_dict(multiple_name_dict, multiple_name, gene_id)

            # create gene_id_set_dict, full name and other names
            gene_id_set_dict[gene_id] = gene_id_set
            gene_id_multiple_dict[gene_id] = gene_id_multiple

    return symbol_dict, synonym_dict, multiple_name_dict, gene_id_set_dict, gene_id_multiple_dict


# 把name對應gene_id加到dict中
# d: dict, li: list, gene_id: str
def create_dict(d, li, gene_id):
    if li == []:
        return d
    for i in li:
        if i != '':
            if d.get(i):
                if gene_id not in d[i]:
                    d[i].append(gene_id)
            else:
                d[i] = [gene_id]
    return d


# 把所有非數字英文字母的字元都換成' '，並把數字與字母間都加入' '
# input and output are strings
def name_normalize(name):
    name = re.sub('(?<=\D)(?=\d)|(?<=\d)(?=\D)', ' ', name)
    new_name = re.sub('[^0-9a-z]+',' ',name.lower())
    if new_name[-1] == ' ':
        new_name = new_name[:-1]
    if new_name[-4:] == 'gene':
        new_name = new_name[:-4]
    return new_name


# 找test是否有在ref_dict中
# 比較test字首有無h
# 有找到就回傳gene_id，無就回傳[]
def simple_find(test_name, ref_dict):
    test_woh = remove_h(test_name) #有些字首的h代表human
    current = ref_dict.get(test_name)
    if current:
        return current

    else:
        current = ref_dict.get(test_woh)
        if current:
            return current
        else:
            return []


# 移除字首的h，有些可能代表human
def remove_h(name):
    if name != '':
        if name[0] == 'h':
            return name[1:]
    return name


# 找test是否有在ref_dict
# test: string
# ref_dict: dict of {ref:gene}
# ref_split_dict: dict of {ref: split ref}
# key is the key_list of ref_dict
def wordbag_find(test_name, gene_id_set_dict, gene_id_key):
    test_set = set(test_name.split(' ')) - {'', 'of', 'the'}
    test_set_reduced = set(re.sub('[0-9\s]+',' ', test_name).split(' ')) - set('abcdefghijklmnopqrstuvwxyz') - {'', 'of', 'an', 'the', 'and', 'gene', 'pseudogene', 'protein', 'receptor'}

    max_count1 = 0
    possible_id1 = []
    max_count2 = 0
    possible_id2 = []

    for gene_id in gene_id_key:
        count1 = len(test_set & gene_id_set_dict[gene_id])
        if count1 > max_count1:
            max_count1 = count1
            possible_id1 = [gene_id]
        elif (count1 == max_count1) & (count1 > 0):
            possible_id1.append(gene_id)

        count2 = len(test_set_reduced & gene_id_set_dict[gene_id])
        if count2 > max_count2:
            max_count2 = count2
            possible_id2 = [gene_id]
        elif (count2 == max_count2) & (count2 > 0):
            possible_id2.append(gene_id)

    if 0 < len(possible_id1) < 10:
        return possible_id1
    else:
        return possible_id2


def multiple_name_find(test_name, candidate, gene_id_multiple_dict):
    test_set = set(test_name.split(' ')) - {''}
    max_count = 0
    possible_id = []

    for gene_id in candidate:
        for ref_split in gene_id_multiple_dict[gene_id]:
            ref_set = set(ref_split)
            count = len(test_set & ref_set)
            if count > max_count:
                max_count = count
                possible_id = [gene_id]
            elif (count == max_count) & (count > 0):
                possible_id.append(gene_id)

    if possible_id == []:
        return candidate
    else:
        return possible_id


# 找最相近的gene name
# 這裡的test_list是沒有title & abstract
def find_near(test, test_list):
    start = int(test[0])
    end = int(test[1])
    before = 0
    after = 10000
    before_near = ''
    after_near = ''
    for ref in test_list:
        ref_start = int(ref[0])
        ref_end = int(ref[1])
        if ref_end < start:
            if ref_end > before:
                before = ref_end
                before_near = [ref_start, ref_end]
        if ref_start > end:
            if ref_start < after:
                after = ref_start
                after_near = [ref_start, ref_end]
    if ((start-before) < (after-end)) & ((start-before) < 5):
        return before_near
    elif ((start-before) >= (after-end)) & ((after-end) < 5):
        return after_near
    else:
        return ''



class GeneNormalizer:
    def __init__(self, ref_file='/app/models/human_gene_data.csv'):
        self.dicts = ref_dict(ref_file)

    def normalize_one(self, text):
        gene_id = self.normalize(text, [(0, len(text))])[0]
        if gene_id:
            return int(gene_id)
        return None

    def normalize(self, text, test_list):
        results = self.answer(text, '', test_list)
        tmp = []
        for ret in results:
            if not ret:
                tmp.append(None)
            elif isinstance(ret, list):
                tmp.append(int(ret[0]))
            else:
                tmp.append(int(ret))
        return tmp


    def answer(self, title, abstract, test_list):
        # reference
        symbol_dict, synonym_dict, multiple_name_dict, gene_id_set_dict, gene_id_multiple_dict = self.dicts
        gene_id_key = list(gene_id_set_dict.keys())

        article = (title + '\n' + abstract).lower()
        current_dict = {}
        mulid_test_list = []
        noans_test_list = []
        final_data = []

        i = -1
        for test in test_list:
            i += 1

            test_start = test[0]
            test_end = test[1]
            test_name = article[test_start: test_end]
            test_nor = name_normalize(test_name)
            test_ns = re.sub('\s','',test_nor)

            # step 1: 從本篇paper中已經找好的test name著手
            possible_id = current_dict.get(test_name) # 若有，possible_id為string
            if possible_id:
                final_data.append(possible_id)
                continue

            # step 2: 從symbol_dict找完全對應者
            possible_id = simple_find(test_ns, symbol_dict) # return a string or []
            if possible_id != []:
                current_dict[test_name] = possible_id
                final_data.append(possible_id)
                continue

            # step 3: 從synonym_dict找完全對應者
            possible_id = simple_find(test_ns, synonym_dict) # return a list
            length = len(possible_id)
            if length == 1:
                current_dict[test_name] = possible_id[0]
                final_data.append(possible_id[0])
                continue
            if length > 1:
                final_data.append('')
                mulid_test_list.append([test_start, test_end, test_name, possible_id, i]) # 存進mulid_test_list
                continue

            # step 4: 從multiple_name_dict找完全對應者
            possible_id = simple_find(test_ns, multiple_name_dict) # return a list
            length = len(possible_id)
            if length == 1:
                current_dict[test_name] = possible_id[0]
                final_data.append(possible_id[0])
                continue
            if length > 1:
                mulid_test_list.append([test_start, test_end, test_name, possible_id, i]) # 存進mulid_test_list
                final_data.append('')
                continue

            # step 5: 用word bag找最相關
            possible_id = wordbag_find(test_nor, gene_id_set_dict, gene_id_key)
            length = len(possible_id)
            if length == 1:
                current_dict[test_name] = possible_id[0]
                final_data.append(possible_id[0])
                continue
            if length > 1:
                # step 6: 每個name去比對，找出最相關的name（比word bag範圍還小）
                possible_id2 = multiple_name_find(test_nor, possible_id, gene_id_multiple_dict)
                length = len(possible_id2)
                if length == 1:
                    current_dict[test_name] = possible_id2[0]
                    final_data.append(possible_id2[0])
                    continue
                if length > 1:
                    mulid_test_list.append([test_start, test_end, test_name, possible_id, i]) # 存進mulid_test_list
                    final_data.append('')
                    continue

            # 剩餘的先視為no ans，最後再處理
            noans_test_list.append([test_start, test_end, test_name, '', i])
            final_data.append('')


        # 處理multiple
        current_gene_id = list(current_dict.values())
        mulid_test_list2 = []

        for test in mulid_test_list:
            test_name = test[2]
            candidate = test[3]
            position = test[4]
            possible_id1 = []
            possible_id2 = []

            # step 7: 比對現有的gene_id
            for gene_id in candidate:
                if gene_id in current_gene_id:
                    possible_id1.append(gene_id)

            length = len(possible_id1)
            if length == 1:
                current_dict[test_name] = possible_id1[0]
                final_data[position] = possible_id1[0]
                continue
            if length > 1:
                mulid_test_list2.append([test[0], test[1], test_name, possible_id1, position])
                continue

            # step 8: 比對文章和gene_id_set_dict的相關性
            max_similar = 0
            article_set = set(name_normalize(article).split(' ')) - {'','the','of'}
            for gene_id in candidate:
                gene_id_set = gene_id_set_dict[gene_id]
                similar = len(article_set & gene_id_set)
                if similar > max_similar:
                    max_similar = similar
                    possible_id2 = [gene_id]
                elif (similar == max_similar) & (similar != 0):
                    possible_id2.append(gene_id)

            length = len(possible_id2)
            if length == 1:
                current_dict[test_name] = possible_id2[0]
                final_data[position] = possible_id2[0]
                continue
            if length > 1:
                mulid_test_list2.append([test[0], test[1], test_name, possible_id2, position])
                continue

        # step 9: 處理剩餘的mulid
        for test in mulid_test_list2:
            test_name = test[2]
            position = test[4]
            test_near = find_near(test, test_list)
            if test_near != '':
                test_near_name = article[test_near[0]: test_near[1]]
                if current_dict.get(test_near_name):
                    possible_id = current_dict[test_near_name]
                    final_data[position] = possible_id
                else:
                    final_data[position] = test[3]
            else:
                final_data[position] = test[3]

        # step 10: 處理no ans
        # noans_test_list: [test_start, test_end, test_name, '', i]
        for test in noans_test_list:
            test_name = test[2]
            position = test[4]
            test_near = find_near(test, test_list)
            if test_near != '':
                test_near_name = article[test_near[0]: test_near[1]]
                if current_dict.get(test_near_name):
                    possible_id = current_dict[test_near_name]
                    final_data[position] = possible_id
                else:
                    final_data[position] = ''
            else:
                final_data[position] = ''

        return final_data


def main():
    g = GeneNormalizer()
    title = 'TIF1gamma, a novel member of the transcriptional intermediary factor 1 family.'
    article = (
        'We report the cloning and characterization of a novel member of '
        'the Transcriptional Intermediary Factor 1 (TIF1) gene family, '
        'human TIF1gamma. Similar to TIF1alpha and TIF1beta, the structure of '
        'TIF1beta is characterized by multiple domains: RING finger, '
        'B boxes, Coiled coil, PHD/TTC, and bromodomain. Although structurally '
        'related to TIF1alpha and TIF1beta, TIF1gamma presents several '
        'functional differences. In contrast to TIF1alpha, but like TIF1beta, '
        'TIF1 does not interact with nuclear receptors in yeast two-hybrid '
        'or GST pull-down assays and does not interfere with retinoic acid '
        'response in transfected mammalian cells. Whereas TIF1alpha and '
        'TIF1beta were previously found to interact with the KRAB silencing domain '
        'of KOX1 and with the HP1alpha, MODI (HP1beta) and MOD2 (HP1gamma) '
        'heterochromatinic proteins, suggesting that they may participate in a complex '
        'involved in heterochromatin-induced gene repression, '
        'TIF1gamma does not interact with either the KRAB domain of KOX1 or the '
        'HP1 proteins. Nevertheless, TIF1gamma, like TIF1alpha and TIF1beta, '
        'exhibits a strong silencing activity when tethered to a promoter. Since '
        'deletion of a novel motif unique to the three TIF1 proteins, called '
        'TIF1 signature sequence (TSS), abrogates transcriptional repression by '
        'TIF1gamma, this motif likely participates in TIF1 dependent repression.'
    )
    offsets = [
        (0, 9), (211, 220), (233, 242), (247, 255), (274, 282), (415, 424),
        (429, 437), (439, 448), (505, 514), (525, 533), (716, 725), (730, 738),
        (807, 811), (825, 833), (835, 839), (841, 848), (854, 858), (860, 868),
        (1001, 1010), (1060, 1064), (1100, 1109), (1116, 1125), (1130, 1138),
        (1351, 1360)
    ]
    final_data = g.answer(title, article, offsets)
    print(final_data)


if __name__ == '__main__':
    main()
