from collections import Counter, defaultdict
import logging

from .utils import Result

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def process_table(file_idx, tables, body_gene, caption_gene, table_gene, table_var):
    """assign genes to variants in tables
    """
    gene_cnt, gene_name_dict = Counter(), dict()
    for mention, _, gene_id in body_gene:
        gene_cnt[gene_id] += 1
        gene_name_dict[gene_id] = mention
    if not gene_cnt:
        max_body_gene_id = 0
        max_body_gene_name = ''
    else:
        max_body_gene_id = max(gene_cnt.keys(), key=gene_cnt.get)
        max_body_gene_name = gene_name_dict[max_body_gene_id]

    caption_gene_dict = defaultdict(list)
    for mention, table_idx, (start, end), gene_id in caption_gene:
        caption_gene_dict[table_idx].append((gene_id, mention))

    gene_dict = dict()
    for mention, (table_idx, row, col), (start, end), gene_id in table_gene:
        gene_dict[(table_idx, row, col)] = (gene_id, mention)

    var_dict = defaultdict(list)
    for mention, (table_idx, row, col), (start, end), var_json in table_var:
        var_dict[(table_idx, row, col)].append((var_json, mention, start, end))

    results = []
    for table_idx, table in enumerate(tables):
        if not table['cells']:
            continue

        if len(caption_gene_dict[table_idx]) == 1:
            caption_gene_id, caption_gene_name = caption_gene_dict[table_idx][0]
        else:
            caption_gene_id, caption_gene_name = 0, ''

        n, m = len(table['cells']), len(table['cells'][0])
        prev = [(0, '') for _ in range(m)]

        for i in range(n):
            current = [(0, '') for _ in range(m)]
            for j in range(m):
                gene_id, gene = gene_dict.get((table_idx, i, j), (None, ''))
                if gene_id:
                    current[j] = (gene_id, gene)
                elif j > 0:
                    current[j] = current[j - 1]

            for j in range(m):
                if current[j][0] == 0:
                    current[j] = prev[j]

            for j in range(m):
                for var_json, var, start, end in var_dict[(table_idx, i, j)]:
                    gene_id, gene = current[j]
                    # reason = 'table'
                    if gene_id == 0:
                        if caption_gene_id != 0:
                            gene_id = caption_gene_id
                            gene = caption_gene_name
                            # reason = 'caption'
                        elif max_body_gene_id != 0:
                            gene_id = max_body_gene_id
                            gene = max_body_gene_name
                            # reason = 'body'

                    result = Result(variant=var, gene=gene,
                                    var_json=var_json, gene_id=gene_id,
                                    file_idx=file_idx, in_table=True,
                                    table_idx=table_idx, row=i, col=j,
                                    start=start, end=end,
                                    gene_start=-1, gene_end=-1)
                    results.append(result)
            prev = current

    return results
