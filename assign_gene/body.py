"""assign genes to variants in body text
"""
import bisect
import logging

from .utils import Result

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_nearest_gene(genes, start, end):
    """return nearest gene
    """
    if not genes:
        return 0, '', -1, -1

    mean_pos = (int(start) + int(end)) / 2
    idx = bisect.bisect_left(genes, (mean_pos, 0, ''))
    if idx == 0:
        return genes[0][1:]
    elif idx == len(genes):
        return genes[-1][1:]
    if mean_pos - genes[idx - 1][0] < genes[idx][0] - mean_pos:
        return genes[idx - 1][1:]
    return genes[idx][1:]


def get_sent_idxes(body):
    """get indexes of sentences
    """
    ret, offset = [0], 0
    for sent in body.split('\n'):
        offset += len(sent) + 1
        ret.append(offset)
    return ret


def check_in_same_sent(idxes, start, gene_start):
    """check if the variant and the gene are in the same sentence
    """
    var_index = bisect.bisect_right(idxes, start)
    gene_index = bisect.bisect_right(idxes, gene_start)
    return var_index == gene_index


def process_body(idx, body, body_gene, body_var):  # pylint: disable=too-many-locals
    """assign genes to variants in body text
    """
    sent_idxes = get_sent_idxes(body)

    genes = []
    for mention, (start, end), gene_id in body_gene:
        pos = (int(start) + int(end)) / 2
        genes.append((pos, int(gene_id), mention, start, end))
    genes.sort()

    results = []
    for var_mention, (start, end), var_json in body_var:
        gene_id, gene_mention, gene_start, gene_end = get_nearest_gene(genes, start, end)
        if not check_in_same_sent(sent_idxes, start, gene_start):
            gene_start, gene_end = -1, -1

        result = Result(variant=var_mention, gene=gene_mention,
                        var_json=var_json, gene_id=gene_id,
                        file_idx=idx, in_table=False,
                        table_idx=-1, row=-1, col=-1,
                        start=start, end=end,
                        gene_start=gene_start, gene_end=gene_end)
        results.append(result)
    return results
