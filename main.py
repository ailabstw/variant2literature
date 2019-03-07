# pylint: disable=protected-access
"""index worker on azure kubernetes
"""
import argparse
import time
import os
import multiprocessing
import queue
import logging
import traceback

from var_utils import VarNormalizer
import parse_data
import var_ner
import gene_ner
import assign_gene
import normalize_var
from utils import timeout

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """
    Returns:
        arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--n-process", type=int, default=1)
    parser.add_argument("--loglevel", type=str, default='INFO')
    parser.add_argument("--input", type=str, default='/app/input')

    parser.set_defaults(nxml_only=False)
    parser.add_argument("--nxml-only", action='store_true', dest='nxml_only')

    parser.set_defaults(table_detect=True)
    parser.add_argument("--no-table-detect", action='store_false', dest='table_detect')

    args = parser.parse_args()
    return args


@timeout(180)
def extract(_id, idx, data, var_extr, gene_extr):  # pylint: disable=too-many-locals
    """process a paper
    """
    body_var, table_var = var_ner.process(_id, data.body, data.tables, var_extr)
    body_gene, caption_gene, table_gene = gene_ner.process(_id, data.body, data.tables, gene_extr)

    body_results = assign_gene.process_body(idx, data.body, body_gene, body_var)
    table_results = assign_gene.process_table(idx, data.tables, body_gene, caption_gene,
                                              table_gene, table_var)
    results = body_results + table_results
    return results


def worker(que, args):  # pylint: disable=too-many-locals
    """worker for one process
    """
    var_extr = var_ner.pytmvar.Extractor()
    gene_extr = gene_ner.pygnormplus.Extractor()
    var_normalizer = VarNormalizer()

    logger.info('init OK')

    for _ in range(10):
        try:
            msg = que.get(timeout=10)
        except queue.Empty:
            return

        _id, dir_path = msg

        logger.info('start processing %s ...', _id)
        t0 = time.time()

        try:
            parsed_data = parse_data.process(_id, dir_path,
                                             nxml_only=args.nxml_only,
                                             table_detect=args.table_detect,
                                             save_data=False)

            results = []
            for idx, filename, data in parsed_data:
                try:
                    results += extract(_id, idx, data, var_extr, gene_extr)
                except TimeoutError:
                    logger.info(f'timeout {_id} {filename}')

            normalize_var.process(results, _id, var_normalizer)
        except Exception:
            traceback.print_exc()

        logger.info('end processing {}: {:.3f} secs'.format(_id, time.time() - t0))


def main():
    """main function
    """
    args = parse_args()

    if args.loglevel:
        logging.getLogger().setLevel(getattr(logging, args.loglevel))

    que = multiprocessing.Queue()

    for pmid in os.listdir(args.input):
        dir_path = os.path.join(args.input, pmid)
        if os.path.isdir(dir_path):
            pmid = os.path.basename(dir_path)
            que.put((pmid, dir_path))

    ps = []
    for _ in range(args.n_process):
        p = multiprocessing.Process(target=worker, args=(que, args), daemon=True)
        p.start()
        ps.append(p)

    while que.qsize() > 0:
        new_ps = []
        for p in ps:
            if p.is_alive():
                new_ps.append(p)
            else:
                p.terminate()
                p = multiprocessing.Process(target=worker, args=(que, args), daemon=True)
                p.start()
                new_ps.append(p)
        ps = new_ps
        time.sleep(10)

    for p in ps:
        p.join()


if __name__ == '__main__':
    main()
