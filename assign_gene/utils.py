"""utils for assigning gene
"""
from typing import NamedTuple


class Result(NamedTuple):
    """ner result
    """
    variant: str
    gene: str
    var_json: str
    gene_id: int
    file_idx: int
    in_table: bool
    table_idx: int
    row: int
    col: int
    start: int
    end: int
    gene_start: int
    gene_end: int
