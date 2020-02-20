# -*- coding: utf-8 -*-

"""Code to get the fold change dictionary from given GEO expression dataset."""

from typing import Dict

import pandas as pd
from .constants import *


def edit_csv(file_path: str, delimiter: str) -> Dict:
    """Convert the gene expression data to fold-change dictionary."""
    gene_exp = pd.read_csv(file_path, sep=delimiter)
    gene_exp.dropna(inplace=True)
    gene_exp.drop(GENEEXPCOL, axis=1, inplace=True)

    #  return fold change dict
    fold_change_dict = {}
    for row in gene_exp.itertuples():
        idx, _, fc, gene = row
        fold_change_dict[gene.upper()] = round(fc, 3)

    return fold_change_dict
