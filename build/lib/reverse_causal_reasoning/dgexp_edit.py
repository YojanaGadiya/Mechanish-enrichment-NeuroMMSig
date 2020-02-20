from typing import Dict

import pandas as pd


def edit_csv(file_path: str) -> Dict:
    """Convert the gene expression data to fold-change dictionary."""
    if file_name.endswith('.csv'):
        gene_exp = pd.read_csv(file_path, sep='\t')
        gene_exp.dropna(inplace=True)
        gene_exp.drop(['ID', 'P.Value', 't', 'B', 'Gene.title'], axis=1, inplace=True)

        #  return fold change dict
        fold_change_dict = {}
        for row in gene_exp.itertuples():
            idx, _, fc, gene = row
            fold_change_dict[gene.upper()] = round(fc, 3)

        return fold_change_dict
    else:
        raise ValueError('Please input a CSV gene expression file.')