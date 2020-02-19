import pandas as pd


def edit_csv(file_path: str) -> dict:
    gene_exp = pd.read_csv(file_path, header=None)
    gene_exp.columns = ['Gene_1', 'Fold_change', 'p_val']


