import pandas as pd


def edit_csv(file_path: str) -> dict:
    """Convert the gene expression data to fold-change dictionary."""
    gene_exp = pd.read_csv(file_path, sep='\t')
    gene_exp.dropna(inplace=True)
    gene_exp.drop(['ID', 'P.Val', 't', 'B', 'Gene.title'], axis=1, inplace=True)

    #  removal of insignificant gene
    fold_change_dict = {}
    for row in gene_exp.itertuples():
        idx, p_val, fc, gene = row
        if p_val < 0.05:
            fold_change_dict[gene] = round(fc, 3)

    return fold_change_dict


