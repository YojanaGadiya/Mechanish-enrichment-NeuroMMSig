# -*- coding: utf-8 -*-

"""Command line interface code for reverse causal reasoning (RCR)."""

import click

from .reverse_causal_reasoning import rcr_main


@click.command()
@click.option(
    '--file-path',
    required=True,
    type=click.Path(dir_okay=False, file_okay=True),
    help='Path to networkx pickle file'
)
@click.option(
    '--file-sep',
    required=False,
    default='\t',
    type=str,
    help='Separator for the network file (CSV: ,)'
)
@click.option(
    '--gene-exp-data',
    required=True,
    type=click.Path(dir_okay=False, file_okay=True),
    help='Path for gene expression data.'
)
@click.option(
    '--gene-exp-data-sep',
    required=False,
    default='\t',
    type=str,
    help='Separator for gene expression data file (CSV : ,)'
)
@click.option(
    '--threshold',
    required=False,
    default=0.5,
    type=float,
    help='Threshold value for fold-change'
)
@click.option(
    '--output-path',
    required=True,
    type=click.Path(dir_okay=False,  file_okay=True),
    help='Path to save CSV file')
def cli(
        file_path: str,
        file_sep: str,
        gene_exp_data: str,
        gene_exp_data_sep: str,
        threshold: int,
        output_path: str
):
    click.echo('RCR concordance count initialized')  # equivalent to print
    rcr_main(file_path, file_sep, gene_exp_data, gene_exp_data_sep, threshold, output_path)
