import click

from .reverse_causal_reasoning import rcr_main


@click.command()
@click.option('--file-path', required=True, type=click.Path(dir_okay=False, file_okay=True),
              help='Path to networkx pickle file')
@click.option('--gene-exp-data', required=True, type=click.Path(dir_okay=False, file_okay=True),
              help='Path for gene expression data.')
@click.option('--threshold', required=False, default=2.0, type=float,
              help='Threshold value for fold-change.')
@click.option('--output-path', required=True, type=click.Path(dir_okay=False, file_okay=True),
              help='Path to save TSV file.')
def cli(graph_path: str, gene_exp_data: str, threshold: int, output_path: str):
    click.echo('Running the task.')  # equivalent to print
    rcr_main(graph_path, gene_exp_data, threshold, output_path)


if __name__ == '__main__':
    cli()
