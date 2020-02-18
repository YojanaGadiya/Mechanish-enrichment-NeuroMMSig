import click

from .reverse_causal_reasoning import rcr_main

@click.command()
@click.option('--file-path', required=True, type=click.Path(dir_okay=False, file_okay=True), help='Path to networkx pickle file')
@click.option('--fold-change', required=True, type=dict, help='Dictionary containing the log change fold.')
@click.option('--output-path', required=True, type=click.Path(dir_okay=False, file_okay=True), help='Path to save TSV file.')
def cli(graph_path: str, fold_change: dict, output_path: str):
    rcr_main(graph_path, fold_change, output_path)


if __name__ == '__main__':
    cli()