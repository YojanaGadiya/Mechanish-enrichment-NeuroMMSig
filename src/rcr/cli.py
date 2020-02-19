import click

from .reverse_causal_reasoning import rcr_main


@click.command()
@click.option('--file-path', required=True, type=click.Path(dir_okay=False, file_okay=True),
              help='Path to networkx pickle file')
@click.option('--fold-change', required=True, type=dict, help='Dictionary containing the log change fold.')
@click.option('--threshold', required=False, default=2.0, type=float, help='Threshold value for fold-change')
@click.option('--output-path', required=True, type=click.Path(dir_okay=False, file_okay=True),
              help='Path to save TSV file.')
def cli(graph_path: str, fold_change: dict, threshold: int, output_path: str):
    click.echo('Running the task.')  # equivalent to print
    rcr_main(graph_path, fold_change, threshold, output_path)


if __name__ == '__main__':
    cli()
