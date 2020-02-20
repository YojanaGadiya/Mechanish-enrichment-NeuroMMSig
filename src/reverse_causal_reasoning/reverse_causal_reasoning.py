# -*- coding: utf-8 -*-

"""Code to carry out reverse causal reasoning (RCR)."""

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests
import os
from tqdm import tqdm
from .network import network_to_file
from .dgexp_edit import edit_csv


def overlay(
        graph: nx.Graph,
        fold_change_dict: dict,
        threshold: float
) -> nx.Graph:
    """Return the overlayed graph with fold-change data.

    :param graph : Graph of the knowledge.
    :param fold_change_dict : The fold-change dictionary from the gene expression data.
    :param threshold : The threshold value to differentiate up-regulated and down-regulated genes from the ambiguous one.
    """
    for i in graph.nodes():
        # check if it exists
        if i.upper() not in fold_change_dict:
            print('{} not found in graph'.format(i))
        else:
            if fold_change_dict[i.upper()] > threshold:  # increased expression
                graph.add_node(i, change=1)
            elif fold_change_dict[i.upper()] < -threshold:  # decreased expression
                graph.add_node(i, change=-1)
            else:  # no change
                graph.add_node(i, change=0)
    return graph


def shortest_path(
        graph: nx.Graph,
        hyp_node: str
) -> dict:
    """Return the shortest path of the hype nodes with all the other nodes."""
    return nx.single_source_shortest_path(graph, hyp_node)


def edge_label_value(
        graph: nx.Graph,
        path_list: list
):
    """Return the product of the edges value of the path.

    :param graph : The graph of the knowledge created.
    :param path_list : The
    """
    if len(path_list) == 1:
        return 0

    edge = {'increase': 1, 'decrease': -1}
    edge_list = []

    for i in range(len(path_list) - 1):
        k = graph.edges[path_list[i], path_list[i + 1]]  # edge dictionary attribute

        # check if exist
        if 'Relation' not in k:
            raise ValueError('Relation value missing for {} and {}'.format(path_list[i], path_list[i + 1]))
        edge_list.append(edge[k['Relation']])

    return np.prod(edge_list, dtype=np.int32)


def node_label_value(
        graph: nx.Graph,
        path_list: list
):
    """Return the product of the starting and end node of the path.

    :param graph : The graph of the knowledge created.
    :param path_list : The list of the path of one node to another.
    """
    if len(path_list) == 1:
        return 0
    node_list = []
    node_attr_dict = nx.get_node_attributes(graph, 'change')

    # calculating product of first and last node of graph.
    start_node = path_list[0]
    end_node = path_list[-1]

    node_list.append(node_attr_dict[start_node])
    node_list.append(node_attr_dict[end_node])

    return np.prod(node_list, dtype=np.int8)


def p_value(
        concordance_count: int,
        nodes: int,
        p: float
) -> float:
    """Return the p-value.

    :param concordance_count : The number of successful prediction.
    :param nodes : The total number of downstream nodes excluding the ambiguous one.
    :param p : Probability of achieving a result.
    """
    return binom.pmf(concordance_count, nodes, p)


def p_val_correction(
        p: list
) -> list:
    """Return corrected p-value using Benjamini and Hochberg correction

    :param p : List of all p-values.
    """

    return multipletests(p, alpha=0.05, method='fdr_bh')


def calculate_concordance(
        graph: nx.Graph,
        hyp_node: str
) -> tuple:
    """Calculation of concordance, non-concordance and p-value for the data

    :param graph : Graph to be used for calculating the concordance.
    :param hyp_node : Starting point of the graph for calculation.
    """
    if hyp_node not in graph:
        raise ValueError('Node not preset in graph.')
    concordance_count = 0
    non_concordance_count = 0
    ambiguous = 0

    path_dict = shortest_path(graph, hyp_node)
    node_num = len(path_dict) - 1  # to remove the node path with itself.

    for i in tqdm(path_dict):
        path = path_dict[i]  # path to travel
        edge_val = edge_label_value(graph, path)  # edge product value
        node_val = node_label_value(graph, path)  # node product value

        # concordance and non-concordance count
        if edge_val != 0 or node_val != 0:
            if edge_val == node_val:
                concordance_count += 1
            else:
                non_concordance_count += 1
        elif node_val == 0:
            ambiguous += 1

    p_val = p_value(concordance_count, node_num - ambiguous, 0.5)
    return node_num, concordance_count, non_concordance_count, p_val


def rcr_main(
        file_path: str,
        file_sep: str,
        gene_exp_path: str,
        gene_file_sep: str,
        threshold: float,
        output_file: str
):
    """Main method for RCR algorithm.

    :param file_path: The knowledge file location used for graph creation.
    :param file_sep: The separator the knowledge file.
    :param gene_exp_path: The gene expression data.
    :param gene_file_sep: The separator for the gene expression data file.
    :param threshold: The fold-change threshold.
    :param output_file: The name of the output file.
    :return:
    """
    # getting graph from pathway data.
    pathway = network_to_file(file_path, file_sep)

    # getting fold change dictionary from gene expression data.
    fold_change = edit_csv(gene_exp_path, gene_file_sep)

    print('Calculating convergence..')
    overlay_graph = overlay(pathway, fold_change, threshold)
    concordance_dict = {}

    for i in overlay_graph.nodes():
        concordance_dict[i] = calculate_concordance(overlay_graph, i)

    # remove leaf nodes or nodes with no downstream nodes.
    remove = []
    for i in concordance_dict:
        node_num, _, _, _ = concordance_dict[i]
        if node_num == 0.0:
            remove.append(i)

    for node in remove:
        del concordance_dict[node]

    # forming dataframe
    concordance_df = pd.DataFrame.from_dict(concordance_dict)
    concordance_df = concordance_df.transpose()
    concordance_df.columns = ['No_of_Nodes', 'Concordance', 'Non-concordance', 'p-value']

    for i in ['No_of_Nodes', 'Concordance', 'Non-concordance']:
        concordance_df[i] = pd.to_numeric(concordance_df[i])

    # p-value correction
    p_val = list(concordance_df['p-value'])
    corrected_p_val = p_val_correction(p_val)[1]
    concordance_df['Corrected p-val'] = corrected_p_val

    # output file saving
    dir = os.path.dirname(__file__)
    output_path = os.path.join(dir, output_file)
    print('Saving to file {}'.format(output_path))
    concordance_df.to_csv(output_file, index=True)


if __name__ == "__main__":
    rcr_main()
