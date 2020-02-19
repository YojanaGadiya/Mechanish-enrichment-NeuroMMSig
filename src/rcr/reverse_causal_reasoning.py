import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests
from .network import network_to_file
from .dgexp_edit import edit_csv


def overlay(graph: nx.Graph, fold_change_dict: dict, threshold: float) -> nx.Graph:
    """Overlaying the fold-change data onto the graph."""
    for i in graph.nodes():
        # check if it exists
        if i not in graph.nodes():
            print('{} not found in graph', i)
            continue
        else:
            if fold_change_dict[i] > threshold:  # increased expression
                graph.add_node(i, change=1)
            elif fold_change_dict[i] < threshold:  # decreased expression
                graph.add_node(i, change=-1)
            else:  # no change
                graph.add_node(i, change=0)
    return graph


def shortest_path(graph: nx.Graph, hyp_node: str) -> dict:
    """Returning the shortest path of the hype nodes with all the other nodes."""
    return nx.single_source_shortest_path(graph, hyp_node)


def edge_label_value(graph: nx.Graph, path_list: list):
    """Return the product of the edges value of the path."""
    if len(path_list) == 1:
        return 0
    else:
        edge = {'increase': 1, 'decrease': -1}
        edge_list = []

        for i in range(len(path_list) - 1):
            k = graph.edges[path_list[i], path_list[i + 1]]  # edge dictionary attribute

            # check if exist
            if 'Relation' not in k:
                print('Relation value missing for {0} and {1}'.format(path_list[i], path_list[i + 1]))
                raise ValueError('Relation value missing')
            else:
                edge_list.append(edge[k['Relation']])
        return np.prod(edge_list, dtype=np.int8)


def node_label_value(graph: nx.Graph, path_list: list):
    """Returns the product of the starting and end node of the path."""
    if len(path_list) == 1:
        return 0
    else:
        node_list = []
        node_attr_dict = nx.get_node_attributes(graph, 'change')
        for i in node_attr_dict:
            if node_attr_dict[i] == 0:
                print('The ')
            if i in path_list:
                node_list.append(node_attr_dict[i])
        return np.prod(node_list, dtype=np.int8)


def p_value(concordance_count: int, nodes: int, p: float) -> float:
    """Return the p-value.

    @:param concordance_count : the number of successful prediction
    @:param total_nodes : the total number of downstream nodes
    @:param p : probability of achieving a result
    """
    return binom.pmf(concordance_count, nodes, p)


def p_val_correction(p: list) -> list:
    """Uses Benjamini and Hochberg p-value correction."""
    return multipletests(p, alpha=0.05, method='fdr_bh')


def calculate_concordance(graph: nx.Graph, hyp_node: str) -> tuple:
    """Calculation of concordance, non-concordance and p-value for the data."""
    if hyp_node not in graph:
        raise ValueError('Node not preset in graph.')
    else:
        concordance_count = 0
        non_concordance_count = 0
        ambiguous = 0

        path_dict = shortest_path(graph, hyp_node)
        node_num = len(path_dict) - 1  # to remove the node path with itself.

        for i in path_dict:
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


def rcr_main(file_path: str, gene_exp_path: str, threshold: float, output_file: str):
    # getting graph from pathway data.
    path = network_to_file(file_path)
    G = nx.read_graphml(path)

    # getting fold change dictionary from gene expression data.
    fold_change = edit_csv(gene_exp_path)

    overlay_graph = overlay(G, fold_change, threshold)
    concordance_dict = {}

    for i in overlay_graph.nodes():
        concordance = calculate_concordance(overlay_graph, i)
        concordance_dict[i] = concordance

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
    concordance_df.to_csv(output_file, index=True)
    return


if __name__ == "__main__":
    rcr_main()
