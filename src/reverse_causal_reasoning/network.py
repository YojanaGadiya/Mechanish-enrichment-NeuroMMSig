# -*- coding: utf-8 -*-

"""Code to create graph from given file."""

import pandas as pd
import networkx as nx
import random


def network_to_file(file_name: str, delimiter: str = ',') -> nx.Graph:
    df = pd.read_csv(file_name,sep=delimiter,header=None)
    df.columns = ['Prt_1', 'Relation', 'Prt_2']

    # changing relation values
    for relation_value in df['Relation']:
        if relation_value not in ['increase', 'decrease']:
            df.replace(relation_value, random.choice(['increase', 'decrease']), inplace=True)

    # create directed acyclic network
    network = nx.DiGraph()
    network.add_nodes_from(df['Prt_1'])
    network.add_nodes_from(df['Prt_2'])
    for _, row in df.iterrows():
        network.add_edge(row['Prt_1'], row['Prt_2'], **{"Relation": row['Relation']})

    return network

