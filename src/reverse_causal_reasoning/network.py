# -*- coding: utf-8 -*-

"""Code to create graph from given file."""

import pandas as pd
import networkx as nx
import random
from .constants import *


def network_to_file(file_name: str, delimiter: str) -> nx.Graph:
    """Create graph from triplets (Protein 1, Interaction, Protein 2)

    :param file_name: Path of file where the triplets are stored
    :param delimiter: Separator used in the file (Default - tab)
    """
    df = pd.read_csv(file_name, sep=delimiter, header=None)
    df.columns = NETWORKCOL # ['Protein1 , Relation, Protein2]

    # changing relation values
    for relation_value in df[NETWORKCOL[1]]:
        if relation_value not in RELATIONVAL:
            df.replace(relation_value, random.choice(RELATIONVAL), inplace=True)

    # create directed acyclic network
    network = nx.DiGraph()
    network.add_nodes_from(df[NETWORKCOL[0]])
    network.add_nodes_from(df[NETWORKCOL[2]])
    for _, row in df.iterrows():
        network.add_edge(row[NETWORKCOL[0]], row[NETWORKCOL[2]], **{RELATION: row[NETWORKCOL[1]]})

    return network
