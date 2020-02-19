import pandas as pd
import networkx as nx
import random


def network_to_file(file_name: str) -> str:
    if file_name.endswith('.csv'):
        df = pd.read_csv(file_name, header=None)
        df.columns = ['Prt_1', 'Relation', 'Prt_2']
        # changing relation values
        for relation_value in df['Relation']:
            if relation_value not in ['increase', 'decrease']:
                df.replace(relation_value, random.choice(['increase', 'decrease']), inplace=True)

        # create network
        DG = nx.DiGraph()
        DG.add_nodes_from(df['Prt_1'])
        DG.add_nodes_from(df['Prt_2'])
        for _, row in df.iterrows():
            DG.add_edge(row['Prt_1'], row['Prt_2'], **{"Relation": row['Relation']})

        path = "test.graphml"
        nx.write_graphml(DG, path)
        return path
    else:
        raise ValueError("Please pass a CSV file without header.")
