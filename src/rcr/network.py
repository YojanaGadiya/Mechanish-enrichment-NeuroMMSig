import pandas as pd
import networkx as nx

def network_to_file(file_name: str): -> str
    df = pd.read_csv(file_name, sep='\t', header=None)
    df.columns = ['Prt_1', 'Relation', 'Prt_2']

    #change relation values
    relations = {'controls-phosphorylation-of': 'increase', 'controls-state-change-of': 'decrease',
                 'in-complex-with':'increase'}
    for i in relations:
        df.replace(i,relations[i], inplace=True)

    #create new dataframe with changed relationship values.
    new_df = df.loc[df['Relation'].isin(relations.values())]

    #creat network
    DG = nx.DiGraph()
    DG.add_nodes_from(new_df['Prt_1'])
    DG.add_nodes_from(new_df['Prt_2'])
    for _,row in new_df.iterrows():
        DG.add_edge(row['Prt_1'], row['Prt_2'], **{"Relation":row['Relation']})

    path = "test.graphml"
    nx.write_graphml(DG, path)
    return path
