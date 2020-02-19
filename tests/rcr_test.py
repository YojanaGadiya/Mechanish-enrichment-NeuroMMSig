import unittest
import os
import networkx as nx
from rcr import reverse_causal_reasoning

class TestRcr(unittest.TestCase):
    """Checks for function with the following graph data:
                -1         +1
             |-----|C(-1)------>F(+1)
    A(+1) ---|                |---------->E(-1)
             |----->B(+1)-----|     +1
                +1            |-------->D(-1)
                                    +1
    """
    def test_network_to_file(self):
        HERE = os.path.dirname(__file__)
        MY_TEST_DATA_PATH_1 = os.path.join(HERE, 'test_data.csv')
        MY_TEST_DATA_PATH_2 = os.path.join(HERE, 'test_data.tsv')

        path_1 = reverse_causal_reasoning.network_to_file(MY_TEST_DATA_PATH_1)
        self.assertEqual('test.graphml',path_1)

    def test_overlay_graph(self):
        fold_change = {'A': 1, 'B': 1, 'C': -1, 'D': -1, 'E': -1, 'F': 1}
        G = nx.read_graphml('test.graphml')
        overlay_graph = reverse_causal_reasoning.overlay(G,fold_change,0)
        node_attribute = nx.get_node_attributes(overlay_graph,'change')
        self.assertEqual(fold_change, node_attribute)

    def test_concordance(self):
        fold_change = {'A': 1, 'B': 1, 'C': -1, 'D': -1, 'E': -1, 'F': 1}
        G = nx.read_graphml('test.graphml')
        overlay_graph = reverse_causal_reasoning.overlay(G, fold_change, 0)

        concordance_dict = {}
        for i in overlay_graph.nodes():
            node_num, concord, non_concord, p_val = reverse_causal_reasoning.calculate_concordance(overlay_graph, i)
            concordance_dict[i] = (node_num, concord, non_concord)

        expected = {'A': (5, 3, 2), 'B': (2, 0, 2), 'C': (1, 0, 1),
                    'D': (0, 0, 0), 'E': (0, 0, 0), 'F': (0, 0, 0)}
        self.assertEqual(expected, concordance_dict)


