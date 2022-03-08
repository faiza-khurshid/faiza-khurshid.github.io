"""Network module tests."""

import os
import pytest
import pandas as pd
import networkx as nx

from pathlib import Path
from collections import Counter

from plab2.startup import DATA_DIR
from plab2.network import Network, Analyzer, Statistics

from .constants import PPI_FILE_PATH, NODE_LIST_PATH, EDGE_LIST_PATH

test_graph_image = Path(DATA_DIR, "graph.pdf")


class TestNetwork:
    """Unit tests for the Network class."""

    def test_import_ppi(self):
        """Checks whether the Network class correctly imports a PPI file."""
        net = Network(ppi_file=PPI_FILE_PATH)
        assert isinstance(net.nodes, dict)
        assert len(net.nodes) == 54
        assert net.graph.number_of_edges() == 38

    def test_import_node_edge_lists(self):
        """Checks whether the Network class correctly imports node and edge lists."""
        net = Network(node_list=NODE_LIST_PATH, edge_list=EDGE_LIST_PATH)
        assert isinstance(net.nodes, dict)
        assert len(net.nodes) == 54
        assert net.graph.number_of_edges() == 38

    def test_enrich_graph(self):
        """Tests for enrich_graph method."""
        net = Network(ppi_file=PPI_FILE_PATH)
        assert len(net.nodes) == 54
        assert net.graph.number_of_edges() == 38  # Without enrichment

        net.enrich_graph()
        assert len(net.nodes) == 54 * 3  # Should be 3 times as many nodes
        molecules = [x['molecule'] for x in net.nodes.values()]
        molecule_counts = Counter(molecules)
        assert molecule_counts['protein'] == 54
        assert molecule_counts['rna'] == 54
        assert molecule_counts['dna'] == 54

        edge_attrs = nx.get_edge_attributes(net.graph, 'rel_type')  # Returns edges with rel_type as dict
        edge_counts = Counter(edge_attrs.values())  # Counts rel_types
        assert edge_counts['translated'] == 54
        assert edge_counts['transcribed'] == 54


class TestAnalyzer:
    """Unit tests for the Analyzer class."""

    def test_generate_graph_image(self):
        """Checks whether the Analyzer class can generate a PDF."""
        if os.path.isfile(test_graph_image):
            os.remove(test_graph_image)
        assert not os.path.isfile(test_graph_image)

        Analyzer(ppi_file=PPI_FILE_PATH).generate_graph_image(graph_output_path=test_graph_image)
        assert os.path.isfile(test_graph_image)

    def test_shortest_paths(self):
        """Checks whether shortest paths can be correctly identified."""
        source, target = "SYMPK", "TRPM7"
        analyze = Analyzer(ppi_file=PPI_FILE_PATH)
        test_symbols = ["SYMPK", "KAT6A", "H3-4", "TRPM7"]

        sps = analyze.shortest_paths(source=source, target=target)
        sp = [analyze.nodes[node_id]['symbol'] for node_id in sps[0]]
        assert analyze.paths is not None
        assert len(sps) == 1
        assert sp == test_symbols


class TestStatistics:
    """Unit tests for the Statistics class."""

    def test_summary_statistics(self):
        """Checks the summary statistics method."""
        stat = Statistics(ppi_file=PPI_FILE_PATH)
        stat.summary_statistics()
        assert isinstance(stat.sum_stats, pd.DataFrame)

        for _, row in stat.sum_stats.iterrows():
            stat = row['Stat']
            val = row['Value']
            if stat == "Nodes":
                assert val == 54
            elif stat == "Graph Density":
                assert val < 0.03
            elif stat == "physical_association":
                assert val == 27
            elif stat == "direct_interaction":
                assert val == 6
            elif stat == "colocalization":
                assert val == 4
            elif stat == "association":
                assert val == 1
            else:
                assert val is not None

    def test_export_stats(self):
        """Checks the export_stats method."""
        stat = Statistics(ppi_file=PPI_FILE_PATH)
        stat.summary_statistics()
        for ext in [".tsv", ".csv", ".json"]:
            export_path = os.path.join(DATA_DIR, "test_export" + ext)
            stat.export_stats(output_path=export_path)
            assert os.path.isfile(export_path)
            os.remove(export_path)
            assert not os.path.isfile(export_path)

        bad_path = os.path.join(DATA_DIR, "test_export.mrbean")
        with pytest.raises(ValueError) as err:
            stat.export_stats(output_path=bad_path)
            assert "Stats output path not csv, tsv, or json!" in str(err.value)
