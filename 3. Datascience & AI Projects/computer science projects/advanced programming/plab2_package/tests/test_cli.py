"""Tests for `plab2` package."""

from pathlib import Path

from click.testing import CliRunner

from plab2.cli import compile
from plab2.startup import PROJECT_DIR

TEST_FOLDER = Path(__file__).parent
PPI_FILE_PATH = TEST_FOLDER.joinpath("data/ppis.csv")


class TestCli:
    """Class for testing the CLI commands."""

    def test_compile(self):
        """Test the compile CLI command."""
        node_path = Path(PROJECT_DIR, "nodes.tsv")
        edge_path = Path(PROJECT_DIR, "edges.tsv")
        runner = CliRunner()
        result = runner.invoke(compile, [str(PPI_FILE_PATH), str(node_path), str(edge_path)])
        assert result.exit_code == 0
        assert node_path.is_file()
        assert edge_path.is_file()
        node_path.unlink()
        edge_path.unlink()
        assert not node_path.is_file()
        assert not edge_path.is_file()

        help_result = runner.invoke(compile, ['--help'])
        assert "Enrich the graph with RNA and DNA molecules." in help_result.output
