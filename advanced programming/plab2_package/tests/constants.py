"""Test module constants."""
from pathlib import Path

TEST_FOLDER = Path(__file__).parent
TEST_DATA_DIR = TEST_FOLDER.joinpath("data")
PPI_FILE_PATH = TEST_DATA_DIR.joinpath("ppis.csv")
NODE_LIST_PATH = TEST_DATA_DIR.joinpath("node_list.tsv")
EDGE_LIST_PATH = TEST_DATA_DIR.joinpath("edge_list.tsv")
IMPORT_PPI_PATH = TEST_DATA_DIR.joinpath("import_ppi.csv")