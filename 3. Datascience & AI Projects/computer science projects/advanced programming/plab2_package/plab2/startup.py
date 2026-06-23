import os
import logging

from pathlib import Path
from sqlalchemy import create_engine

from plab2.constants import HGNC, UNIPROT

logger = logging.getLogger(__name__)

# Default Directory Paths
# home_dir = os.path.expanduser("~")  # Also acceptable
home_dir = Path.home()
PROJECT_DIR = home_dir.joinpath(".bschultz")
LOG_DIR = PROJECT_DIR.joinpath("logs")
DATA_DIR = PROJECT_DIR.joinpath("data")
DB_PATH = PROJECT_DIR.joinpath("plab2.db")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(os.path.join(DATA_DIR, HGNC), exist_ok=True)
os.makedirs(os.path.join(DATA_DIR, UNIPROT), exist_ok=True)

# Logging Configuration
LOG_FILE_PATH = os.path.join(LOG_DIR, "bruces_log.log")

fh = logging.FileHandler(LOG_FILE_PATH)
fh.setLevel(logging.DEBUG)

sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
sh.setFormatter(formatter)

logger.addHandler(fh)
logger.addHandler(sh)

# SQLite init
CONN_STRING = f"sqlite:///{DB_PATH}"
engine = create_engine(CONN_STRING)
