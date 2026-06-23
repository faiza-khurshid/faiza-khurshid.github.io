"""Unit tests for utils.py."""

import os
import pytest

from pathlib import Path
from sqlalchemy import create_engine, inspect, select

from plab2.network import Network
from plab2.startup import DATA_DIR, PROJECT_DIR
from plab2.utils import ApiInterface, Database
from plab2.constants import HGNC, UNIPROT
from plab2.models import Hgnc, Uniprot

from .constants import IMPORT_PPI_PATH

TEST_SYMBOL = "TNF"
TEST_ACCESSION = 'P01375'
test_hgnc_cache_file_path = Path(DATA_DIR, HGNC, TEST_SYMBOL.lower() + ".json")
test_uniprot_cache_file_path = Path(DATA_DIR, UNIPROT, TEST_ACCESSION + ".xml")
interface = ApiInterface(symbol=TEST_SYMBOL)

TEST_DB = Path(PROJECT_DIR, 'test.db')
TEST_DB_CONN = f"sqlite:///{TEST_DB}"
test_engine = create_engine(TEST_DB_CONN)


class TestDatabase:
    """Class comprising of unit tests for Database methods found in utils.py."""

    @pytest.fixture(scope='module')
    def test_db(self):
        """Create test DB and drop after."""
        db = Database(db_engine=test_engine)
        inspector = inspect(test_engine)

        print('Build Test Database')
        db.build_database()
        assert TEST_DB.is_file()  # DB is created

        # Check tables and columns
        tables = inspector.get_table_names()
        hgnc_cols = ("id", "hgnc_id", "ensembl", "symbol")
        uniprot_cols = ("id", "accession", "tax_id", "name", "fullname", "hgnc")
        assert all([x in tables for x in (HGNC, UNIPROT)])  # Correct tables
        assert all([x in Hgnc.__table__.columns for x in hgnc_cols])  # Correct HGNC columns
        assert all([x in Uniprot.__table__.columns for x in uniprot_cols])  # Correct Uniprot columns

        yield db

        print('Delete Test Database')
        Path.unlink(TEST_DB)
        assert not TEST_DB.is_file()

    def test_import(self, test_db):
        """Test that importer works."""
        net = Network(ppi_file=IMPORT_PPI_PATH, enrich=True, db_engine=test_engine)  # Enrich to populate DB

        hgnc_stmt = select(Hgnc.hgnc_id)
        up_stmt = select(Uniprot.accession)
        hgnc_entries = net.session.execute(hgnc_stmt).all()
        up_entries = net.session.execute(up_stmt).all()

        assert len(hgnc_entries) == 3
        assert len(up_entries) == 3

        assert all([hgnc_id[0] in ('HGNC:10298', 'HGNC:6001', 'HGNC:11892') for hgnc_id in hgnc_entries])
        assert all([acc_num[0] in ('P27635', 'P60568', 'P01375') for acc_num in up_entries])

        net.session.close()  # Close connection to DB


class TestApiInterface:
    """Class comprising of unit tests for ApiInterface methods found in utils.py."""

    def test_download_info(self):
        """Unit tests for the download_hgnc_info method."""
        interface.get_identifiers()
        interface.raw_data = dict()
        assert not interface.raw_data

        interface.download_info(HGNC)
        assert interface.raw_data[HGNC]
        assert interface.raw_data[HGNC]['numFound'] == 1

        interface.download_info(UNIPROT)
        assert interface.accession == TEST_ACCESSION
        assert interface.raw_data[UNIPROT]

    def test_cache_file_exists(self):
        """Unit tests for the cache_file_exists method."""
        interface.get_identifiers()
        pairs = ((HGNC, test_hgnc_cache_file_path), (UNIPROT, test_uniprot_cache_file_path))
        for pair in pairs:
            db, cf_path = pair

            if Path.is_file(cf_path):
                Path.unlink(cf_path)

            assert interface.cache_file_exists(cf_path) is False
            interface.read_data(db)
            assert interface.cache_file_exists(cf_path) is True

    def test_get_identifiers(self):
        """Unit tests for the get_identifiers method."""
        ids = interface.get_identifiers()
        expected = {
            'hgnc': {
                'ensembl': 'ENSG00000232810', 'hgnc_id': 'HGNC:11892',  'uniprot': ['P01375'],
            },
            'uniprot': {
                'accession': 'P01375', 'full_name': 'Tumor necrosis factor',
                'name': 'TNFA_HUMAN', 'symbol': 'TNF', 'tax_id': '9606',
            }
        }
        assert ids == expected
