import json
import time
import logging
import requests
import xmltodict

from pathlib import Path
from typing import Optional, Union
from sqlalchemy import select, inspect
from sqlalchemy.orm import Session
from sqlalchemy_utils import database_exists

from plab2.startup import DATA_DIR, engine
from plab2.models import Base, Hgnc, Uniprot
from plab2.constants import HGNC, UNIPROT, HGNC_API, UNIPROT_API

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Database:
    """For interfacing with the relational database."""

    def __init__(self, db_engine=engine):
        self.engine = db_engine
        self.session = Session(bind=self.engine)

        # When DB doesn't exist
        if not database_exists(self.engine.url):
            self.build_database()

        # When DB exists, but tables not made yet
        tables = inspect(self.engine).get_table_names()
        if not tables:
            self.build_database()

    def get_identifiers(self, symbol: str) -> Optional[dict]:
        """Get identifiers for a given HGNC symbol from the relational database."""
        stmt = select(
            Hgnc.hgnc_id,
            Hgnc.ensembl,
            Uniprot.accession,
            Uniprot.tax_id
        ).filter_by(symbol=symbol).join(Uniprot)

        entries = self.session.execute(stmt).first()

        if entries:
            return dict(entries)

        else:
            profiler = ApiInterface(symbol=symbol)
            identifiers = profiler.get_identifiers()
            if identifiers:
                self.add_data_to_db(symbol=symbol, identifiers=identifiers)

            return identifiers

    def rebuild_database(self) -> None:
        """Burn everything and builds the database."""
        self.drop_database()
        self.build_database()

    def build_database(self) -> None:
        """Build the tables of the database."""
        logger.info("Building database...")
        Base.metadata.create_all(bind=self.engine)

    def drop_database(self) -> None:
        """Drop all of the associated tables in the database."""
        logger.warning("Dropping database...")
        Base.metadata.drop_all(bind=self.engine)

    def add_data_to_db(self, symbol: str, identifiers: dict) -> None:
        """Adds new entries to database."""
        # Add HGNC entry to get ID
        hgnc_data = identifiers[HGNC]
        hgnc_entry = Hgnc(hgnc_id=hgnc_data["hgnc_id"], ensembl=hgnc_data["ensembl"], symbol=symbol)
        self.session.add(hgnc_entry)
        self.session.commit()
        logger.info(f"Added {symbol} to HGNC table")
        hgnc_table_id = hgnc_entry.id

        # Add UniProt if there is data
        if UNIPROT in identifiers:
            up_data = identifiers[UNIPROT]
            up_entry = Uniprot(
                accession=up_data["accession"],
                name=up_data["name"],
                fullname=up_data["full_name"],
                tax_id=up_data["tax_id"],
                hgnc=hgnc_table_id
            )
            self.session.add(up_entry)
            self.session.commit()
            logger.info(f"Added {up_data['accession']} to UniProt table")


class ApiInterface:
    """Download data for a given HGNC symbol from HGNC and UniProt and import into database."""

    def __init__(self, symbol: str):
        """Init method.

        Parameters
        ----------
        symbol : str
            HGNC symbol of interest.
        """
        self.symbol = symbol
        self.accession = None
        self.raw_data = dict()

    @staticmethod
    def cache_file_exists(path: Union[str, Path]) -> bool:
        """Checks if cache file exists."""
        return path.is_file() if isinstance(path, Path) else Path.is_file(path)

    def get_identifiers(self) -> dict:
        """Gathers HGNC ID, EnSembl Gene ID, and UniProt IDs from raw HGNC data.

        Returns
        -------
        dict, None
            Returns identifiers as dict if data available, else None.
        """
        self.read_data(database=HGNC)
        identifiers = dict()
        if self.raw_data:
            if HGNC in self.raw_data and self.raw_data[HGNC] and self.raw_data[HGNC]['numFound'] > 0:
                identifiers[HGNC] = self.__parse_hgnc_content()

                # Now that uniprot accession IDs extracted, get first and collect UniProt data if accession # present
                if identifiers[HGNC][UNIPROT]:
                    self.accession = identifiers[HGNC][UNIPROT][0]
                    self.read_data(database=UNIPROT)

            if UNIPROT in self.raw_data and self.raw_data[UNIPROT]:
                identifiers[UNIPROT] = self.__parse_uniprot_content()

        return identifiers

    def __parse_uniprot_content(self) -> dict:
        """Parse the raw UniProt content and returns relevant identifiers."""
        up_response = self.raw_data[UNIPROT]
        entry_dict = xmltodict.parse(up_response)['uniprot']['entry']

        symbol = None
        primary_gene = entry_dict["gene"][0] if isinstance(entry_dict["gene"], list) else entry_dict["gene"]
        gene_name_entry = primary_gene["name"]
        gene_names = gene_name_entry if isinstance(gene_name_entry, list) else [gene_name_entry]
        for gene_name in gene_names:
            if gene_name["@type"] == "primary":
                symbol = gene_name["#text"]
                break

        # Handle full name issues
        full_name_entry = entry_dict["protein"]["recommendedName"]["fullName"]
        if not isinstance(full_name_entry, str):
            full_name_entry = full_name_entry["#text"]

        up_identifiers = {
            'accession': entry_dict["accession"][0],
            'name': entry_dict["name"],
            'full_name': full_name_entry,
            'tax_id': entry_dict["organism"]["dbReference"]["@id"],
            "symbol": symbol,
        }

        if up_identifiers['accession'] != self.accession:
            logger.debug(f"Queried {self.accession} but parsed {up_identifiers['accession']}")

        return up_identifiers

    def __parse_hgnc_content(self) -> dict:
        """Parse the raw HGNC content and returns relevant identifiers."""
        raw_hgnc_data = self.raw_data[HGNC]
        if len(raw_hgnc_data['docs']) > 1:
            logger.debug(f"{self.symbol} has more than 1 'docs'")
        content = raw_hgnc_data['docs'][0]
        hgnc_identifiers = {
            'hgnc_id': content["hgnc_id"] if "hgnc_id" in content else None,
            'ensembl': content["ensembl_gene_id"] if "ensembl_gene_id" in content else None,
            UNIPROT: content["uniprot_ids"] if "uniprot_ids" in content else None
        }
        return hgnc_identifiers

    def read_data(self, database: str):
        """Retrieve the raw data for the given database either by reading from cached file or downloading it."""
        identifier = self.symbol.lower() if database == HGNC else self.accession
        extension = "json" if database == HGNC else "xml"
        cache_path = Path(DATA_DIR, database, f"{identifier}.{extension}")
        if self.cache_file_exists(cache_path):
            with open(cache_path, 'r') as cachefile:
                if database == HGNC:
                    content = json.load(cachefile)
                else:  # Uniprot
                    content = cachefile.read()
            self.raw_data[database] = content

        else:
            downloaded = self.download_info(database=database)
            if downloaded:
                time.sleep(0.1)
                with open(cache_path, 'w') as cachefile:
                    if database == HGNC:
                        json.dump(self.raw_data[HGNC], cachefile, indent=2)
                    else:  # Uniprot
                        cachefile.write(self.raw_data[UNIPROT])
                logger.info(f"Cache file for {identifier} successfully written to {cache_path}")

    def download_info(self, database: str) -> bool:
        """Downloads data using HGNC or UniProt API for a given HGNC symbol.

        Parameters
        ----------
        database : str
            Either "uniprot" or "hgnc"

        Returns
        -------
        bool
            True if data was downloaded, False if not
        """
        if database == HGNC:
            api_query = HGNC_API.format(self.symbol)
            header = {"Accept": "application/json"}

        elif database == UNIPROT:
            api_query = UNIPROT_API.format(self.accession)
            header = {"Accept": "text/xml"}

        else:
            raise ValueError("database must but either 'hgnc' or 'uniprot'!")

        r = requests.get(api_query, headers=header)

        if r.ok is False:
            logger.error(f"{api_query} returned bad status code: {r.status_code}")
            return False

        response = r.json()['response'] if database == HGNC else r.text  # Other key is responseHeader
        if response:  # Found a result
            self.raw_data[database] = response
            return True

        else:
            logger.warning(f"No results found for {self.symbol}")
            return False
