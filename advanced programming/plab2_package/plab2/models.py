from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import declarative_base

Base = declarative_base()

# Define tables using class definitions
class Hgnc(Base):
    __tablename__ = 'hgnc'

    id = Column(Integer, primary_key=True, autoincrement=True)
    hgnc_id = Column(String)
    ensembl = Column(String)
    symbol = Column(String)

class Uniprot(Base):
    __tablename__ = 'uniprot'

    id = Column(Integer, primary_key=True, autoincrement=True)
    accession = Column(String)
    tax_id = Column(Integer)
    name = Column(String)
    fullname = Column(String)
    hgnc = Column(Integer, ForeignKey(Hgnc.id))
