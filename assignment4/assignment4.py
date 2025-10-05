#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO   # to read bio datasets
import sqlalchemy
from sqlalchemy import create_engine, text 

# as same as the given snippet code to import the correct version
if sqlalchemy.__version__.startswith('1.4'):
    from sqlalchemy import create_engine
    from sqlalchemy.engine import make_url, URL
    from sqlalchemy.sql import text
elif sqlalchemy.__version__.startswith('2'):
    from sqlalchemy import create_engine
    from sqlalchemy.url import make_url, URL
    from sqlalchemy.sql import text

# as for the assignment mentioned, I read it in hardcode way, otherwise I can read it from CLI
file_path = '/homes/aakhlaghi/archaea.1.genomic.gbff'

'''(optional:CLI)
import argparse
def parse_args():
    p = argparse.ArgumentParser(description="Parse GenBank files and store basic stats in MariaDB/MySQL.")
    p.add_argument("files", nargs="+", help="Path(s) to .gb or .gbff file(s)")  # nargs to add one or multiple files
    return p.parse_args()
'''

engine = create_engine("mysql+mysqldb:///?read_default_file=/homes/aakhlaghi/.my.cnf")

with engine.connect() as conn:
    result = conn.execute(text("SELECT DATABASE();"))
    print("Connected successfully to:", result.scalar())

ddl_statements = [
    """CREATE TABLE IF NOT EXISTS genomes (
         acce_ver VARCHAR(40) NOT NULL,
         organism VARCHAR(255) NOT NULL,
         gene_size BIGINT UNSIGNED NOT NULL,
         gene_count INT UNSIGNED NOT NULL,
         protein_count INT UNSIGNED NOT NULL,
         taxon_id INT UNSIGNED NULL,
         PRIMARY KEY (acce_ver)
       ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;""",

    """CREATE TABLE IF NOT EXISTS proteins (
         id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
         acce_ver VARCHAR(40) NOT NULL,
         prot_id VARCHAR(64),
         product TEXT,
         locus_tag VARCHAR(64),
         gene_name VARCHAR(128),
         loc_start INT UNSIGNED,
         loc_end INT UNSIGNED,
         strand ENUM('+','-'),
         PRIMARY KEY (id),
         FOREIGN KEY (acce_ver) REFERENCES genomes(acce_ver)
           ON DELETE CASCADE ON UPDATE CASCADE
       ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;""",

    """CREATE TABLE IF NOT EXISTS protein_ec (
         prot_id_ref BIGINT UNSIGNED NOT NULL,
         ec VARCHAR(32) NOT NULL,
         PRIMARY KEY (prot_id_ref, ec),
         FOREIGN KEY (prot_id_ref) REFERENCES proteins(id)
           ON DELETE CASCADE ON UPDATE CASCADE
       ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;""",

    """CREATE TABLE IF NOT EXISTS protein_go (
         prot_id_ref BIGINT UNSIGNED NOT NULL,
         go_id VARCHAR(32) NOT NULL,
         PRIMARY KEY (prot_id_ref, go_id),
         FOREIGN KEY (prot_id_ref) REFERENCES proteins(id)
           ON DELETE CASCADE ON UPDATE CASCADE
       ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;"""
]

with engine.begin() as conn:
    for stmt in ddl_statements:
        conn.execute(text(stmt))

print("All tables created successfully.")