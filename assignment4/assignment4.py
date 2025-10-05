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

import time

print("Starting data insertion with executemany...")

# -------------------------------
# تابع استخراج داده‌ها از رکوردها
# -------------------------------
def extract_record_fields(record):
    """Extract genome-level information from GenBank record."""
    acce = record.annotations.get("accessions", [str(record.id)])[0]
    version = record.annotations.get("sequence_version")
    acce_ver = f"{acce}.{version}" if version else acce
    organism = record.annotations.get("organism", "unknown")
    gene_size = len(record.seq)

    gene_count = 0
    protein_count = 0
    for feat in record.features:
        if feat.type == "gene":
            gene_count += 1
        elif feat.type == "CDS":
            protein_count += 1

    taxon_id = None
    for feat in record.features:
        if feat.type == "source":
            for x in feat.qualifiers.get("db_xref", []):
                if x.startswith("taxon:"):
                    taxon_id = int(x.split(":")[1])

    return {
        "acce_ver": acce_ver,
        "organism": organism,
        "gene_size": gene_size,
        "gene_count": gene_count,
        "protein_count": protein_count,
        "taxon_id": taxon_id
    }


def extract_protein_data(record, acce_ver):
    """Extract CDS features (proteins, ECs, GOs)."""
    proteins = []
    ec_map = []
    go_map = []

    for feat in record.features:
        if feat.type == "CDS":
            q = feat.qualifiers
            prot_id = q.get("protein_id", [None])[0]
            product = q.get("product", [None])[0]
            locus_tag = q.get("locus_tag", [None])[0]
            gene_name = q.get("gene", [None])[0]
            loc_start = int(feat.location.start)
            loc_end = int(feat.location.end)
            strand = "+" if feat.location.strand == 1 else "-"

            proteins.append({
                "acce_ver": acce_ver,
                "prot_id": prot_id,
                "product": product,
                "locus_tag": locus_tag,
                "gene_name": gene_name,
                "loc_start": loc_start,
                "loc_end": loc_end,
                "strand": strand
            })

            # EC numbers
            for ec in q.get("EC_number", []):
                ec_map.append((prot_id, ec))
            # GO terms
            for go in q.get("GO_function", []) + q.get("GO", []):
                go_map.append((prot_id, go))

    return proteins, ec_map, go_map

# batch size to avoid limitation
BATCH_SIZE = 200

records = list(SeqIO.parse(file_path, "genbank"))
total = len(records)
print(f"Total records found: {total}")

start_time = time.time()
count = 0

for i in range(0, total, BATCH_SIZE):
    batch = records[i:i + BATCH_SIZE]

    with engine.begin() as conn:
        genome_rows = []
        protein_rows = []
        ec_rows = []
        go_rows = []

        # preparation for each batch 
        for record in batch:
            g = extract_record_fields(record)
            genome_rows.append(g)

            proteins, ecs, gos = extract_protein_data(record, g["acce_ver"])
            for p in proteins:
                protein_rows.append(p)
            for ec in ecs:
                ec_rows.append(ec)
            for go in gos:
                go_rows.append(go)

        # executemany
        if genome_rows:
            conn.execute(text("""
                INSERT INTO genomes (acce_ver, organism, gene_size, gene_count, protein_count, taxon_id)
                VALUES (:acce_ver, :organism, :gene_size, :gene_count, :protein_count, :taxon_id)
                ON DUPLICATE KEY UPDATE
                    organism=VALUES(organism),
                    gene_size=VALUES(gene_size),
                    gene_count=VALUES(gene_count),
                    protein_count=VALUES(protein_count),
                    taxon_id=VALUES(taxon_id);
            """), genome_rows)

        if protein_rows:
            conn.execute(text("""
                INSERT INTO proteins (acce_ver, prot_id, product, locus_tag, gene_name, loc_start, loc_end, strand)
                VALUES (:acce_ver, :prot_id, :product, :locus_tag, :gene_name, :loc_start, :loc_end, :strand);
            """), protein_rows)

        count += len(batch)
        print(f"Committed {count}/{total} records ({time.time() - start_time:.1f}s elapsed)")

print(f"All {count} records inserted successfully in {time.time() - start_time:.1f}s.")

