#!/usr/bin/env python3
import time
from Bio import SeqIO
from sqlalchemy import create_engine, text

# FILE_PATH = "/path/to/input.gbff" 

def extract_genome(record):
    acce = record.annotations.get("accessions", [str(record.id)])[0]
    version = record.annotations.get("sequence_version")
    acce_ver = f"{acce}.{version}" if version else acce

    organism = record.annotations.get("organism", "unknown")
    genome_size = len(record.seq)

    gene_count = sum(1 for f in record.features if f.type == "gene")
    protein_count = sum(1 for f in record.features if f.type == "CDS")

    taxon_id = None
    for f in record.features:
        if f.type == "source":
            for x in f.qualifiers.get("db_xref", []):
                if x.startswith("taxon:"):
                    try:
                        taxon_id = int(x.split(":", 1)[1])
                    except ValueError:
                        taxon_id = None
            break

    return {
        "acce_ver": acce_ver,
        "organism": organism,
        "genome_size": genome_size,
        "gene_count": gene_count,
        "protein_count": protein_count,
        "taxon_id": taxon_id,
    }

def extract_proteins(record, acce_ver):
    proteins = []
    ecs = []
    gos = []

    for feat in record.features:
        if feat.type != "CDS":
            continue

        q = feat.qualifiers
        prot_id = q.get("protein_id", [None])[0]
        product = q.get("product", [None])[0]
        locus_tag = q.get("locus_tag", [None])[0]
        gene_ref = q.get("gene", [None])[0]

        loc_start = int(feat.location.start)  # 0-based
        loc_end = int(feat.location.end)      # end-exclusive
        strand = "+" if feat.location.strand == 1 else "-"

        # fallback key if protein_id missing
        if not prot_id:
            prot_id = f"{acce_ver}:{locus_tag}:{loc_start}-{loc_end}"

        proteins.append({
            "prot_id": prot_id,
            "acce_ver": acce_ver,
            "product": product,
            "locus_tag": locus_tag,
            "gene_ref": gene_ref,
            "loc_start": loc_start,
            "loc_end": loc_end,
            "strand": strand,
        })

        for ec in q.get("EC_number", []):
            ecs.append({"prot_id": prot_id, "ec": ec})

        go_terms = []
        go_terms += q.get("GO", [])
        go_terms += q.get("GO_function", [])
        for go in go_terms:
            gos.append({"prot_id": prot_id, "go_id": go})

    return proteins, ecs, gos

def create_tables(conn):
    conn.execute(text("""
    CREATE TABLE IF NOT EXISTS genomes (
        acce_ver VARCHAR(40) NOT NULL,
        organism VARCHAR(255) NOT NULL,
        genome_size BIGINT UNSIGNED NOT NULL,
        gene_count INT UNSIGNED NOT NULL,
        protein_count INT UNSIGNED NOT NULL,
        taxon_id INT UNSIGNED NULL,
        PRIMARY KEY (acce_ver)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
    """))

    conn.execute(text("""
    CREATE TABLE IF NOT EXISTS proteins (
        prot_id VARCHAR(80) NOT NULL,
        acce_ver VARCHAR(40) NOT NULL,
        product TEXT,
        locus_tag VARCHAR(80),
        gene_ref VARCHAR(128),
        loc_start INT UNSIGNED,
        loc_end INT UNSIGNED,
        strand ENUM('+','-'),
        PRIMARY KEY (prot_id),
        INDEX (acce_ver),
        CONSTRAINT fk_proteins_genomes
            FOREIGN KEY (acce_ver) REFERENCES genomes(acce_ver)
            ON DELETE CASCADE ON UPDATE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
    """))

    conn.execute(text("""
    CREATE TABLE IF NOT EXISTS protein_ec (
        prot_id VARCHAR(80) NOT NULL,
        ec VARCHAR(32) NOT NULL,
        PRIMARY KEY (prot_id, ec),
        CONSTRAINT fk_ec_protein
            FOREIGN KEY (prot_id) REFERENCES proteins(prot_id)
            ON DELETE CASCADE ON UPDATE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
    """))

    conn.execute(text("""
    CREATE TABLE IF NOT EXISTS protein_go (
        prot_id VARCHAR(80) NOT NULL,
        go_id VARCHAR(64) NOT NULL,
        PRIMARY KEY (prot_id, go_id),
        CONSTRAINT fk_go_protein
            FOREIGN KEY (prot_id) REFERENCES proteins(prot_id)
            ON DELETE CASCADE ON UPDATE CASCADE
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
    """))

def sanity_checks(conn):
    for t in ["genomes", "proteins", "protein_ec", "protein_go"]:
        n = conn.execute(text(f"SELECT COUNT(*) FROM {t}")).scalar()
        print(f"{t}: {n}")

    # join check: should be <= proteins count
    n = conn.execute(text("""
        SELECT COUNT(*) FROM proteins p
        JOIN genomes g ON g.acce_ver = p.acce_ver
    """)).scalar()
    print("proteins JOIN genomes:", n)

def main():
    engine = create_engine("mysql+mysqldb:///?read_default_file=/homes/aakhlaghi/.my.cnf")

    t0 = time.time()
    with engine.begin() as conn:
        create_tables(conn)
    print("tables ready")

    with engine.begin() as conn:
        for record in SeqIO.parse(FILE_PATH, "genbank"):
            genome = extract_genome(record)
            proteins, ecs, gos = extract_proteins(record, genome["acce_ver"])

            # 1) insert genome (parent)
            conn.execute(text("""
                INSERT INTO genomes (acce_ver, organism, genome_size, gene_count, protein_count, taxon_id)
                VALUES (:acce_ver, :organism, :genome_size, :gene_count, :protein_count, :taxon_id)
                ON DUPLICATE KEY UPDATE
                    organism=VALUES(organism),
                    genome_size=VALUES(genome_size),
                    gene_count=VALUES(gene_count),
                    protein_count=VALUES(protein_count),
                    taxon_id=VALUES(taxon_id)
            """), genome)

            # 2) insert proteins (child)
            if proteins:
                conn.execute(text("""
                    INSERT INTO proteins (prot_id, acce_ver, product, locus_tag, gene_ref, loc_start, loc_end, strand)
                    VALUES (:prot_id, :acce_ver, :product, :locus_tag, :gene_ref, :loc_start, :loc_end, :strand)
                    ON DUPLICATE KEY UPDATE
                        acce_ver=VALUES(acce_ver),
                        product=VALUES(product),
                        locus_tag=VALUES(locus_tag),
                        gene_ref=VALUES(gene_ref),
                        loc_start=VALUES(loc_start),
                        loc_end=VALUES(loc_end),
                        strand=VALUES(strand)
                """), proteins)

            # 3) insert EC/GO (bridge tables)
            if ecs:
                conn.execute(text("""
                    INSERT IGNORE INTO protein_ec (prot_id, ec) VALUES (:prot_id, :ec)
                """), ecs)

            if gos:
                conn.execute(text("""
                    INSERT IGNORE INTO protein_go (prot_id, go_id) VALUES (:prot_id, :go_id)
                """), gos)

    with engine.begin() as conn:
        sanity_checks(conn)

    print(f"done in {time.time() - t0:.1f}s")

if __name__ == "__main__":
    main()
