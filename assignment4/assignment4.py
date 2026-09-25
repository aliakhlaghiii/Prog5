#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from sqlalchemy import create_engine, text


FILE_PATH = "/data/datasets/NCBI/refseq/ftp.ncbi.nlm.nih.gov/refseq/release/archaea/archaea.1.genomic.gbff"

db_file = "/homes/aakhlaghi/exam.db"


def load_data(FILE_PATH):

    rows = []

    for record in SeqIO.parse(FILE_PATH, "genbank"):

        acce = record.annotations.get("accessions", [str(record.id)])[0]

        version = record.annotations.get("sequence_version")

        acce_ver = f"{acce}.{version}" if version else acce

        organism = record.annotations.get("organism", "unknown")

        genome_size = len(record.seq)

        gene_count = sum(
            1 for f in record.features
            if f.type == "gene"
        )

        protein_count = sum(
            1 for f in record.features
            if f.type == "CDS"
        )

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


        for feat in record.features:

            if feat.type != "CDS":
                continue

            q = feat.qualifiers

            prot_id = q.get("protein_id", [None])[0]

            product = q.get("product", [None])[0]

            locus_tag = q.get("locus_tag", [None])[0]

            gene_ref = q.get("gene", [None])[0]

            loc_start = int(feat.location.start)

            loc_end = int(feat.location.end)

            strand = "+" if feat.location.strand == 1 else "-"


            if not prot_id:

                prot_id = f"{acce_ver}:{locus_tag}:{loc_start}-{loc_end}"


            ecs = q.get("EC_number", [])

            ec = ";".join(ecs) if ecs else None


            go_terms = []

            go_terms += q.get("GO", [])

            go_terms += q.get("GO_function", [])

            go_id = ";".join(go_terms) if go_terms else None


            rows.append({

                "acce_ver": acce_ver,

                "organism": organism,

                "genome_size": genome_size,

                "gene_count": gene_count,

                "protein_count": protein_count,

                "taxon_id": taxon_id,

                "prot_id": prot_id,

                "product": product,

                "locus_tag": locus_tag,

                "gene_ref": gene_ref,

                "loc_start": loc_start,

                "loc_end": loc_end,

                "strand": strand,

                "ec": ec,

                "go_id": go_id
            })


    df = pd.DataFrame(rows)

    return df



def extract_manual(df, cols_a, cols_b, *, key, strict=True):

    for cols, name in ((cols_a, 'A'), (cols_b, 'B')):

        if key not in cols:

            raise ValueError(f"{name} is missing the key {key}")


    if strict:

        missing_a = [c for c in cols_a if c not in df.columns]

        missing_b = [c for c in cols_b if c not in df.columns]


        if missing_a or missing_b:

            raise ValueError(
                f"A missing {missing_a}, B missing {missing_b}"
            )


        extra = (set(cols_a) & set(cols_b)) - {key}


        if extra:

            raise ValueError(
                f"non-key overlap: {sorted(extra)}"
            )


    return df[cols_a].copy(), df[cols_b].copy()



def drop_tables(conn):

    conn.execute(
        text("DROP TABLE IF EXISTS proteins;")
    )

    conn.execute(
        text("DROP TABLE IF EXISTS genomes;")
    )



def records(frame):

    clean = frame.astype(object).where(
        frame.notna(),
        None
    )

    return clean.to_dict(
        orient='records'
    )



def insert_one(conn, table_name, frame):

    cols = list(frame.columns)

    col_list = ",".join(cols)

    val_list = ",".join(
        f":{c}" for c in cols
    )


    conn.execute(

        text(
            f"""
            INSERT INTO {table_name}
            ({col_list})
            VALUES ({val_list});
            """
        ),

        records(frame),
    )



def insert_to_table(
    conn,
    table1_name,
    table1_df,
    table2_name,
    table2_df
):

    insert_one(
        conn,
        table1_name,
        table1_df
    )

    insert_one(
        conn,
        table2_name,
        table2_df
    )



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

        ec TEXT,

        go_id TEXT,

        PRIMARY KEY (prot_id),

        INDEX (acce_ver),

        CONSTRAINT fk_proteins_genomes

            FOREIGN KEY (acce_ver)

            REFERENCES genomes(acce_ver)

            ON DELETE CASCADE

            ON UPDATE CASCADE

    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
    """))



def query(conn):

    rows = conn.execute(text("""

        SELECT

            g.organism,

            g.acce_ver,

            g.genome_size,

            g.gene_count,

            g.protein_count,

            g.taxon_id,

            p.prot_id,

            p.product,

            p.loc_start,

            p.loc_end,

            p.strand,

            p.locus_tag,

            p.gene_ref,

            p.ec,

            p.go_id

        FROM genomes g

        JOIN proteins p

        ON g.acce_ver = p.acce_ver

    """)).fetchall()


    print("result:", rows)



def main():

    engine = create_engine(
        f"sqlite:///{db_file}"
    )


    df = load_data(FILE_PATH)


    genomes_cols = [

        'acce_ver',

        'organism',

        'genome_size',

        'gene_count',

        'protein_count',

        'taxon_id'
    ]


    proteins_cols = [

        'prot_id',

        'acce_ver',

        'product',

        'locus_tag',

        'gene_ref',

        'loc_start',

        'loc_end',

        'strand',

        'ec',

        'go_id'
    ]


    genomes_df, proteins_df = extract_manual(

        df,

        genomes_cols,

        proteins_cols,

        key='acce_ver'
    )


    genomes_df = genomes_df.drop_duplicates()


    with engine.begin() as conn:

        drop_tables(conn)

        create_tables(conn)

        insert_to_table(

            conn,

            'genomes',

            genomes_df,

            'proteins',

            proteins_df
        )


        query(conn)



if __name__ == "__main__":

    main()