#!/usr/bin/env python3
"""Assignment 6: dbNSFP classifier counts and normalized MariaDB export.

Edit the settings in main() for another CSV/TSV dataset or machine.
A codebook containing metadata is consulted by the student, not joined as rows.
"""

import configparser
from functools import reduce

from pyspark.sql import SparkSession
from pyspark.sql import functions as F


# Reusable exam functions ---------------------------------------------------
def create_spark(app_name="PySpark", master=None, jdbc_jar=None):
    builder = SparkSession.builder.appName(app_name)
    if master:
        builder = builder.master(master)
    if jdbc_jar:
        builder = builder.config("spark.jars", jdbc_jar)
    return builder.getOrCreate()


def read_data(spark, data_path, sep=",", header=True, infer_schema=False):
    return (
        spark.read.option("header", header)
        .option("sep", sep)
        .option("inferSchema", False)
        .csv(data_path)
    )


def recasting(df, col_types, col_rename_map=None, verbose=False):
    """Rename known columns, then cast requested columns; fail on collisions."""
    applied = {}
    for old, new in (col_rename_map or {}).items():
        if old not in df.columns:
            if verbose:
                print(f"missing column for rename: {old}")
            continue
        if old != new:
            if new in df.columns:
                raise ValueError(f"Cannot rename {old} to existing column {new}")
            df = df.withColumnRenamed(old, new)
        applied[old] = new
    for name, dtype in col_types.items():
        name = applied.get(name, name)
        if name not in df.columns:
            if verbose:
                print(f"missing column for cast: {name}")
            continue
        df = df.withColumn(name, F.col(name).cast(dtype))
    return df


def check_dataframe(df, n=5):
    df.printSchema()
    df.show(n, truncate=False)


def join_data(df1, df2, key, how="inner"):
    """Use only for a second row-level dataset with a real shared key."""
    return df1.join(df2, on=key, how=how)


# Assignment 6 queries ----------------------------------------------------
CLASSIFIER_SUFFIXES = ("_converted_rankscore", "_rankscore", "_score", "_pred")
MISSING_TEXT = ("", ".", "NA", "N/A", "null", "nan")


def classifier_groups(df):
    """Group exact classifier suffixes; do not treat rank scores as predictions."""
    groups = {}
    for name in df.columns:
        suffix = next((s for s in CLASSIFIER_SUFFIXES if name.endswith(s)), None)
        if suffix:
            groups.setdefault(name[:-len(suffix)], []).append(name)
    if not groups:
        raise ValueError("No classifier columns found: inspect column names and suffixes")
    return groups


def is_prediction(name):
    """A valid output is nonempty and not a documented missing marker."""
    value = F.trim(F.col(name).cast("string"))
    return F.col(name).isNotNull() & ~F.lower(value).isin(*[s.lower() for s in MISSING_TEXT])


def prediction_flags(groups):
    """One flag per classifier per row, even if it has score and pred columns."""
    flags = {}
    for classifier, columns in groups.items():
        # Rank scores only order predictions; prefer pred, otherwise score.
        candidates = [c for c in columns if c.endswith(("_pred", "_score"))]
        if candidates:
            flags[classifier] = F.when(
                reduce(lambda a, b: a | b, [is_prediction(c) for c in candidates]),
                1,
            ).otherwise(0)
    if not flags:
        raise ValueError("No _pred or _score columns found for classifier counts")
    return flags


def q1_classifier_counts(df, flags):
    counts = df.agg(*[F.sum(flag).alias(name) for name, flag in flags.items()]).first()
    result = {name: int(counts[name] or 0) for name in flags}
    for name, count in sorted(result.items()):
        print(f"{name}: {count}")
    return result


def q2_top5_classifiers(df, groups, counts, keep_columns, n=5):
    top = sorted(counts, key=lambda name: (-counts[name], name))[:n]
    print("Top classifiers:", [(name, counts[name]) for name in top])
    missing = [c for c in keep_columns if c not in df.columns]
    if missing:
        raise ValueError(f"Required columns missing: {missing}")
    selected = keep_columns + [c for name in top for c in groups[name]]
    return df.select(*dict.fromkeys(selected)), top


def q3_position_id(df, chr_col, pos_col, result_col="position_id"):
    """Use string chromosomes: X, Y and MT are valid chromosome labels."""
    if result_col in df.columns:
        raise ValueError(f"Column already exists: {result_col}")
    return df.withColumn(
        result_col,
        F.when(
            F.col(chr_col).isNotNull() & F.col(pos_col).isNotNull(),
            F.concat_ws(":", F.col(chr_col), F.col(pos_col)),
        ),
    )


def top_prediction_group(df, group_col, flags, n=1):
    """Count available classifier outputs, not merely rows in a group."""
    total = reduce(lambda a, b: a + b, flags.values())
    return (
        df.where(F.col(group_col).isNotNull())
        .groupBy(group_col)
        .agg(F.sum(total).alias("prediction_count"))
        .orderBy(F.desc("prediction_count"), F.asc(group_col))
        .limit(n)
    )


def q4_position(df, flags):
    result = top_prediction_group(df, "position_id", flags)
    result.show(truncate=False)
    return result


def q5_proteins(df, flags):
    result = top_prediction_group(df, "Ensembl_proteinid", flags)
    result.show(truncate=False)
    return result


def normalize_3nf(df, selected_groups):
    """Split positions, variants, genes, proteins and classifier results.

    Assumes one atomic gene/protein identifier per cell. The variant key uses
    chromosome, coordinate and alleles; genome assembly must be consistent.
    """
    classifier_cols = [c for cols in selected_groups.values() for c in cols]
    df = df.withColumn(
        "variant_id",
        F.sha2(F.to_json(F.struct(
            F.col("chr"), F.col("pos(1-based)"), F.col("ref"), F.col("alt")
        )), 256),
    )
    return {
        "positions": df.select("position_id", "chr", "pos(1-based)").distinct(),
        "variants": df.select("variant_id", "position_id", "ref", "alt").distinct(),
        "genes": df.select("Ensembl_geneid", "genename").distinct(),
        "proteins": df.select("Ensembl_proteinid", "Ensembl_geneid").distinct(),
        "predictions": df.select("variant_id", "Ensembl_proteinid", *classifier_cols).distinct(),
    }


def read_mysql_config(path):
    """Read the [client] section of a MariaDB/MySQL defaults file."""
    config = configparser.ConfigParser(interpolation=None)
    if not config.read(path) or not config.has_section("client"):
        raise ValueError(f"Missing [client] section in {path}")
    section = config["client"]
    return section["database"], section["user"], section["password"]


def save_to_db(tables, jdbc_url, user, password, mode="overwrite"):
    for name, table in tables.items():
        (table.write.format("jdbc")
         .option("driver", "org.mariadb.jdbc.Driver")
         .option("url", jdbc_url)
         .option("dbtable", name)
         .option("user", user)
         .option("password", password)
         .mode(mode).save())
        print(f"Saved {name}")


def main():
    data_path = "/data/datasets/dbNSFP/snpEff/data/dbNSFP4.9a.txt.gz.SMALL"
    data_sep = "\t"
    data_header = True
    master = "spark://spark.bin.bioinf.nl:7077"
    jdbc_jar = "/usr/share/java/mariadb-java-client.jar"
    mysql_config = "/homes/aakhlaghi/.my.cnf"
    jdbc_host = "mariadb.bin.bioinf.nl"
    save_results = True 

    # For another dataset, edit its column names/types here.
    col_rename_map = {"#chr": "chr"}
    col_types = {"pos(1-based)": "long"}  # Keep chromosome labels as strings.
    keep_columns = [
        "chr", "pos(1-based)", "ref", "alt", "genename",
        "Ensembl_proteinid", "Ensembl_geneid",
    ]

    spark = create_spark("assignment6_Ali", master, jdbc_jar)
    try:
        df = read_data(spark, data_path, sep=data_sep, header=data_header)
        df = recasting(df, col_types, col_rename_map, verbose=True)
        check_dataframe(df)

        groups = classifier_groups(df)
        flags = prediction_flags(groups)
        counts = q1_classifier_counts(df, flags)
        df, top = q2_top5_classifiers(df, groups, counts, keep_columns)
        selected_groups = {name: groups[name] for name in top}
        selected_flags = prediction_flags(selected_groups)

        df = q3_position_id(df, "chr", "pos(1-based)")
        
        # if I want to save answer of a query to file:
        result4 = q4_position(df, selected_flags)

        with open("assignment6_results.txt", "a", encoding="utf-8") as file:
            file.write(f"Q4 - position with most predictions: {result4.first()}\n")

        q5_proteins(df, selected_flags)

        if save_results:
            db_name, user, password = read_mysql_config(mysql_config)
            tables = normalize_3nf(df, selected_groups)
            save_to_db(tables, f"jdbc:mariadb://{jdbc_host}/{db_name}", user, password)
    finally:
        spark.stop()


if __name__ == "__main__":
    main()
