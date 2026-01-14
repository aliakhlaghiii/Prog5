#!/usr/bin/python3
"""Process dbNSFP data using Spark and save to MariaDB."""

import sys
import subprocess
sys.path.append('/opt/spark/python')
sys.path.append('/opt/spark/python/lib/py4j-0.10.9.7-src.zip')

from sqlalchemy import create_engine, text
from pyspark.sql import SparkSession
from pyspark.sql.functions import (col, concat_ws,
                                    monotonically_increasing_id)  
from pyspark.sql.types import IntegerType, FloatType

# sql default setting
def get_mysql_config(key):
    """Get MySQL config value from .my.cnf."""
    result = subprocess.run(
        ['grep', f'^{key}', '/homes/aakhlaghi/.my.cnf'],
        capture_output=True,
        text=True,
        check=False
    )
    if result.stdout:
        return result.stdout.split('=')[1].strip()
    return None

# connect to data base
def test_mysql():
    """Test MySQL connection and return database name."""
    print("Test MySQL connection...", flush=True)
    engine = create_engine(
        "mysql+mysqldb:///?read_default_file="
        "/homes/aakhlaghi/.my.cnf"
    )
    with engine.connect() as conn:
        result = conn.execute(text("SELECT DATABASE();"))
        db_name = result.scalar()
        print(f"Connected to: {db_name}", flush=True)
        return db_name

# spark default setting
def init_spark():
    """Initialize and return Spark session."""
    print("\nInitialize Spark with JDBC driver", flush=True)
    spark = (SparkSession.builder
        .appName("assignment6_Ali")
        .master("spark://spark.bin.bioinf.nl:7077")
        .config("spark.jars",
                "/usr/share/java/mariadb-java-client.jar")
        .getOrCreate())
    print(f"Spark: {spark.sparkContext.master}", flush=True)
    return spark


def load_data(spark):
    file_path = "/data/datasets/dbNSFP/snpEff/data/dbNSFP4.9a.txt.gz.SMALL"
    df = spark.read.csv(
        file_path,
        sep="\t",
        header=True,
        inferSchema=False
    )
    return df


def apply_schema(spark_df):
    """Cast columns to proper types to apply logic/mathematical works"""
    print("\nApply schema transformations...", flush=True)
    int_cols = [
        "#chr", "pos(1-based)", "hg19_chr", "hg19_pos(1-based)",
        "hg18_chr", "hg18_pos(1-based)", "aapos", "codonpos"
    ]
    for int_col in int_cols:
        spark_df = spark_df.withColumn(
            int_col,
            col(int_col).cast(IntegerType())
        )
    # cast classifier scores columns
    score_columns = [
        c for c in spark_df.columns
        if '_score' in c.lower() or '_rankscore' in c.lower()
    ]
    for score_col in score_columns:
        spark_df = spark_df.withColumn(
            score_col,
            col(score_col).cast(FloatType())
        )
    print("Schema applied!", flush=True)
    return spark_df

# extract classifiers columns
def get_classifier_cols(spark_df):
    """Get all classifier column names."""
    return [c for c in spark_df.columns if '_score' in c or '_pred' in c]


def q1_classifier_counts(spark_df, classifier_cols):
    """Q1: Count predictions per classifier."""
    print("\nQ1: Classifier counts", flush=True)
    predictions_count = {}   #a dict to keep the classifier with its number of predictions
    for col_name in classifier_cols:
        count = spark_df.filter(   # number of predictions per classifier
            col(col_name).isNotNull()   # we have to use not-null values since in codebook some columns have Nan
        ).count()
        predictions_count[col_name] = count   # allocate number of predictions to each classifier
        print(f"{col_name}: {count}", flush=True)
    return predictions_count


def q2_top5_classifiers(classifier_cols, predictions_count):
    """Q2: Get top 5 classifiers and their columns to keep."""
    print("\nQ2: Top 5 classifiers", flush=True)
    classifier_groups = {}
    for col_name in classifier_cols:
        base_name = (col_name
                     .replace('_score', '')   # used replace so that all the names are consistent
                     .replace('_pred', '')
                     .replace('_converted_rankscore', '')
                     .replace('_rankscore', ''))
        if base_name not in classifier_groups:
            classifier_groups[base_name] = []
        classifier_groups[base_name].append(col_name)

    classifier_totals = {}
    for base_name, cols in classifier_groups.items():
        total = sum(
            predictions_count.get(col_name, 0)
            for col_name in cols
        )
        classifier_totals[base_name] = total

    sorted_classifiers = sorted(
        classifier_totals.items(),
        key=lambda x: x[1],   # as same as key=dictname.get()
        reverse=True
    )
    top5 = sorted_classifiers[:5]
    for name, count in top5:
        print(f"  {name}: {count}", flush=True)

    return top5, classifier_groups


def filter_top5_cols(spark_df, top5, classifier_groups):
    """Keep only top 5 classifier columns."""
    top5_names = [name for name, _ in top5]
    columns_to_keep = [   # These columns must be kept under all circumstances
        'chr', 'pos(1-based)', 'ref', 'alt', 'genename',
        'Ensembl_proteinid', 'Ensembl_geneid'
    ]
    for base_name in top5_names:
        columns_to_keep.extend(classifier_groups[base_name])   # extend in spark is like append in py

    spark_df = spark_df.select(   # rebuild the dataframe
        [c for c in spark_df.columns if c in columns_to_keep]
    )
    print(f"Kept {len(spark_df.columns)} columns", flush=True)
    return spark_df


def q3_q4_positions(spark_df):
    """Q3-Q4: Create position ID and find top position based on its redundancy."""
    print("\nQ3: Position identifiers", flush=True)
    spark_df = spark_df.withColumn(   # withColumn:make a new column named anything(here is position_id)
        "position_id",
        concat_ws(":", col("chr"), col("pos(1-based)"))   # concat 2 cols with ':' as new column values
    )

    print("\nQ4: Top position", flush=True)
    position_counts = (spark_df.groupBy("position_id").count()
                       .orderBy("count", ascending=False))   #orderby in spark is as same as sort in py
    top_position = position_counts.first()   # it returns the first row in spark like row[0]
    print(
        f"Position: {top_position['position_id']} "
        f"({top_position['count']})",
        flush=True
    )
    return spark_df


def q5_proteins(spark_df):
    """Q5: Find protein with most predictions."""
    print("\nQ5: Top protein", flush=True)
    protein_counts = (spark_df.groupBy("Ensembl_proteinid")
                      .count()
                      .orderBy("count", ascending=False))
    top_protein = protein_counts.first()
    print(
        f"Protein: {top_protein['Ensembl_proteinid']} "
        f"({top_protein['count']})",
        flush=True
    )


def normalize_3nf(spark_df):
    """Normalize to 3NF: create genes and variants tables."""
    print("\nNormalizing to 3NF...", flush=True)

    genes_df = spark_df.select(
        "Ensembl_geneid",
        "genename",
        "Ensembl_proteinid"
    ).distinct()

    genes_df = genes_df.withColumn(
        "gene_id",
        monotonically_increasing_id()
    )
    print(f"Genes: {genes_df.count()}", flush=True)

    spark_df = spark_df.join(
        genes_df.select("Ensembl_geneid", "gene_id"),
        on="Ensembl_geneid",
        how="left"
    )

    variants_df = spark_df.select(
        "chr", "pos(1-based)", "ref", "alt", "position_id",
        "gene_id",
        *[c for c in spark_df.columns
          if '_score' in c or '_pred' in c]
    )
    print(f"Variants: {variants_df.count()}", flush=True)

    return genes_df, variants_df


def save_to_db(genes_df, variants_df, db_name):
    """Save tables to MariaDB."""
    print("\nSaving to database...", flush=True)
    db_user = get_mysql_config('user')
    db_password = get_mysql_config('password')
    print(f"Using user: {db_user}", flush=True)

    jdbc_url = (
        f"jdbc:mariadb://mariadb.bin.bioinf.nl/{db_name}"
    )

    genes_df.write.format("jdbc") \
        .option("driver", "org.mariadb.jdbc.Driver") \
        .option("url", jdbc_url) \
        .option("dbtable", "genes") \
        .option("user", db_user) \
        .option("password", db_password) \
        .mode("overwrite") \
        .save()

    variants_df.write.format("jdbc") \
        .option("driver", "org.mariadb.jdbc.Driver") \
        .option("url", jdbc_url) \
        .option("dbtable", "variants") \
        .option("user", db_user) \
        .option("password", db_password) \
        .mode("overwrite") \
        .save()

    print('Saved successfully!', flush=True)


def main():
    """Main execution."""
    db_name = test_mysql()
    spark = init_spark()

    spark_df = load_data(spark)
    spark_df = apply_schema(spark_df)

    classifier_cols = get_classifier_cols(spark_df)
    predictions_count = q1_classifier_counts(spark_df, classifier_cols)

    top5, classifier_groups = q2_top5_classifiers(
        classifier_cols, predictions_count
    )
    spark_df = filter_top5_cols(spark_df, top5, classifier_groups)

    spark_df = q3_q4_positions(spark_df)
    q5_proteins(spark_df)

    genes_df, variants_df = normalize_3nf(spark_df)
    save_to_db(genes_df, variants_df, db_name)


if __name__ == "__main__":
    main()