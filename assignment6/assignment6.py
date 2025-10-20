#!/usr/bin/python3
"""Process dbNSFP data using Spark and save to MariaDB."""

import sys
import subprocess
from pyspark.sql import SparkSession
from pyspark.sql.functions import (col, concat_ws,
                                    monotonically_increasing_id)
from pyspark.sql.types import IntegerType, FloatType
from sqlalchemy import create_engine, text

sys.path.append('/opt/spark/python')
sys.path.append('/opt/spark/python/lib/py4j-0.10.9.7-src.zip')


def get_mysql_config(key):
    """Extract MySQL configuration from .my.cnf file."""
    result = subprocess.run(
        ['grep', f'^{key}', '/homes/aakhlaghi/.my.cnf'],
        capture_output=True,
        text=True,
        check=False
    )
    return (result.stdout.split('=')[1].strip()
            if result.stdout else None)


def test_mysql_connection():
    """Test connection to MySQL database."""
    print("Test MySQL connection...", flush=True)
    engine = create_engine(
        "mysql+mysqldb:///?read_default_file="
        "/homes/aakhlaghi/.my.cnf"
    )
    with engine.connect() as conn:
        result = conn.execute(text("SELECT DATABASE();"))
        db_name = result.scalar()
        print(f"Connected successfully to: {db_name}",
              flush=True)
        return db_name


def initialize_spark():
    """Initialize Spark session with JDBC driver."""
    print("\nInitialize Spark with JDBC driver", flush=True)
    spark = (SparkSession.builder
        .appName("assignment6_Ali")
        .master("spark://spark.bin.bioinf.nl:7077")
        .config("spark.jars",
                "/usr/share/java/mariadb-java-client.jar")
        .getOrCreate())

    print(f"Spark initialized: {spark.sparkContext.master}",
          flush=True)
    return spark


def load_and_prepare_data(spark):
    """Load data and apply schema transformations."""
    print("\nRead the data...", flush=True)
    file_path = (
        '/data/datasets/dbNSFP/snpEff/data/dbNSFP4.9a.txt.gz'
    )

    spark_df = spark.read.csv(
        file_path,
        sep='\t',
        header=True,
        inferSchema=False
    )

    spark_df = spark_df.limit(100)
    print(
        f"Data loaded: {spark_df.count()} rows, "
        f"{len(spark_df.columns)} columns",
        flush=True
    )

    print("\nApply schema transformations...", flush=True)

    int_cols = [
        "chr", "pos(1-based)", "hg19_chr",
        "hg19_pos(1-based)",
        "hg18_chr", "hg18_pos(1-based)", "aapos",
        "codonpos"
    ]
    for int_col in int_cols:
        spark_df = spark_df.withColumn(
            int_col,
            col(int_col).cast(IntegerType())
        )

    score_columns = [
        c for c in spark_df.columns
        if '_score' in c.lower() or '_rankscore' in c.lower()
    ]
    for score_col in score_columns:
        spark_df = spark_df.withColumn(
            score_col,
            col(score_col).cast(FloatType())
        )

    print("Schema applied successfully!", flush=True)
    return spark_df


def analyze_classifiers(spark_df):
    """Analyze classifier predictions and keep top 5."""
    print("\nQ1: Classifier prediction counts", flush=True)

    classifier_cols = [
        c for c in spark_df.columns
        if '_score' in c or '_pred' in c
    ]

    predictions_count = {}
    for classifier_col in classifier_cols:
        count = spark_df.filter(
            col(classifier_col).isNotNull()
        ).count()
        predictions_count[classifier_col] = count
        print(f"{classifier_col}: {count} predictions",
              flush=True)

    print("\nQ2: Top 5 classifiers", flush=True)

    classifier_groups = {}
    for classifier_col in classifier_cols:
        base_name = (classifier_col
                     .replace('_score', '')
                     .replace('_pred', '')
                     .replace('_converted_rankscore', '')
                     .replace('_rankscore', ''))
        if base_name not in classifier_groups:
            classifier_groups[base_name] = []
        classifier_groups[base_name].append(classifier_col)

    classifier_totals = {}
    for base_name, cols in classifier_groups.items():
        total = sum(
            predictions_count.get(classifier_col, 0)
            for classifier_col in cols
        )
        classifier_totals[base_name] = total

    sorted_classifiers = sorted(
        classifier_totals.items(),
        key=lambda x: x[1],
        reverse=True
    )
    top5 = sorted_classifiers[:5]

    for name, count in top5:
        print(f"  {name}: {count}", flush=True)

    top5_names = [name for name, _ in top5]
    columns_to_keep = [
        'chr', 'pos(1-based)', 'ref', 'alt', 'genename',
        'Ensembl_proteinid', 'Ensembl_geneid'
    ]
    for base_name in top5_names:
        columns_to_keep.extend(classifier_groups[base_name])

    kept_cols = [
        col_name for col_name in spark_df.columns
        if col_name in columns_to_keep
    ]
    spark_df = spark_df.select(kept_cols)
    print(f"Kept {len(spark_df.columns)} columns", flush=True)

    return spark_df, classifier_groups


def analyze_positions(spark_df):
    """Q3-Q4: Create position IDs and find most
    predicted positions."""
    print("\nQ3: Creating position identifiers", flush=True)

    spark_df = spark_df.withColumn(
        "position_id",
        concat_ws(":", col("chr"), col("pos(1-based)"))
    )

    print("\nQ4: Position with most predictions", flush=True)

    position_counts = (spark_df.groupBy("position_id").count()
                       .orderBy("count", ascending=False))
    top_position = position_counts.first()
    print(
        f"Position: {top_position['position_id']} "
        f"({top_position['count']} predictions)",
        flush=True
    )

    return spark_df


def analyze_proteins(spark_df):
    """Q5: Find protein with most predictions."""
    print("\nQ5: Protein with most predictions", flush=True)

    protein_counts = (spark_df.groupBy("Ensembl_proteinid")
                      .count()
                      .orderBy("count", ascending=False))
    top_protein = protein_counts.first()
    print(
        f"Protein: {top_protein['Ensembl_proteinid']} "
        f"({top_protein['count']} predictions)",
        flush=True
    )


def normalize_and_save(spark_df, classifier_groups, db_name):
    """Normalize to 3NF and save to database."""
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
    print(f"Genes table: {genes_df.count()} unique genes",
          flush=True)

    spark_df = spark_df.join(
        genes_df.select("Ensembl_geneid", "gene_id"),
        on="Ensembl_geneid",
        how="left"
    )

    score_pred_cols = [
        c for c in spark_df.columns
        if '_score' in c or '_pred' in c
    ]
    variant_cols = (
        ["chr", "pos(1-based)", "ref", "alt",
         "position_id", "gene_id"] +
        score_pred_cols
    )
    variants_df = spark_df.select(*variant_cols)
    print(f"Variants table: {variants_df.count()} variants",
          flush=True)

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
    """Main execution function."""
    db_name = test_mysql_connection()
    spark = initialize_spark()

    spark_df = load_and_prepare_data(spark)
    spark_df, classifier_groups = analyze_classifiers(spark_df)
    spark_df = analyze_positions(spark_df)
    analyze_proteins(spark_df)

    normalize_and_save(spark_df, classifier_groups, db_name)


if __name__ == "__main__":
    main()