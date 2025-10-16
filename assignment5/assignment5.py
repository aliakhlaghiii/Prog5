#!/usr/bin/python3
"""
Assignment 5 - PySpark GenBank Parser
Reads a GBFF file and answers questions using Spark DataFrames.
Only processes 4 feature types: CDS, ncRNA, rRNA, and gene (without CDS).
"""

import sys
import argparse
from Bio import SeqIO

# Spark setup - using cluster mode
sys.path.append('/opt/spark/python')
sys.path.append('/opt/spark/python/lib/py4j-0.10.9.7-src.zip')
from pyspark.sql import SparkSession, Row
from pyspark.sql.types import (
    StructType, StructField, StringType, IntegerType, LongType
)

print("Initializing Spark...", flush=True)
spark = (SparkSession.builder
         .appName("assignment5_Ali")
         .master("spark://spark.bin.bioinf.nl:7077")
         .getOrCreate())
sc = spark.sparkContext
print(f"Spark initialized: {spark.sparkContext.master}", flush=True)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse GenBank file with Spark")
    parser.add_argument("file", help="Path to .gbff file")
    return parser.parse_args()


def extract_features_from_record(record):
    """
    Extract ONLY relevant features from a GenBank record:
    1. CDS (coding)
    2. ncRNA (non-coding)
    3. rRNA (non-coding)
    4. gene without CDS (non-coding/cryptic genes)

    Returns a list of Row objects (one per feature).
    """
    rows = []

    # Get genome-level info
    acce = record.annotations.get("accessions", [str(record.id)])[0]
    version = record.annotations.get("sequence_version")
    acce_ver = f"{acce}.{version}" if version else acce
    organism = record.annotations.get("organism", "unknown")
    genome_size = len(record.seq)

    # Get taxon_id
    taxon_id = None
    for feat in record.features:
        if feat.type == "source":
            for xref in feat.qualifiers.get("db_xref", []):
                if xref.startswith("taxon:"):
                    taxon_id = int(xref.split(":", maxsplit=1)[1])
                    break
            break

    # Collect all CDS locus_tags to identify cryptic genes
    cds_locus_tags = set()
    for feat in record.features:
        if feat.type == "CDS":
            locus = feat.qualifiers.get("locus_tag", [None])[0]
            if locus:
                cds_locus_tags.add(locus)

    # Extract ONLY relevant features
    for feat in record.features:
        feature_type = feat.type

        # Filter: only process CDS, ncRNA, rRNA, and gene
        if feature_type not in ["CDS", "ncRNA", "rRNA", "gene"]:
            continue

        # For gene: skip if it has a corresponding CDS (not cryptic)
        if feature_type == "gene":
            locus = feat.qualifiers.get("locus_tag", [None])[0]
            if locus and locus in cds_locus_tags:
                continue  # This gene has a CDS, so skip it

        # Get qualifiers
        qualifiers = feat.qualifiers
        prot_id = qualifiers.get("protein_id", [None])[0]
        product = qualifiers.get("product", [None])[0]
        locus_tag = qualifiers.get("locus_tag", [None])[0]
        gene_name = qualifiers.get("gene", [None])[0]

        # Get location
        loc_start = int(feat.location.start)
        loc_end = int(feat.location.end)
        strand = "+" if feat.location.strand == 1 else "-"

        # Get EC and GO
        ec_numbers = ",".join(qualifiers.get("EC_number", []))
        go_terms = ",".join(qualifiers.get("GO_function", []) +
                           qualifiers.get("GO", []))

        # Determine if coding or non-coding
        # CDS = coding, everything else = non-coding
        is_coding = (feature_type == "CDS")

        # Create row
        row = Row(
            acce_ver=acce_ver,
            organism=organism,
            genome_size=genome_size,
            taxon_id=taxon_id,
            feature_type=feature_type,
            is_coding=is_coding,
            prot_id=prot_id,
            product=product,
            locus_tag=locus_tag,
            gene_name=gene_name,
            loc_start=loc_start,
            loc_end=loc_end,
            strand=strand,
            ec_numbers=ec_numbers if ec_numbers else None,
            go_terms=go_terms if go_terms else None
        )
        rows.append(row)

    return rows


def print_question_1(dataframe):
    """
    Q1: How many features does an Archaeal genome have on average?
    """
    # Count features per genome (acce_ver)
    features_per_genome = dataframe.groupBy("acce_ver").count()

    # Calculate average
    avg_features = features_per_genome.agg({"count": "avg"}).collect()[0][0]

    print(f"\nQ1) Average features per Archaeal genome: {int(avg_features)}")
    return int(avg_features)


def print_question_2(dataframe):
    """Q2: What is the ratio between coding and non-coding features?"""
    # Count coding features
    coding_count = dataframe.filter(dataframe.is_coding).count()

    # Count non-coding features
    non_coding_count = dataframe.filter(~dataframe.is_coding).count()

    # Calculate ratio
    if non_coding_count == 0:
        ratio = float('inf')
    else:
        ratio = coding_count / non_coding_count

    print(f"\nQ2) Ratio of coding to non-coding features:")
    print(f"  Coding features: {coding_count}")
    print(f"  Non-coding features: {non_coding_count}")
    print(f"  Ratio (coding/non-coding): {ratio:.2f}")

    return ratio


def print_question_3(dataframe):
    """Q3: Minimal and maximal number of proteins in a genome for all organisms?"""

    # Count proteins (CDS) per genome
    proteins_per_genome = (dataframe.filter(dataframe.feature_type == "CDS")
                           .groupBy("acce_ver")
                           .count()
                           .withColumnRenamed("count", "protein_count"))

    # Get organism name for each genome
    genome_info = dataframe.select("acce_ver", "organism").distinct()
    proteins_with_organism = proteins_per_genome.join(genome_info, "acce_ver")

    # Find the genome with MINIMUM proteins
    min_row = proteins_with_organism.orderBy("protein_count", ascending=True).first()

    # Find the genome with MAXIMUM proteins
    max_row = proteins_with_organism.orderBy("protein_count", ascending=False).first()

    print(f"\nQ3) Minimal and maximal number of proteins per genome:")
    print(f"Minimum: {min_row['protein_count']} proteins")
    print(f"  Genome: {min_row['acce_ver']}")
    print(f"  Organism: {min_row['organism']}")
    print(f"Maximum: {max_row['protein_count']} proteins")
    print(f"  Genome: {max_row['acce_ver']}")
    print(f"  Organism: {max_row['organism']}")

    return min_row['protein_count'], max_row['protein_count']


def print_question_4(dataframe):
    """
    Q4: What is the average length of a feature?
    """
    # Compute feature length as (loc_end - loc_start)
    df_len = dataframe.withColumn("feat_len", dataframe.loc_end - dataframe.loc_start)

    # Average length over all features
    avg_len = df_len.agg({"feat_len": "avg"}).collect()[0][0]

    print(f"\nQ4) Average feature length (bp): {avg_len:.2f}")
    return float(avg_len)


def print_question_5(dataframe):
    """
    Q5: Remove all non-coding (RNA) features and save the cleaned-up version.
    """
    # Filter to keep only coding features
    df_coding_only = dataframe.filter(dataframe.is_coding)

    coding_count = df_coding_only.count()

    print(f"\nQ5) Cleaned-up DataFrame created (coding features only)")
    print(f"  Total coding features: {coding_count}")

    return df_coding_only


def print_question_6(df_cleaned):
    """
    Q6: What is the average length of a feature in the cleaned-up version?
    """
    # Compute feature length as (loc_end - loc_start)
    df_len = df_cleaned.withColumn("feat_len", df_cleaned.loc_end - df_cleaned.loc_start)

    # Average length over all features
    avg_len = df_len.agg({"feat_len": "avg"}).collect()[0][0]

    print(f"\nQ6) Average feature length in cleaned-up DataFrame (bp): {avg_len:.2f}")
    return float(avg_len)


def main():
    """Main execution function."""
    print("="*60, flush=True)
    print("Assignment 5 - Spark GenBank Parser", flush=True)
    print("="*60, flush=True)

    args = parse_args()
    file_path = args.file
    print(f"Input file: {file_path}", flush=True)

    # Check if file exists
    try:
        with open(file_path, encoding='utf-8') as file:
            pass
    except FileNotFoundError:
        print(f"ERROR: File not found: {file_path}", flush=True)
        sys.exit(1)

    print("\nParsing GenBank file...", flush=True)
    print("Processing only: CDS, ncRNA, rRNA, and cryptic genes", flush=True)

    # Parse file and extract features
    all_rows = []
    record_count = 0

    try:
        with open(file_path, encoding='utf-8') as handle:
            for record in SeqIO.parse(handle, "genbank"):
                record_count += 1
                rows = extract_features_from_record(record)
                all_rows.extend(rows)

                if record_count % 100 == 0:
                    print(f"Processed {record_count} records...", flush=True)
    except Exception as error:
        print(f"ERROR while parsing file: {error}", flush=True)
        sys.exit(1)

    print(f"\nTotal records parsed: {record_count}", flush=True)
    print(f"Total features extracted: {len(all_rows)}", flush=True)

    if len(all_rows) == 0:
        print("ERROR: No features extracted from file!", flush=True)
        sys.exit(1)

    # Create Spark DataFrame
    print("\nCreating Spark DataFrame...", flush=True)

    # Define schema
    schema = StructType([
        StructField("acce_ver", StringType(), True),
        StructField("organism", StringType(), True),
        StructField("genome_size", LongType(), True),
        StructField("taxon_id", IntegerType(), True),
        StructField("feature_type", StringType(), True),
        StructField("is_coding", StringType(), True),
        StructField("prot_id", StringType(), True),
        StructField("product", StringType(), True),
        StructField("locus_tag", StringType(), True),
        StructField("gene_name", StringType(), True),
        StructField("loc_start", IntegerType(), True),
        StructField("loc_end", IntegerType(), True),
        StructField("strand", StringType(), True),
        StructField("ec_numbers", StringType(), True),
        StructField("go_terms", StringType(), True),
    ])

    # Create DataFrame
    dataframe = spark.createDataFrame(all_rows, schema=schema)

    print(f"DataFrame created with {dataframe.count()} rows", flush=True)

    # Show schema
    print("\nDataFrame Schema:", flush=True)
    dataframe.printSchema()

    # Show sample data
    print("\nSample data (first 5 rows):", flush=True)
    dataframe.show(5, truncate=False)

    # Answer questions
    print("\n" + "="*60, flush=True)
    print("ANSWERS", flush=True)
    print("="*60, flush=True)

    print_question_1(dataframe)
    print_question_2(dataframe)
    print_question_3(dataframe)
    print_question_4(dataframe)
    df_cleaned = print_question_5(dataframe)
    print_question_6(df_cleaned)

    print("\n" + "="*60, flush=True)
    print("Done!", flush=True)
    print("="*60, flush=True)

    # Stop Spark
    spark.stop()


if __name__ == "__main__":
    main()