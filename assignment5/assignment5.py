#!/usr/bin/python3
"""
Assignment 5 - PySpark GenBank Parser
Reads a GBFF file and answers questions using Spark DataFrames.
"""

from Bio import SeqIO
import argparse
import sys

# Spark setup - using cluster mode
sys.path.append('/opt/spark/python')
sys.path.append('/opt/spark/python/lib/py4j-0.10.9.7-src.zip')
from pyspark.sql import SparkSession
from pyspark.sql.types import (
    StructType, StructField, StringType, IntegerType, LongType
)
from pyspark.sql import Row

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
    Extract all features from a GenBank record into a flat structure.
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
    
    # Extract each feature
    for feat in record.features:
        feature_type = feat.type
        
        # Get qualifiers
        q = feat.qualifiers
        prot_id = q.get("protein_id", [None])[0]
        product = q.get("product", [None])[0]
        locus_tag = q.get("locus_tag", [None])[0]
        gene_name = q.get("gene", [None])[0]
        
        # Get location
        loc_start = int(feat.location.start)
        loc_end = int(feat.location.end)
        strand = "+" if feat.location.strand == 1 else "-"
        
        # Get EC and GO
        ec_numbers = ",".join(q.get("EC_number", []))
        go_terms = ",".join(q.get("GO_function", []) + q.get("GO", []))
        
        # Create row
        row = Row(
            acce_ver=acce_ver,
            organism=organism,
            genome_size=genome_size,
            taxon_id=taxon_id,
            feature_type=feature_type,
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


def Print_Question_1(df):
    """
    Q1: How many features does an Archaeal genome have on average?
    """
    # Count features per genome (acce_ver)
    features_per_genome = df.groupBy("acce_ver").count()
    
    # Calculate average
    avg_features = features_per_genome.agg({"count": "avg"}).collect()[0][0]
    
    print(f"\nQ1) Average features per Archaeal genome: {int(avg_features)}")
    return int(avg_features)

def Print_Question_2(df):
    """Q2: What is the ratio between coding and non-coding features?"""
    # Count coding features (CDS)
    coding_count = df.filter(df.feature_type == "CDS").count()
    
    # Count non-coding features (everything except CDS)
    non_coding_count = df.filter(df.feature_type != "CDS").count()
    
    # Calculate ratio
    if non_coding_count == 0:
        ratio = float('inf')
    else:
        ratio = coding_count / non_coding_count
    
    print(f"\nQ2) Ratio of coding to non-coding features:")
    print(f"Ratio (coding/non-coding): {ratio:.2f}")
    
    return ratio

def Print_Question_3(df):
    """Q3: Minimal and maximal number of proteins in a genome for all organisms?"""
    
    # count proteins (CDS) per genome
    proteins_per_genome = (df.filter(df.feature_type == "CDS")
                           .groupBy("acce_ver")
                           .count()
                           .withColumnRenamed("count", "protein_count"))
    
    # get organism name for each genome
    genome_info = df.select("acce_ver", "organism").distinct()
    proteins_with_organism = proteins_per_genome.join(genome_info, "acce_ver")
    
    # find the genome with MINIMUM proteins
    # sort beshe kam be ziad, avalin ro begirim
    min_row = proteins_with_organism.orderBy("protein_count", ascending=True).first()
    
    # find the genome with MAXIMUM proteins
    max_row = proteins_with_organism.orderBy("protein_count", ascending=False).first()
    
    print(f"\nQ3) Minimal and maximal number of proteins per genome:")
    print(f"Minimum: {min_row['protein_count']} proteins")
    print(f"Genome: {min_row['acce_ver']}")
    print(f"Organism: {min_row['organism']}")
    print(f"Maximum: {max_row['protein_count']} proteins")
    print(f"Genome: {max_row['acce_ver']}")
    print(f"Organism: {max_row['organism']}")
    
    return min_row['protein_count'], max_row['protein_count']

def Print_Question_4(df):
    """
    Q4: What is the average length of a feature?
    """
    # Compute feature length as (loc_end - loc_start)
    df_len = df.withColumn("feat_len", df.loc_end - df.loc_start)
    
    # Average length over all features
    avg_len = df_len.agg({"feat_len": "avg"}).collect()[0][0]
    
    print(f"\nQ4) Average feature length (bp): {avg_len:.2f}")
    return float(avg_len)

def Print_Question_5(df):
    df_non_coding = df.filter(df.feature_type != "CDS")
    return df_non_coding

def Print_Question_6(df):
    """
    Q4: What is the average length of a non_coding feature?
    """
    df_non_coding = Print_Question_5(df)

    # Compute feature length as (loc_end - loc_start)
    df_non_coding_len = df_non_coding.withColumn("feat_len", df_non_coding.loc_end - df_non_coding.loc_start)
    
    # Average length over all features
    avg_len_non_coding = df_non_coding_len.agg({"feat_len": "avg"}).collect()[0][0]
    
    print(f"\nQ4) Average feature length (bp): {avg_len_non_coding:.2f}")
    return float(avg_len_non_coding)

def main():
    """Main execution function."""
    print("="*60, flush=True)
    print("Assignment 5 - Spark GenBank Parser", flush=True)
    print("="*60, flush=True)
    
    args = parse_args()
    file_path = args.file
    print(f"Input file: {file_path}", flush=True)
    
    print("\nParsing GenBank file...", flush=True)
    
    # Parse file and extract features
    all_rows = []
    record_count = 0
    
    with open(file_path) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            record_count += 1
            rows = extract_features_from_record(record)
            all_rows.extend(rows)
            
            if record_count % 100 == 0:
                print(f"Processed {record_count} records...", flush=True)
    
    print(f"\nTotal records parsed: {record_count}", flush=True)
    print(f"Total features extracted: {len(all_rows)}", flush=True)
    
    # Create Spark DataFrame
    print("\nCreating Spark DataFrame...", flush=True)
    
    # Define schema
    schema = StructType([
        StructField("acce_ver", StringType(), True),
        StructField("organism", StringType(), True),
        StructField("genome_size", LongType(), True),
        StructField("taxon_id", IntegerType(), True),
        StructField("feature_type", StringType(), True),
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
    df = spark.createDataFrame(all_rows, schema=schema)
    
    print(f"DataFrame created with {df.count()} rows", flush=True)
    
    # Show schema
    print("\nDataFrame Schema:", flush=True)
    df.printSchema()
    
    # Show sample data
    print("\nSample data (first 5 rows):", flush=True)
    df.show(5, truncate=False)
    
    # Answer questions
    print("\n" + "="*60, flush=True)
    print("ANSWERS", flush=True)
    print("="*60, flush=True)
    
    Print_Question_1(df)
    Print_Question_2(df)
    Print_Question_3(df)
    Print_Question_4(df)
    Print_Question_6(df)

    print("\n" + "="*60, flush=True)
    print("Done!", flush=True)
    print("="*60, flush=True)
    
    # Stop Spark
    spark.stop()


if __name__ == "__main__":
    main()