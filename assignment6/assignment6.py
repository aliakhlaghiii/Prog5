#!/usr/bin/python3
import sys
sys.path.append('/opt/spark/python')
sys.path.append('/opt/spark/python/lib/py4j-0.10.9.7-src.zip')

from sqlalchemy import create_engine, text
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, concat_ws, count, isnull

# MySQL Connection Test
print("Step 1: Test MySQL connection...", flush=True)
engine = create_engine(
    "mysql+mysqldb:///?read_default_file=/homes/aakhlaghi/.my.cnf"
)
with engine.connect() as conn:
    result = conn.execute(text("SELECT DATABASE();"))
    db_name = result.scalar()
    print(f"Connected successfully to: {db_name}", flush=True)

# Spark Session with JDBC
print("\nStep 2: Initialize Spark with JDBC driver", flush=True)
spark = (SparkSession.builder
    .appName("assignment6_Ali")
    .master("spark://spark.bin.bioinf.nl:7077")
    .config("spark.jars", "/usr/share/java/mariadb-java-client.jar")
    .getOrCreate())

sc = spark.sparkContext
print(f"Spark initialized: {spark.sparkContext.master}", flush=True)

# review the dataset in bash
# zcat /data/datasets/dbNSFP/snpEff/data/dbNSFP4.9a.txt.gz | head -n 5 | column -t -s $'\t' | less -S

# Read Data
print("\nStep 3: Read the data...", flush=True)
file_path = '/data/datasets/dbNSFP/snpEff/data/dbNSFP4.9a.txt.gz'
spark_df = spark.read.csv(
    file_path,
    sep='\t',
    header=True,
    inferSchema=False
)

# Test with first 100 rows
spark_df = spark_df.limit(100)
print(spark_df)
print(f"Data loaded: {spark_df.count()} rows, {len(spark_df.columns)} columns", flush=True)
