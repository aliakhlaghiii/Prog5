#!/usr/bin/python3
import sys
sys.path.append('/opt/spark/python')
sys.path.append('/opt/spark/python/lib/py4j-0.10.9.7-src.zip')

from sqlalchemy import create_engine, text
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, concat_ws, count, isnull

# MySQL Connection Test
print("Test MySQL connection...", flush=True)
engine = create_engine(
    "mysql+mysqldb:///?read_default_file=/homes/aakhlaghi/.my.cnf"
)
with engine.connect() as conn:
    result = conn.execute(text("SELECT DATABASE();"))
    db_name = result.scalar()
    print(f"Connected successfully to: {db_name}", flush=True)

# Spark Session with JDBC
print("\nInitialize Spark with JDBC driver", flush=True)
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
print("\nRead the data...", flush=True)
file_path = '/data/datasets/dbNSFP/snpEff/data/dbNSFP4.9a.txt.gz'
spark_df = spark.read.csv(
    file_path,
    sep='\t',
    header=True,
    inferSchema=False
)

# Test with first 100 rows
spark_df = spark_df.limit(100)
print(f"Data loaded: {spark_df.count()} rows, {len(spark_df.columns)} columns", flush=True)

# Schema Design (only for data types transformation)
print("\nApply schema transformations...", flush=True)

from pyspark.sql.types import IntegerType, FloatType

# integers
spark_df = spark_df \
    .withColumn("chr", col("#chr").cast(IntegerType())) \
    .withColumn("pos(1-based)", col("pos(1-based)").cast(IntegerType())) \
    .withColumn("hg19_chr", col("hg19_chr").cast(IntegerType())) \
    .withColumn("hg19_pos(1-based)", col("hg19_pos(1-based)").cast(IntegerType())) \
    .withColumn("hg18_chr", col("hg18_chr").cast(IntegerType())) \
    .withColumn("hg18_pos(1-based)", col("hg18_pos(1-based)").cast(IntegerType())) \
    .withColumn("aapos", col("aapos").cast(IntegerType())) \
    .withColumn("codonpos", col("codonpos").cast(IntegerType()))

# floats
score_columns = [c for c in spark_df.columns if '_score' in c.lower() or '_rankscore' in c.lower()]
for score_col in score_columns:
    spark_df = spark_df.withColumn(score_col, col(score_col).cast(FloatType()))

print("Schema applied successfully!", flush=True)

# Q1: How many predictions each of the 43 classifiers makes?
# find all the columns related to classifier
classifier_cols = [c for c in spark_df.columns if '_score' in c or '_pred' in c]

# find non-null values for each classifier
predictions_count = {}
for col_name in classifier_cols:
    count = spark_df.filter(col(col_name).isNotNull()).count()
    predictions_count[col_name] = count
    print(f"{col_name}: {count} predictions", flush=True)
    
# Q2: Make a top five of classifiers based on the total prediction it makes, 
# and drop all the columns associated with the others.

# extract all classifiers and their prediction
classifier_groups = {}
for col_name in classifier_cols:
    # remove the extra suffix or prefix in column name
    base_name = col_name.replace('_score', '').replace('_pred', '').replace('_converted_rankscore', '').replace('_rankscore', '')
    if base_name not in classifier_groups:
        classifier_groups[base_name] = []
    classifier_groups[base_name].append(col_name)

# counting all the predictions for each classifier
classifier_totals = {}
for base_name, cols in classifier_groups.items():
    total = 0
    for col_name in cols:
        if col_name in predictions_count:
            total += predictions_count[col_name]
    classifier_totals[base_name] = total

# sorting and select top5
sorted_classifiers = sorted(classifier_totals.items(), key=lambda x: x[1], reverse=True)
top5 = sorted_classifiers[:5]
print("Top 5 classifiers:", flush=True)
for name, count in top5:
    print(f"  {name}: {count}", flush=True)

# save top5
top5_names = [name for name, _ in top5]
columns_to_keep = ['chr', 'pos(1-based)', 'ref', 'alt', 'genename', 'Ensembl_proteinid', 'Ensembl_geneid']
for base_name in top5_names:
    columns_to_keep.extend(classifier_groups[base_name])

# remove the others
spark_df = spark_df.select([c for c in spark_df.columns if c in columns_to_keep])
print(f"Kept {len(spark_df.columns)} columns", flush=True)

# Q3: Merge the chr and pos columns to yield a unique identifier for each position in the genome.
spark_df = spark_df.withColumn(
    "position_id",
    concat_ws(":", col("chr"), col("pos(1-based)")))

# Q4: What position has the most predictions associated with it?
position_counts = spark_df.groupBy("position_id").count().orderBy("count", ascending=False)
top_position = position_counts.first()
print(f"Position with most predictions: {top_position['position_id']} ({top_position['count']} predictions)", flush=True)
