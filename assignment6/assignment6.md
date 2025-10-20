## Questions

How many predictions does each of the 43 classifiers make?
Make a top five of classifiers based on the total predictions they make, and drop all columns associated with the others.
Merge the chr and pos columns to yield a unique identifier for each position in the genome.
What position has the most predictions associated with it?
What protein (Ensembl_proteinid) has the most predictions associated with it?
Answers (based on the run: 100 rows × 458 columns)
# Note: These counts were obtained from a 100-row sample of the full dbNSFP file.

1) Predictions per classifier (summary):
The *_pred columns for most tools had 100 non-missing values:
SIFT_pred, SIFT4G_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, etc.
The following *_score columns had nonzero counts:
MetaSVM_score = 60, MetaLR_score = 60, MetaRNN_score = 83, LIST-S2_score = 83, AlphaMissense_score = 74, DANN_score = 100, fathmm-MKL_coding_score = 100, GenoCanyon_score = 100, integrated_fitCons_score = 100, GM12878_fitCons_score = 100, H1-hESC_fitCons_score = 100, and HUVEC_fitCons_score = 100.
All other *_score columns had 0 predictions.

2) Top 5 classifiers by total predictions:
fathmm-MKL_coding — 200
MetaRNN — 183
LIST-S2 — 183
AlphaMissense — 174
MetaSVM — 160
(After filtering, 17 columns were retained in the working dataframe.)

3) Unique genomic identifier format:
Each genomic position was represented as chr:pos
—for example, 1:65565

4) Position with the most predictions:
1:65565 — 3 predictions (within the retained top-5 classifiers).

5) Protein with the most predictions:
ENSP00000493376