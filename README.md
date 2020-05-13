# AMASC

AMASC (automated marker analysis for single-cell RNA-seq) is a workflow for identifying a gene sets that can be applied in cell-type annotation for single-cell RNA-seq data.  

# Prerequisites

* R >= 3.6

* Python > 3.7

* XGBoost (Python)

* scikit-learn

* numpy

* pandas

* matplotlib

# Usage

1. Prepare CITE-Seq data sets in the CSV format as `mydataset1_rna.csv` and `mydataset1_pe.csv`.
2. Modify the paths in `AMASC.R`
3. Run `AMASC.R`
4. The selected features will be in file `features_time.txt`

Please note that the process is stochastic. 
