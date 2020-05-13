# AMASC

AMASC (Automated Marker Analysis for Single-Cell RNA-seq) is a workflow for identifying a gene sets that can be applied in cell-type annotation for single-cell RNA-seq data.  

# Prerequisites

* R >= 3.6
* Python > 3.7
* XGBoost (Python)
* scikit-learn
* numpy
* pandas
* matplotlib

# Usage

1. Prepare preprocessed (library-size-adjusted, log1p-transformed) CITE-Seq data sets in the CSV format as `<MY_DATASET>_rna.csv` and `<MY_DATASET>_pe.csv`
2. Modify the paths in `AMASC.R`
3. Run `AMASC.R`
4. The selected features will be in file `features_<TIME>.txt`

Please note that the process is stochastic and the parameters may require adjustments by data set.

# Contributors

AMASC was developed by Tai-Hsien Ou Yang, Wei-Yi Cheng, and James Cai
