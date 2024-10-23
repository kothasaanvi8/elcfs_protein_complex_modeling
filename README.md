# README

## Overview
This pipeline predicts protein complexes by modeling protein-protein interactions and assigning scores to various
protein pairs using available association and expression data. It constructs weighted graphs and then identifies
clusters characterized by high node connectivity and subgraph cohesion.

### Key files of interest
- <span style="font-weight: normal;">[pair probability scores](https://drive.google.com/file/d/1Cn5zRitrl4i1S_WU2DBKCpR1Aj0_bLmb/view?usp=sharing)</span>
- <span style="font-weight: normal;">[predicted clusters](https://drive.google.com/file/d/1-INlBd1bCZiYBd7vwJqvjYiBB3QpEoNK/view?usp=sharing)</span>
- <span style="font-weight: normal;">[Entrez to UniProt ID mapper](https://drive.google.com/file/d/1tT28krOYbgxQ2S3RtSvuA6cG2SA91AN0/view?usp=sharing)</span>
- <span style="font-weight: normal;">[NCI 60 non- and cell-specific features](https://drive.google.com/file/d/1-A66udgX7YWC-Oh7luKY5SzJuhKjPDKR/view?usp=sharing)</span>

### [GitHub repository](https://github.com/garywilkins/elcfs_protein_complex_modeling): <span style="font-weight: normal; font-size: 1em;">garywilkins/elcfs_protein_complex_modeling</span>

### [Complete repository directory](https://drive.google.com/drive/folders/1OEa0NVYCYnhOZmRa_pycIlpNx9ySJIea?usp=drive_link): <span style="font-weight: normal; font-size: 1em;">contains Python script libraries, Colab/jupyter notebooks, and all generated files</span>

### Directory organization: 
- [utils](https://drive.google.com/drive/folders/1uTxp5JFcyc3ZuDj0Y2dj8H61KWbdavvt?usp=drive_link): Python scripts
- [google_colab_notebooks](https://drive.google.com/drive/folders/1t_LPl1zaOkumbwkWCjw3zZ4TTroMkzWH?usp=drive_link): Jupyter notebooks to run the pipeline
- [CFinder-2.0.6--1448](https://drive.google.com/drive/folders/1v2vkPTisRzkhYpw0knofvix4nYfVWH5g?usp=drive_link): Markov Clustering Algorithm (MCL) script for clustering protein interactions
- [cluster_one-1.0.jar](https://drive.google.com/file/d/1v5qzD4x9a3Z_nhyYMqhfay43R5fTFfE1/view?usp=sharing): ClusterONE script for clustering protein interactions
- [srcData](https://drive.google.com/drive/folders/1uRBDZmEyVQ7tMJGwATt3HAvUCA5i3UjD?usp=sharing): Raw data for constructing the aggregate feature matrix
- [featureData](https://drive.google.com/drive/folders/1uKeFzRYbmWqB-Uj12LWMzTCsoh9X7hFd?usp=drive_link): Generated features from protein- association and expression data
- [dataPartitions](https://drive.google.com/drive/folders/1uKAA-lPKHZ16NfGQ3Daff1_SpDXyXIJ3?usp=drive_link): Five-fold partitioned datasets generated for cross-validation training and testing
- [models](https://drive.google.com/drive/folders/1vluMdzWFwDj3amFcsAuWMODswQzPxXtA?usp=drive_link): Trained Random Forest SciKit-Learn models, pickle files
- [modelPerformance](https://drive.google.com/drive/folders/1cv9KvPHR-Ki8yzAqYOqYa31JlGSX3ER0?usp=drive_link): Model evaluation metrics, e.g., precision, recall, and F1 scores
- [predictions](https://drive.google.com/drive/folders/1uGkJ9nt4os_8gi20hgkvY8gj3sh8KKZC?usp=drive_link): Scores for pairwise interactions
- [assemblyData](https://drive.google.com/drive/folders/1157Db3XnN_cW281r5jzjS49BSqDJ8XP7?usp=drive_link): Intermediate files generated to facilitate the prediction of complexes
- [complexPredictions](https://drive.google.com/drive/folders/1toyvOccqIr7IkaQoGo4ZDnpq5UZJf0N2?usp=drive_link): Predicted protein complexes
- [manuscriptGraphics](https://drive.google.com/drive/folders/1ty2RRhuJV2BThuRYmCwO9xOXOvvq8GDh?usp=drive_link): Results: manuscript figures and CSV files (tables)

### Running the pipeline
See [Using Google Colab with GitHub](https://colab.research.google.com/github/googlecolab/colabtools/blob/main/notebooks/colab-github-demo.ipynb) to learn how to load and run public notebooks directly from GitHub.

*Notes*
- *Google provides resources for free, however, runtimes are neither guaranteed nor unlimited; additionally, freely available runtimes may lack the specific configurations necessary for running some of the notebooks due to high memory and computational parallel processing demands.  The purchase of a plan is likely required.*

- *The root directory is set to /content/drive/My Drive.  Please ensure that this directory is correctly set to access your data.  Data storage is temporary; therefore, results should be saved before disconnecting from Colab or terminating a session.*

_Directions for generating specific results/output files_:
1) Copy the GitHub repository directory to either your local drive or preferred cloud storage.  
2) Refer above to *Directory Organization* and copy the input files associated with the desired results as specified below from the complete, publicly available directory.  
3) Finally, select *Run all* from the *Runtime* menu located at the top left of the Google Colab notebook. 

_Directions for running the pipeline from beginning to end_:
1) Copy the GitHub repository directory to either your local drive or preferred cloud storage. 
2) For each notebook in the order shown below, select *Run all* from the *Runtime* menu located at the top left of the Google Colab notebook. 

**1. Construct feature matrix.ipynb**
- **Purpose**: Assemble pairsets, import and generate features from original data sources
- **Inputs**: None
- **Outputs**:
  - Feature matrices, generated for each data source: `featureData/featMat_*tsv.tar.bz2`
  - Aggregated feature matrix, labeled: `featureData/featMat_aggLabeled.tsv.tar.bz2`
  - Aggregated feature matrix, unlabeled: `featureData/featMat_aggUnlabeled.tsv.tar.bz2`

**2. Build models.ipynb**
- **Purpose**: Train Random Forests (parallel)
- **Inputs**: `featureData/featMat_aggLabeled.tsv.tar.bz2`
- **Outputs**:
  - 5CV Data Partitions: `dataPartitions/*/5CV/*/trainingData.tsv`, `trainingLabels.tsv`, `testData.tsv`, `testLabels.tsv`
  - Models: `models/models*5CV*/model*pkl`
  - Training Metrics: `models/models*5CV*/model*selectFeatures_names.txt`, `Partition Specific*Training Performance_*csv` 

**3. Make predictions.ipynb**
- **Purpose**: Generates interaction scores for protein pairs.
- **Inputs**:
  - 5CV Data Partitions: `dataPartitions/*/5CV/*/trainingData.tsv`, `trainingLabels.tsv`, `testData.tsv`, `testLabels.tsv`
  - Models: `models/models*5CV*/model*pkl`
  - Training Metrics: `models/models*5CV*/model*selectFeatures_names.txt`, `Partition Specific*Training Performance_*csv`
  - Feature matrix: `featMat_aggUnlabeled.tsv.tar.bz2`
- **Outputs**:
  - Testing Metrics: `modelPerformance/modelsPerformance*/Partition-Classifier_Testing_Performance*csv`
  - Partitions Available (ea. pair): `modelPerformance/modelsPerformance*/candidatePartitions*csv`
  - Partition Predictions: `modelPerformance/modelsPerformance*/candidatePartitions_preds*csv`
  - Partition Scores: `modelPerformance/modelsPerformance*/candidatePartitions_predsProbs*csv`
  - Prec, recall: `modelPerformance/modelsPerformance*/recall_weighteds*csv`, `precision_weighted*.tsv`, `combined-weighted_results*.csv`

**4. Make predictions-CellSpecific.ipynb**
- **Purpose**: Generates interaction scores using pre-trained models and cell line-specific features.
- **Inputs**:
  - Feature matrices: `featureData/featMat_aggLabeled.tsv.tar.bz2`, `featMat_aggUnlabeled.tsv.tar.bz2`
  - 5CV Data Partitions: `dataPartitions/*/5CV/*/trainingData.tsv`, `trainingLabels.tsv`, `testData.tsv`, `testLabels.tsv`
  - Models: `models/models*5CV*/model*pkl`
  - Training Metrics: `models/models*5CV*/model*selectFeatures_names.txt`, `Partition Specific*Training Performance_*csv`
- **Outputs**:
  - NCI 60 Cell-line specific features: `featureData/cell_specific/nci60Feats_cellSpec_rscIncl.tsv.tar.bz2`
  - Non- and cell line -specific scores: `predictions/predictions.tsv.tar.bz2`

**5. Generate complexes.ipynb**
- **Purpose**: Predicts protein complexes by clustering pairs based on interaction probabilities
- **Inputs**: `predictions/predictions.tsv.tar.bz2`
- **Outputs**: `complexPredictions/*predictions*[nodes]*[size]*[density]*[overlap]*txt`

**6. Evaluate complexes helper.ipynb**
- **Purpose**: Prepares data for evaluating protein complex predictions.
- **Inputs**:
  - `predictions/predictions.tsv.tar.bz2`
  - `complexPredictions/*predictions*[nodes]*[size]*[density]*[overlap]*txt`
- **Outputs**:
  - (multiple param-setting cluster predictions consolidated): 
  `complexPredictions/clusterONE_MCL_predictions_gridSearch*.pkl`

**7. Analyze Results.ipynb**
- **Purpose**: Analyzes and presents results, including figures and tables.
- **Inputs**:
  - Prec, recall: `modelPerformance/modelsPerformance*/recall_weighteds*csv`, `precision_weighted*.tsv`, `combined-weighted_results*.csv`
  - Non- and cell line -specific scores: `predictions/predictions.tsv.tar.bz2`
  - NCI 60 Cell-line specific features: `featureData/cell_specific/nci60Feats_cellSpec_rscIncl.tsv.tar.bz2`
  - (multiple param-setting cluster predictions consolidated): 
  `complexPredictions/clusterONE_MCL_predictions_gridSearch*.pkl`
- **Outputs**:
  - Final results: `manuscriptGraphics/*`

### References
**Data and Software**:
- Protein Association Data:
    - [BioPlex 3.0](https://bioplex.hms.harvard.edu/data/) (Huttlin et al., 2017)
    - [hu.MAP 1.0](http://hu.proteincomplexes.org/download/) (Drew et al., 2017)
    - [hu.MAP 2.0](http://humap2.proteincomplexes.org/download) (Huttlin et al., 2021)


- Protein Expression Data:
    - [FANTOM5](https://www.proteinatlas.org/download/rna_tissue_fantom.tsv.zip) (Lizio et al., 2015)
    - [GTEx (v7)](https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip) (GTEx Consortium, 2017)
    - [NCI-60 Cell-Line Panel](https://www.proteomicsdb.org/proteomicsdb/#projects/35) (Scherf et al., 2000)
    - [SubCellBarCode](https://lehtio-lab.se/subcellbarcode/) (Orre et al., 2019)
    - [Uhlén et al., 2015](https://www.proteinatlas.org/download/rna_transcript.tsv.gz)


- Code:
    - [ClusterOne](https://paccanarolab.org/cluster-one/) (T. Nepusz et al., 2012)
    - [hu.MAP 2.0](https://github.com/marcottelab/protein_complex_maps) (Huttlin et al., 2021)
    - [MCL -- Markov Clustering Algorithm](https://micans.org/mcl/) (S. Van Dongen, 2008)
    - [Ouyang et al., 2019](https://modelzoo.cellprofiling.org)