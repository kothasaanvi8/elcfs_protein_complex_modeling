# README

## Overview
This pipeline predicts protein complexes by modeling protein-protein interactions and assigning scores to various protein pairs using available association and expression data. It constructs weighted graphs and identifies clusters with high node connectivity and subgraph cohesion. The repository contains all necessary files and instructions to reproduce these results. Due to the size of certain generated results, they are stored in this [Google Drive folder](https://drive.google.com/drive/folders/1OEa0NVYCYnhOZmRa_pycIlpNx9ySJIea?usp=drive_link).

### Key files in the Google Drive directory expected to be of use from this project:
- <span style="font-weight: normal;">[predictions.tsv.tar.bz2](https://drive.google.com/file/d/1Cn5zRitrl4i1S_WU2DBKCpR1Aj0_bLmb/view?usp=sharing)</span>: Predicted probabilities for protein pair interactions
- <span style="font-weight: normal;">[clustersELCFS_topParams.pkl](https://drive.google.com/file/d/1-INlBd1bCZiYBd7vwJqvjYiBB3QpEoNK/view?usp=sharing)</span>: Clustered groups of proteins predicted to form complexes
- <span style="font-weight: normal;">[idmappingDF_Entrez->UniprotKB.pkl](https://drive.google.com/file/d/1tT28krOYbgxQ2S3RtSvuA6cG2SA91AN0/view?usp=sharing)</span>: Mapping file to translate Entrez IDs to UniProt IDs
- <span style="font-weight: normal;">[nci60Feats_non+cellSpecific_resc.tsv.tar.bz2](https://drive.google.com/file/d/1-A66udgX7YWC-Oh7luKY5SzJuhKjPDKR/view?usp=sharing)</span>: Cell line- and non-specific comparative features generated from protein and RNA expression in the NCI-60 Human Tumor Cell Lines Screen

### How to run the pipeline
See [Using Google Colab with GitHub](https://colab.research.google.com/github/googlecolab/colabtools/blob/main/notebooks/colab-github-demo.ipynb)
to learn how to load and run public notebooks directly from GitHub.

To run the pipeline from start to finish, select *Run all* from the *Runtime* menu located at the top left of each each notebook listed in the order shown under Pipeline Execution.<br>

Notes:<br>
_*Google provides free resources, but runtimes are temporary, neither guaranteed nor unlimited._<br>
_*Certain configurations may require a paid plan._<br>
_*All data is cleared after disconnecting/releasing the runtime. Save data to either a local drive or preferred cloud storage._<br>
_*Directories absent from the repo, due to GitHub storage constraints, e.g., featureData, are created automatically during execution_.<br>

### Pipeline Execution
**i. Construct feature matrix.ipynb**
- **Purpose**: Assemble pairsets, import and generate features from original data sources
- **Inputs**: None
- **Outputs**:
  - Feature matrices, generated for each data source: `featureData/featMat_*tsv.tar.bz2`
  - Aggregated feature matrix, labeled (feature matrices fused with CORUM-generated labels): `featureData/featMat_aggLabeled.tsv.tar.bz2`
  - Aggregated feature matrix, unlabeled: `featureData/featMat_aggUnlabeled.tsv.tar.bz2`

**ii. Build models.ipynb**
- **Purpose**: Train Random Forests (parallel)
- **Inputs**: `featureData/featMat_aggLabeled.tsv.tar.bz2`
- **Outputs**:
  - 5CV Data Partitions: `dataPartitions/*/5CV/*/trainingData.tsv`, `trainingLabels.tsv`, `testData.tsv`, `testLabels.tsv`
  - Models: `models/models*5CV*/model*pkl`
  - Training Metrics: `models/models*5CV*/model*selectFeatures_names.txt`, `Partition Specific*Training Performance_*csv` 

**iii. Make predictions.ipynb**
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

**iv. Make predictions-CellSpecific.ipynb**
- **Purpose**: Generates interaction scores using pre-trained models and cell line-specific features.
- **Inputs**:
  - Feature matrices: `featureData/featMat_aggLabeled.tsv.tar.bz2`, `featMat_aggUnlabeled.tsv.tar.bz2`
  - 5CV Data Partitions: `dataPartitions/*/5CV/*/trainingData.tsv`, `trainingLabels.tsv`, `testData.tsv`, `testLabels.tsv`
  - Models: `models/models*5CV*/model*pkl`
  - Training Metrics: `models/models*5CV*/model*selectFeatures_names.txt`, `Partition Specific*Training Performance_*csv`
- **Outputs**:
  - NCI 60 Cell-line specific features: `featureData/cell_specific/nci60Feats_cellSpec_rscIncl.tsv.tar.bz2`
  - Non- and cell line -specific scores: `predictions/predictions.tsv.tar.bz2`

**v. Generate complexes.ipynb**
- **Purpose**: Predicts protein complexes by clustering pairs based on interaction probabilities
- **Inputs**: `predictions/predictions.tsv.tar.bz2`
- **Outputs**: `complexPredictions/*predictions*[nodes]*[size]*[density]*[overlap]*txt`

**vi. Evaluate complexes helper.ipynb**
- **Purpose**: Prepares data for evaluating protein complex predictions.
- **Inputs**:
  - `predictions/predictions.tsv.tar.bz2`
  - `complexPredictions/*predictions*[nodes]*[size]*[density]*[overlap]*txt`
- **Outputs**:
  - (multiple param-setting cluster predictions consolidated): 
  `complexPredictions/clusterONE_MCL_predictions_gridSearch*.pkl`

**vii. Analyze Results.ipynb**
- **Purpose**: Analyzes and presents results, including figures and tables.
- **Inputs**:
  - Prec, recall: `modelPerformance/modelsPerformance*/recall_weighteds*csv`, `precision_weighted*.tsv`, `combined-weighted_results*.csv`
  - Non- and cell line -specific scores: `predictions/predictions.tsv.tar.bz2`
  - NCI 60 Cell-line specific features: `featureData/cell_specific/nci60Feats_cellSpec_rscIncl.tsv.tar.bz2`
  - (multiple param-setting cluster predictions consolidated): 
  `complexPredictions/clusterONE_MCL_predictions_gridSearch*.pkl`
- **Outputs**:
  - Final results: `manuscriptGraphics/*`

### Alternatively, each notebook constitutes a module of the pipeline, thus enabling independent execution. To run a specific notebook, follow the subsequent steps.

1) Identify the notebook for which the desired files are the output.  
2) Copy the notebook's corresponding input files, and their directories, if absent from the aforementioned Google Drive folder to the directory specified under Pipeline Execution.  
3) Finally, select *Run all* from the *Runtime* menu located at the top left of the Google Colab notebook.


### Data and Software

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