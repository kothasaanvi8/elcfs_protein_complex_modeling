# README

### Overview
This pipeline predicts protein complexes by modeling protein-protein interactions and assigning scores to various
protein pairs using available association and expression data. It constructs weighted graphs and then identifies
clusters characterized by high node connectivity and subgraph cohesion.

### Key files of interest
- [pair probability scores](https://drive.google.com/file/d/1Cn5zRitrl4i1S_WU2DBKCpR1Aj0_bLmb/view?usp=sharing)
- [predicted clusters](https://drive.google.com/file/d/1-INlBd1bCZiYBd7vwJqvjYiBB3QpEoNK/view?usp=sharing)
- [Entrez to UniProt ID mapper](https://drive.google.com/file/d/1tT28krOYbgxQ2S3RtSvuA6cG2SA91AN0/view?usp=sharing)
- [NCI 60 non- and cell-specific features](https://drive.google.com/file/d/1-A66udgX7YWC-Oh7luKY5SzJuhKjPDKR/view?usp=sharing)

### Directory organization:
- **`utils/`**: Python scripts
- **`google_colab_notebooks/`**: Jupyter notebooks for running the pipeline
- **`CFinder-2.0.6--1448/`**: Markov Clustering Algorithm (MCL) for clustering protein interactions
- **`srcData/`**: Source data for building the feature matrix
- **`featureData/`**: Features generated from protein- association and expression data
- **`dataPartitions/`**: Five-fold partitioned datasets generated for cross-validation training and testing
- **`models/`**: Trained models
- **`modelPerformance/`**: Model evaluation metrics, e.g., precision, recall, and F1 scores for models
- **`predictions/`**: Pairwise interaction scores
- **`assemblyData/`**: Intermediate data generated to facilitate and evaluate protein complex predictions
- **`complexPredictions/`**: Protein complex predictions
- **`manuscriptGraphics/`**: Final results, including figures and CSV files from manuscript

### Run the pipeline (Notebooks, input and output files)
- See [Using Google Colab with GitHub](https://colab.research.google.com/github/googlecolab/colabtools/blob/main/notebooks/colab-github-demo.ipynb) to learn how to load and run public notebooks directly from GitHub
- The root directory is set to `/content/drive/My Drive/`. Please ensure that this directory is correctly set to access your data.
- Colab offers resources for free, but runtimes are neither guaranteed nor unlimited, so purchasing a plan may be necessary.
- Data storage is temporary; therefore, results should be saved before disconnecting from Colab or terminating a session.

For each notebook, select *Run all* from *Runtime* in the top menu. To generate specific results, use the input and output files available in the [pubic Google directory](https://drive.google.com/drive/folders/1OEa0NVYCYnhOZmRa_pycIlpNx9ySJIea?usp=drive_link) as specified below for the notebook of interest.

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