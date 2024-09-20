# ELCFS Pipeline Readme
Scripts for generating manuscript results

## Overview
This pipeline predicts protein complexes by building protein-protein interaction models and assigning scores to as many pairs as possible from available expression or association data.  It then creates weighted graphs from which it identifies “clusters” with high node connectivity and subgraph cohesion.  The pipeline features a set of Jupyter notebooks that utilize “Google Colab” to harness Google’s powerful computing resources, also seamlessly integrated with GitHub, further benefiting greater accessibility.

## Running the pipeline
The results reflect the performance of a voting classifier committee that followed the manuscript’s eponymous approach on different pair sets as its feature-based composition expanded to accommodate more diverse sources of protein association and expression data.

Pair sets:
- D17 (Drew et al., 2017 training and test sets merged)
- LM19 (Lugo-Martinez et al., 2019 expanded pairs’ set)
- PS24 (CORUM 3.0 human core set-generated)

Protein Association Data: 
- BioPlex 3.0 
- hu.MAP 1.0, 2.0

Protein Expression Data:
- FANTOM5
- GTEx (v7)
- Lugo-Martinez et al., 2019 & 2021 
- Ouyang et al., 2019 
- SubCellBarCode 
- Uhlén et al., 2015

### Pairsets D17 & LM19
1. Use the following notebook to model interactions and generate pair scores:
- ./google_colab_notebooks/PPI Prediction Performance Analysis v3.ipynb

Note: For faster execution, n multiple instances of the notebook were generated, such that the parallel processing flags to the model training function were assigned to each integer between 0 and n-1.
 
2. Use the following notebook to predict complex assemblies:
- ./google_colab_notebooks/Generate complexes.ipynb

Note: Initially, threshold values along the interaction model’s F1-measure curve were selected manually to balance accuracy with processing time.

### Pairsets PS24
3. Use the following notebook to assemble the training data:
- ./google_colab_notebooks/Construct feature matrix.ipynb

4. Use the following notebook to create the models:
- ./google_colab_notebooks/Build models.ipynb

5. Use the following notebook to generate predictions: 
- ./google_colab_notebooks/Make predictions.ipynb

6. Use the following notebook to generate complexes:
- ./google_colab_notebooks/Generate complexes.ipynb

7. Use the following notebook to generate precursors for complex evaluation:
- ./google_colab_notebooks/Evaluate complexes helper.ipynb

### Analysis
8. Run the following notebook to generate results, e.g., figures and tables, across the pair sets
- ./google_colab_notebooks/Analyze results.ipynb