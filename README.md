# ELCFS Pipeline Readme


## Overview


This pipeline predicts protein complexes by generating protein-protein interaction models and scoring all available protein pairs. To facilitate user-friendly cloud computing, it includes a library of Python scripts and a collection of Google Colab notebooks, which are also compatible with the Jupyter computing platform.


## Manuscript Results


The results presented in the manuscript were generated over two phases of pipeline development, each corresponding to different sets of proteins. Below are the primary sources of data used:


Protein Association Data: 
  - Drew et al., 2017
  - BioPlex 3.0 (Huttlin et al., 2021)


Protein Expression Data: 
  - Lugo-Martinez et al., 2019 & 2021
  - SubCellBarCode (Orre et al., 2019)
  - FANTOM5 (Ouyang et al., 2014)
  - Uhlén et al., 2015
  - gTEX (Yizhak et al., 2019)


## Pipeline Execution


### Preliminary Results


1. Model Interactions and Score Protein Pairs:
   - Use the PPI Prediction Performance Analysis v3.ipynb notebook to model interactions and score all available protein pairs.
   - For faster execution, enable parallel execution by making copies of the notebook and setting partition intervals for each copy.


   ./google_colab_notebooks/PPI Prediction Performance Analysis v3.ipynb


2. Predict Complex Assemblies:
   - Use the Generate complexes.ipynb notebook to predict complex assemblies.
   - Note: We manually selected thresholds between the interaction model’s F1-measure and values that provided manageable processing times.


   ./google_colab_notebooks/Generate complexes.ipynb


3. Prepare Precursors for Evaluation Functions:
   - Use the complexPrediction_plusAnalysis.ipynb notebook to prepare precursors for evaluation functions.


   ./google_colab_notebooks/complexPrediction_plusAnalysis.ipynb


### Final Results


1. Assemble Training Data:
   - Use the Construct feature matrix.ipynb notebook to assemble the training data.


     ./google_colab_notebooks/Construct feature matrix.ipynb


2. Create Models:
   - Use the Build models.ipynb notebook to create the models.
   - Note: During Phase 1a, we made copies of the modeling notebook for each fold of cross-validation to take advantage of parallel execution while minimizing deviations from total pipeline automation.


   ./google_colab_notebooks/Build models.ipynb


3. Generate Predictions:
   - Use the Make predictions.ipynb notebook to generate predictions.


   ./google_colab_notebooks/Make predictions.ipynb


4. Generate Complexes:
   - Use the Generate complexes.ipynb notebook to generate complexes.


   ./google_colab_notebooks/Generate complexes.ipynb


5. Generate Precursors for Complex Evaluation:
   - Use the Evaluate complexes helper.ipynb notebook to generate precursors for complex evaluation.


   ./google_colab_notebooks/Evaluate complexes helper.ipynb


### Analysis


- Run the Analyze results.ipynb notebook to collect results from both phases and generate the manuscript’s figures and tables.


  ./google_colab_notebooks/Analyze results.ipynb