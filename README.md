ELCFS Pipeline Readme

This pipeline predicts complexes for proteins by generating protein-interaction models and scoring all pairs for which data is available.  It includes a library of Python scripts and a collection of Google Colab notebooks, also available on the Jupyter computing platform, to facilitate user-friendly cloud computing.

We generated the manuscript results in 2 phases over pipeline development, each corresponding to different sets of proteins.  
Pair-level data from Drew et al., 2017, its protein-level feature expansion  from Lugo-Martinez et al., 2019, 2021 and protein-level data from using SubCellBarCode (Orre et al., 2019) 

(added) Pair-level data from BioPlex 3.0 (Huttlin et al., 2021) and protein-level data from using FANTOM5 (2014), Ouyang et al., 2019, Uhlén et al., 2015, and gTEX (Yizhak et al., 2019).


  Preliminary Results

Model interactions and score all available protein-pairs
For faster execution: enable parallel execution, make copies, and for each notebook copy, set the interval of partitions to the model.   
./google_colab_notebooks/PPI Prediction Performance Analysis v3.ipynb

Predict complex assemblies
Note: We manually intervened, selecting thresholds between the interactions’ model F1-measure and values for which processing time was manageable.
./google_colab_notebooks/Generate complexes.ipynb

Make precursors for evaluation functions
./google_colab_notebooks/complexPrediction_plusAnalysis.ipynb

  Final results

Assemble training data
./google_colab_notebooks/Construct feature matrix.ipynb

Create models
Note: See Phase 1a.  We made a copy of the modeling notebook only for each fold of cross-validation to continue taking advantage of parallel execution while minimizing deviation from total pipeline automation.
/google_colab_notebooks/Build models.ipynb

Generated predictions
./google_colab_notebooks/Make predictions.ipynb

Generated complexes 
./google_colab_notebooks/Generate complexes.ipynb

Generated precursors to complex evaluation for analysis
./google_colab_notebooks/Evaluate complexes helper.ipynb

Analysis

Run the analysis notebook to collect results from both phases and generate the manuscript’s figures and tables.
./google_colab_notebooks/Analyze results.ipynb
