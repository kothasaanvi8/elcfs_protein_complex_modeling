# ELCFS Pipeline Readme

## <u>Overview</u>
This pipeline predicts protein complexes by building protein-protein interaction models and assigning scores to as many pairs as possible from available expression or association data.  It then creates weighted graphs from which it identifies “clusters” with high node connectivity and subgraph cohesion.

The results show the performance of a voting classifier committee that followed the manuscript’s eponymous approach on different pair sets as its feature-based composition expanded to accommodate more diverse sources of protein association and expression data.

*Pair sets:*
- D17 (Drew et al., 2017 training and test sets merged)
- LM19 (Lugo-Martinez et al., 2019 expanded pairs’ set)
- PS24 (CORUM 3.0 human core set-generated)

*Protein Association Data:* 
- BioPlex 3.0 
- hu.MAP 1.0, 2.0

*Protein Expression Data:*
- FANTOM5
- GTEx (v7)
- Lugo-Martinez et al., 2019 & 2021 
- Ouyang et al., 2019 
- SubCellBarCode 
- Uhlén et al., 2015

## <u>Resources for cloud computing</u> 
The pipeline includes a series of Jupyter notebooks located in the ‘google_colab_notebooks’ directory, which leverage Google Colab to access Google’s powerful computing resources.  While these notebooks can be used independently, cloud computing is much better suited to the computational demands.
Moreover, Google Colab is seamlessly integrated with GitHub, further benefiting greater accessibility.  Refer to the following for instructions on opening a Google Colab Notebook from GitHub.
1. Locate the GitHub repository.
2. Identify the notebook file and copy its URL
3. Open Google Colab (https://colab.research.google.com/)
4. Select open notebook under the File menu and click the GitHub tab.
5. Paste the URL of the notebook into the search bar (using the GitHub username/repo identifier is also an option) 
6. A list of notebooks will be displayed for the specified repository from which you may choose to both edit and run the notebook.
Note: Any changes must be saved. Furthermore, resources are not guaranteed at the free use level.

## <u>Run the pipeline</u>
*Assemble pair sets, pool data sources and construct the feature matrix.*
1.  Construct feature matrix.ipynb

*Generate RF models.*
2.  Build models.ipynb

*Generate scores for protein pairs.*
3.  Make predictions.ipynb

*Generate cell-line specific scores predictions for protein pairs.*
4.  Make predictions-CellSpecific.ipynb

*Generate protein complex predictions.*
5.  Generate complexes.ipynb

*Generate precursors to evaluation of protein complex predictions.*
6.  Evaluate complexes helper.ipynb

*Run the following notebook to generate results, e.g., figures and tables, across the pair sets.*
7.  Analyze Results.ipynb