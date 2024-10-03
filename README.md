ELCFS Pipeline - README

Overview

The ELCFS Pipeline is designed to predict protein complexes by modeling protein-protein interactions and assigning
scores to various protein pairs using available association and expression data. The pipeline constructs weighted graphs
and identifies clusters characterized by high node connectivity and subgraph cohesion.

The pipeline employs a voting classifier committee approach, as described in the manuscript, to evaluate performance
across different pair sets. As the feature composition grows, the pipeline integrates a wider array of protein
association and expression data sources.

Protein Association Data:

	•	BioPlex 3.0 (Huttlin et al., 2017)
	•	hu.MAP 1.0, 2.0 (Drew et al., 2017; Huttlin et al., 2021)

Protein Expression Data:

	•	FANTOM5 (Lizio et al., 2015)
	•	GTEx (v7) (GTEx Consortium, 2017)
	•	Lugo-Martinez et al., 2019 & 2021 (Lugo-Martinez et al., 2019; 2021)
	•	Ouyang et al., 2019 (Ouyang et al., 2019)
	•	SubCellBarCode (Baran-Gale et al., 2020)
	•	Uhlén et al., 2015 (Uhlén et al., 2015)

Google Colab Instructions

The pipeline includes several Jupyter notebooks located in the google_colab_notebooks directory, which use Google Colab
for computational tasks. Google Colab integrates with GitHub for seamless access to resources. Although these notebooks
can be run independently, cloud computing is ideal for handling the pipeline’s computational needs.

Steps to Run the Pipeline in Google Colab:
	1.	Open Google Colab @ "https://colab.research.google.com/"
	2.	In the File menu, select “Open Notebook” 
    3.  On the subsequent GUI that pops up, select the "GitHub" tab from the sidebar on the left.
	4.	In the GitHub search bar, enter the following repository: GaryWilkins/elcfs_protein_complex_modeling.
	5.	Run the notebooks in order as they are numbered. Descriptions of the notebooks are provided below.

	Notes: 
    i)  Google Colab does not provide persistent storage outputs, thus before disconnecting the runtime, outputs must be saved either locally or to the user's cloud storage, e.g., Google Drivee. 
    ii) Google Colab’s free-tier resources are limited; you may need a subscription for more intensive computations.

Notebooks:
	1.	Construct feature matrix.ipynb
            - Combines pair sets, pools data sources, and constructs the feature matrix.
	2.	Build models.ipynb
            - Trains random forest (RF) models for protein interactions.
	3.	Make predictions.ipynb
            - Generates interaction scores for protein pairs.
	4.	Make predictions-CellSpecific.ipynb
            - Generates cell-specific interaction scores for protein pairs.
	5.	Generate complexes.ipynb
            - Predicts protein complexes based on interaction scores.
	6.	Evaluate complexes helper.ipynb
            - Prepares data for evaluating protein complex predictions.
	7.	Analyze Results.ipynb
            - Generates results, including figures and tables, for the pair sets.