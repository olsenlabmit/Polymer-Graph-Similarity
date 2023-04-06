# Quantifying Pairwise Chemical Similarity for Polymers


This repository supports the following manuscript which has been submitted for peer-review.

Jiale Shi, Nathan J. Rebello, Dylan Walsh, Weizhong Zou, Michael E. Deagen, Bruno Salomao Leao, Debra J. Audus, Bradley D. Olsen, "Quantifying Pairwise Chemical Similarity for Polymers", *Macromolecule*, submitted.

Defining the similarity between chemical entities is an essential task in polymer informatics, enabling ranking, clustering, and classification. Despite its importance, pairwise chemical similarity for polymers remains an open problem. Here, a similarity function for polymers with well-defined backbones is designed based on polymers’ stochastic graph representations generated from canonical BigSMILES, a structurally-based line notation for describing macromolecules. The stochastic graph representations are separated into three parts: repeat units, end groups, and polymer topology. The [earth mover’s distance](https://en.wikipedia.org/wiki/Earth_mover%27s_distance) is utilized to calculate the similarity of the repeat units and end groups, while the [graph edit distance](https://en.wikipedia.org/wiki/Graph_edit_distance) is used to calculate the similarity of the topology. These three values can be linearly or nonlinearly combined to yield an overall pairwise chemical similarity score for polymers that is largely consistent with the chemical intuition of expert users and is adjustable based on the relative importance of different chemical features for a given similarity problem. This method gives a reliable solution to quantitatively calculate the pairwise chemical similarity score for polymers and represents a vital step toward building search engines and quantitative design tools for polymer data.

The repository is intended for the following use cases:

- Illustrate key ideas from the manuscript Method section including earth mover's distance and graph edit distance
- Allow for full reproducibility of the data in the manuscript

## Running the code

### Running notebooks in Google Colab

If you are interested in running one or more notebooks in [Google Colab](https://colab.research.google.com/), first click on the relevant link below.

#### Notebook
[Polymer_Similarity_Case1](./notebook/Polymer_Similarity_Case1.ipynb) [available]
[Polymer_Similarity_Case2]() [not available]
[Polymer_Similarity_Case3]() [not available]
[Polymer_Similarity_Case4]() [not available]
[Polymer_Similarity_Case5]() [not available]
[Polymer_Similarity_Case6]() [not available]

Then open the colab badge <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" width="150" height="30"/> in the notebook.

It will open a colab notebook. Then you can run the notebook as normal. All the required libraries and functions exist in the colab notebook.



## Contact

Jiale Shi, PhD  
Postdoctoral Associate  
Department of Chemical Engineering 
Massachusetts Institute of Technology (MIT) 

Email: jialeshi@mit.edu  
GithubID: shijiale0609  
 

## How to cite

If you use the code, please cite our repository since our manuscript is currently submitted for peer-review:

Jiale Shi, Nathan J. Rebello, Dylan Walsh, Weizhong Zou, Michael E. Deagen, Bruno Salomao Leao, Debra J. Audus, Bradley D. Olsen (2023), "Quantifying Pairwise Chemical Similarity for Polymers", https://github.com/olsenlabmit/Polymer-Graph-Similarity
