# Quantifying Pairwise Chemical Similarity for Polymers


This repository supports the following manuscript, which has been submitted for peer-review.

Jiale Shi, Nathan J. Rebello, Dylan Walsh, Weizhong Zou, Michael E. Deagen, Bruno Salomao Leao, Debra J. Audus, Bradley D. Olsen, "Quantifying Pairwise Chemical Similarity for Polymers", submitted.

In this work, we proposed a reliable method to quantitatively calculate the pairwise chemical similarity score for polymers, where the [earth mover’s distance (EMD)](https://en.wikipedia.org/wiki/Earth_mover%27s_distance) is utilized to calculate the similarity of the repeat units and end groups, while the [graph edit distance (GED)](https://en.wikipedia.org/wiki/Graph_edit_distance) is used to calculate the similarity of the topology. These three values then are combined to yield an overall pairwise chemical similarity score for polymers.

The repository is intended for the following use cases:

- Illustrate key ideas from the manuscript Method section including earth mover's distance and graph edit distance
- Allow for full reproducibility of the data in the manuscript

## Running the code

### Running notebooks in Google Colab

If you are interested in running one or more notebooks in [Google Colab](https://colab.research.google.com/), first click on the relevant link below.

#### Notebook for main text
* [Polymer_Similarity_Method_Section_PolymerA_PolymerB](./notebook/Polymer_Similarity_Method_Section_PolymerA_PolymerB.ipynb)
* [Polymer_Similarity_Case1_Varying_Repeat_Units](./notebook/Polymer_Similarity_Case1_Varying_Repeat_Units.ipynb)
* [Polymer_Similarity_Case2_Varying_Topologies](./notebook/Polymer_Similarity_Case2_Varying_Topologies.ipynb)
* [Polymer_Similarity_Case3_Varying_Both_Repeat_Units_and_Topologies](./notebook/Polymer_Similarity_Case3_Varying_Both_Repeat_Units_and_Topologies.ipynb)
* [Polymer_Similarity_Case4_Graft_Copolymers](./notebook/Polymer_Similarity_Case4_Graft_Copolymers.ipynb) 
* [Polymer_Similarity_Case5_Segmented_Polymers](./notebook/Polymer_Similarity_Case5_Segmented_Polymers.ipynb) 

#### Notebook for supporting information
* [SI-I: Earth Mover’s Distance vs. Simple Sum or Average Method](./notebook/Polymer_Similarity_SI_Earth_Mover_Distance_vs_Simple_Sum_or_Average_Method.ipynb)
* [SI-II: Modifying Weights for Final Ranking Order of  Overall Similarity Score](./notebook/Polymer_Similarity_SI_Modifying_Weights_for_Final_Ranking_Order_of_Overall_Similarity_Score.ipynb)
* SI-III: Overall Pairwise Similarity Scores with Arithmetic Means
  1. [Overall Pairwise Similarity Scores with Arithmetic Means Case 1](./notebook/Polymer_Similarity_SI_Overall_Pairwise_Similarity_Scores_with_Arithmetic_Means_Case1.ipynb)
  2. [Polymer_Similarity_Case2_Varying_Topologies_Arithmetic_Mean](./notebook/Polymer_Similarity_Case2_Varying_Topologies_Arithmetic_Mean.ipynb)
  3. [Polymer_Similarity_Case3_Varying_Both_Repeat_Units_and_Topologies_Arithmetic_Mean](./notebook/Polymer_Similarity_Case3_Varying_Both_Repeat_Units_and_Topologies_Arithmetic_Mean.ipynb)
  4. [Polymer_Similarity_Case4_Graft_Copolymers_Arithmetic_Mean](./notebook/Polymer_Similarity_Case4_Graft_Copolymers_Arithmetic_Mean.ipynb) 
  5. [Polymer_Similarity_Case5_Segmented_Polymers_Arithmetic_Mean](./notebook/Polymer_Similarity_Case5_Segmented_Polymers_Arithmetic_Mean.ipynb) 
* [SI-IV: Similarity Calculation for Random Copolymer Polymer A with Different Composition Ratio of Repeat Units and Diblock Copolymer Polymer B with Different Block Length Ratio](./notebook/Polymer_Similarity_PolymerA_PolymerB_Fine_Calculation.ipynb)
* [SI-V: Similarity Calculation for Graft Polymers with Different Lengths of Side Chains](./notebook/Polymer_Similarity_Case4_Graft_Copolymers_Fine_Calculation.ipynb)
* [SI-VI: Similarity Calculation for Segmented Polymers with Different Fractions and Different Lengths of Macromonomers](./notebook/Polymer_Similarity_Case5_Segmented_Polymers_Fine_Calculation.ipynb)
* [SI-VII: Similarity Calculation including Tacticity](./notebook/Polymer_Similarity_Tacticity.ipynb)

Then open the colab badge <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" width="75" height="15"/> in the notebook.

It will open a colab notebook. Then you can run the notebook as normal. All the required libraries and functions are present in the colab notebook.



## Contact

Jiale Shi, PhD  

Postdoctoral Associate  

Department of Chemical Engineering 

Massachusetts Institute of Technology (MIT) 

Email: jialeshi@mit.edu  

GithubID: shijiale0609  
 

## How to cite

If you use the code, please cite our repository since our manuscript is currently in review:

Jiale Shi, Nathan J. Rebello, Dylan Walsh, Weizhong Zou, Michael E. Deagen, Bruno Salomao Leao, Debra J. Audus, Bradley D. Olsen (2023), "Quantifying Pairwise Chemical Similarity for Polymers", https://github.com/olsenlabmit/Polymer-Graph-Similarity.
