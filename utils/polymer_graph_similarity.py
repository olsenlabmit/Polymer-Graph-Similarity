import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import shutil
import sys
import os.path

#necessary library
#!pip install -q rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs

#if not shutil.which("pyomo"):
#    !pip install -q pyomo
#    assert(shutil.which("pyomo"))

#if not (shutil.which("cbc") or os.path.isfile("cbc")):
#    if "google.colab" in sys.modules:
#        !apt-get install -y -qq coinor-cbc
#        print("install cbc")
#    else:
#        try:
#            !conda install -c conda-forge coincbc 
#            print("install cbc-2")
#        except:
#            print("skip cbc")
#            pass
            

assert(shutil.which("cbc") or os.path.isfile("cbc"))
    
from pyomo.environ import *


#from pyomo.common.fileutils import Executable
def Similarity_Score_EMD(query_smiles_list = None, 
                         query_smiles_level_list = None, 
                         target_smiles_list = None, 
                         target_smiles_level_list = None,
                         #level_weight = False,
                         #level_ratio = 3,
                         embedding_function = 'MorganFingerprint', #
                         similarity_score_function = 'Tanimoto', # Dice, Cosine
                         restrain_emd = False):
  
    #obtain the length of query smiles list and target smiles list
    if query_smiles_list != None:
        query_smiles_list_length = len(query_smiles_list)
    else:
        print ("Missing query smiles list")
        return    

    if target_smiles_list != None:
        target_smiles_list_length = len(target_smiles_list)
    else:
        print ("Missing target smiles list")
        return
    
    if set(query_smiles_list) == set(target_smiles_list):
        query_smiles_array = np.array(query_smiles_list)
        query_smiles_level_array = np.array(query_smiles_level_list)
        inds_query = query_smiles_array.argsort()
        sorted_query_smiles_array = query_smiles_array[inds_query]
        sorted_query_smiles_level_array = query_smiles_level_array[inds_query]

        target_smiles_array = np.array(target_smiles_list)
        target_smiles_level_array = np.array(target_smiles_level_list)
        inds_target = target_smiles_array.argsort()
        sorted_target_smiles_array = target_smiles_array[inds_target]
        sorted_target_smiles_level_array = target_smiles_level_array[inds_target]

        if np.array_equal(sorted_query_smiles_array, sorted_target_smiles_array) and np.array_equal(sorted_query_smiles_level_array, sorted_target_smiles_level_array):
            return 1

        query_smiles_reduced_list = list(set(query_smiles_list))
        query_smiles_reduced_list_number = []
        for i in range(0, len(query_smiles_reduced_list)):
            query_smiles_reduced_list_number_i  = 0
            for j in range(0, len(query_smiles_list)):
                if query_smiles_reduced_list[i] == query_smiles_list[j]:
                    query_smiles_reduced_list_number_i = query_smiles_reduced_list_number_i   + query_smiles_level_list[j]
            query_smiles_reduced_list_number.append(query_smiles_reduced_list_number_i/sum(query_smiles_level_list))

        target_smiles_reduced_list = list(set(target_smiles_list))
        target_smiles_reduced_list_number = []
        for i in range(0, len(target_smiles_reduced_list)):
            target_smiles_reduced_list_number_i  = 0
            for j in range(0, len(target_smiles_list)):
                if target_smiles_reduced_list[i] == target_smiles_list[j]:
                    target_smiles_reduced_list_number_i = target_smiles_reduced_list_number_i   + target_smiles_level_list[j]
            target_smiles_reduced_list_number.append(target_smiles_reduced_list_number_i/sum(target_smiles_level_list))

        if query_smiles_reduced_list_number == target_smiles_reduced_list_number:
            #print("check3")
            return 1.0

    if embedding_function == 'RDKFingerprint':    
        query_mol_list = [Chem.MolFromSmiles(x) for x in query_smiles_list]
        query_fingerprint_list = [Chem.RDKFingerprint(x) for x in query_mol_list]
        target_mol_list = [Chem.MolFromSmiles(x) for x in target_smiles_list]
        target_fingerprint_list = [Chem.RDKFingerprint(x) for x in target_mol_list]

    elif embedding_function == 'MorganFingerprint':    
        query_mol_list = [Chem.MolFromSmiles(x) for x in query_smiles_list]
        query_fingerprint_list = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=2048) for x in query_mol_list]
        target_mol_list = [Chem.MolFromSmiles(x) for x in target_smiles_list]
        target_fingerprint_list = [AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=2048) for x in target_mol_list]

    elif embedding_function == 'MACCSkeys':    
        query_mol_list = [Chem.MolFromSmiles(x) for x in query_smiles_list]
        query_fingerprint_list = [MACCSkeys.GenMACCSKeys(x) for x in query_mol_list]
        target_mol_list = [Chem.MolFromSmiles(x) for x in target_smiles_list]
        target_fingerprint_list = [MACCSkeys.GenMACCSKeys(x) for x in target_mol_list]


    else:
        print(embedding_function + " is not included in the current vision, please choose an available embedding function.");
        return False         

    Demand = {}
    Supply = {}
    T = {}

    if query_smiles_level_list == None or len(set(query_smiles_level_list)) ==1:
        for i in range(0, query_smiles_list_length):
            Demand["P" + str(i+1)] = 1/query_smiles_list_length
    else:
        print("Query smiles list has different levels")
        #query_weight_sum = sum(query_smiles_level_list)
        #level 1, weight = 1; level 2, weight = 3^1;  level n, weight = 3^n-1
        query_weight_sum = 0.0
        for i in range(0, query_smiles_list_length):
            query_weight_sum = query_weight_sum + query_smiles_level_list[i]

        for i in range(0, query_smiles_list_length):
            Demand["P" + str(i+1)] = query_smiles_level_list[i]/query_weight_sum
            #print("P" + str(i+1), Demand["P" + str(i+1)])


    if target_smiles_level_list == None or len(set(target_smiles_level_list)) ==1:       
        for j in range(0,target_smiles_list_length):
            Supply["Q" + str(j+1)] = 1/target_smiles_list_length
    else:
        print("Target smiles list has different levels")
        
        target_weight_sum = 0.0
        for j in range(0, target_smiles_list_length):
            target_weight_sum = target_weight_sum + target_smiles_level_list[j]

        for j in range(0, target_smiles_list_length):
            Supply["Q" + str(j+1)] = target_smiles_level_list[j]/target_weight_sum
            #print("Q" + str(j+1), Supply["Q" + str(j+1)])


    # embedding function and similarity 
    if similarity_score_function == 'Tanimoto':
        for i in range(0,query_smiles_list_length):
            for j in range(0,target_smiles_list_length):
                # calculate the fingerprint similarityscore between query[i],target[j] and input the distance = 1- similarityscore
                T[("P" + str(i+1), "Q" + str(j+1))] = 1 - DataStructs.FingerprintSimilarity(query_fingerprint_list[i],target_fingerprint_list[j])
                print("P" + str(i+1), "->Q" + str(j+1), T[("P" + str(i+1), "Q" + str(j+1))] )

    elif similarity_score_function == 'Dice':
        for i in range(0,query_smiles_list_length):
            for j in range(0,target_smiles_list_length):
                # calculate the fingerprint similarityscore between query[i],target[j] and input the distance = 1- similarityscore
                T[("P" + str(i+1), "Q" + str(j+1))] = 1 - DataStructs.FingerprintSimilarity(query_fingerprint_list[i],target_fingerprint_list[j], metric=DataStructs.DiceSimilarity)
                #print("P" + str(i+1), "->Q" + str(j+1), T[("P" + str(i+1), "Q" + str(j+1))] )

    elif similarity_score_function == 'Cosine':
        for i in range(0,query_smiles_list_length):
            for j in range(0,target_smiles_list_length):
                # calculate the fingerprint similarityscore between query[i],target[j] and input the distance = 1- similarityscore
                T[("P" + str(i+1), "Q" + str(j+1))] = 1 - DataStructs.FingerprintSimilarity(query_fingerprint_list[i],target_fingerprint_list[j], metric=DataStructs.CosineSimilarity)
                #print("P" + str(i+1), "->Q" + str(j+1), T[("P" + str(i+1), "Q" + str(j+1))] )
       
    else:
        print(similarity_score_function + " is not included in the current vision, please choose an available similarity function.");
        return

    

    #print(len(Demand), len(Supply), len(T))
    # Step 0: Create an instance of the model
    model = ConcreteModel()
    model.dual = Suffix(direction=Suffix.IMPORT)

    # Step 1: Define index sets
    CUS = list(Demand.keys())
    SRC = list(Supply.keys())

    # Step 2: Define the decision 
    model.x = Var(CUS, SRC, domain = NonNegativeReals)

    # Step 3: Define Objective
    model.Cost = Objective(
    expr = sum([T[c,s]*model.x[c,s] for c in CUS for s in SRC]),
    sense = minimize)

    # Step 4: Constraints
    model.src = ConstraintList()
    for s in SRC:
        model.src.add(sum([model.x[c,s] for c in CUS]) == Supply[s])
        
    model.dmd = ConstraintList()
    for c in CUS:
        model.dmd.add(sum([model.x[c,s] for s in SRC]) == Demand[c])

    # add restrain to the EMD
    if restrain_emd == True: 
        model.restrain = ConstraintList()
        for i in range(0,query_smiles_list_length):
            model.restrain.add(model.x[CUS[i],SRC[i]] == Supply["Q" + str(i+1)])
    
    results = SolverFactory('cbc').solve(model)


    if 'ok' == str(results.Solver.status):
        #print("EMD(P,Q) = ",model.Cost())
        #print ("\n")
        #print("S(P,Q) = ", 1- model.Cost())
        SimilarityScore = 1- model.Cost()
        return SimilarityScore
        
    else:
        print("No Valid Solution Found")
        return False





# Graph Edit Distance
def Similarity_Score_Graph_Edit_Distance(Graph1 = None, 
                         Graph2 = None, 
                         alpha = 1):
    if Graph1 == None:
        print("Missing Graph1")
        return
    if Graph2 == None:
        print("Missing Graph2")
        return
    
    # Since Graph1 is the subgraph of Graph2, the calculation of graph edit distance can be simpilified.
    Graph1_number = Graph1.number_of_nodes() 
    Graph2_number = Graph2.number_of_nodes() 
    #graph_edit_distance = abs(Graph2_number - Graph1_number)
    graph_edit_distance = nx.graph_edit_distance(Graph1, Graph2)
    print(Graph1_number,Graph2_number, graph_edit_distance )
    # utilize the exponential decay function to turn the graph edit distance to similarity score
    #similarity_score = np.exp(-alpha*graph_edit_distance/(min(Graph1_number, Graph2_number))
    similarity_score = np.exp(-alpha*graph_edit_distance/((Graph1_number+ Graph2_number)/2))

    return similarity_score 



def Combined_Similarity_Score(Repeat_Unit_Similarity_Score = None,
                              Repeat_Unit_Weight = 0.5,
                              Graph_Similarity_Score = None,
                              Graph_Weight = 0.5,
                              End_Group_Similarity_Score = None,
                              End_Group_Weight = 0.0,
                              Mean_Function = 'arithmetic'):
  
    # Verify whether the weight sum is normalized.
    if  abs(Repeat_Unit_Weight + Graph_Weight +  End_Group_Weight -1 ) >=0.000000000001:
        print("Weight Sum is not normalized.")
        return False

    # Not consider the end group
    if End_Group_Similarity_Score == None:

        if Mean_Function == 'arithmetic':
            combined_similarity_score = (Repeat_Unit_Weight * Repeat_Unit_Similarity_Score + Graph_Weight * Graph_Similarity_Score)

        elif Mean_Function == 'geometric':
            combined_similarity_score = pow(Repeat_Unit_Similarity_Score,Repeat_Unit_Weight)*pow(Graph_Similarity_Score,Graph_Weight)

        else:
            print("Your input mean function ", Mean_Function, " is not implemented, please choose those implemented mean function, like arithmetic, geometric")
    
    # consider the end group
    else:
      
        if Mean_Function == 'arithmetic':
            combined_similarity_score = Repeat_Unit_Weight * Repeat_Unit_Similarity_Score + Graph_Weight * Graph_Similarity_Score + End_Group_Weight * End_Group_Similarity_Score

        elif Mean_Function == 'geometric':
            combined_similarity_score = pow(Repeat_Unit_Similarity_Score,Repeat_Unit_Weight)*pow(Graph_Similarity_Score,Graph_Weight)*pow(End_Group_Similarity_Score, End_Group_Weight)

        else:
            print("Your input mean function ", Mean_Function, " is not implemented, please choose those implemented mean function, like arithmetic, geometric")
            
    return combined_similarity_score




def Similarity_Score_Two_Polymer(query = None,
                                 target = None,
                                 #level_weight = True,
                                 #level_ratio = 3,
                                 embedding_function = 'RDKFingerprint', #Embedding function
                                 similarity_score_function = 'Tanimoto', # Similarity function for two vectors
                                 restrain_emd = False, # Whether to restrain the emd
                                 alpha=1, #reduced parameter for the exponential decay function
                                 Repeat_Unit_Weight=0.5,
                                 Graph_Weight=0.5,
                                 End_Group_Weight = 0.0,
                                 Mean_Function = 'geometric',
                                 details_print = False):

    if query == None or target == None:
        print ("Either query polymer or target polymer is missing! Please check your input.")
        return False
    S_repeat_unit = Similarity_Score_EMD(query_smiles_list = query.repeat_unit_smiles_list, 
                     query_smiles_level_list = query.repeat_unit_smiles_level_list, 
                     target_smiles_list = target.repeat_unit_smiles_list, 
                     target_smiles_level_list = target.repeat_unit_smiles_level_list,
                     #level_weight = level_weight,
                     #level_ratio = level_ratio,	
                     embedding_function = embedding_function,
                     similarity_score_function = similarity_score_function,
                     restrain_emd = restrain_emd)
    
    S_graph = Similarity_Score_Graph_Edit_Distance(Graph1=query.graph_representation, 
                                                   Graph2=target.graph_representation, 
                                                   alpha=alpha)
    
    if End_Group_Weight == 0.0: 
        S_combined = Combined_Similarity_Score(Repeat_Unit_Similarity_Score=S_repeat_unit ,
                                       Repeat_Unit_Weight=Repeat_Unit_Weight,
                                       Graph_Similarity_Score=S_graph,
                                       Graph_Weight=Graph_Weight,
                                       End_Group_Similarity_Score = None,
                                       End_Group_Weight = End_Group_Weight,
                                       Mean_Function = Mean_Function)
        if details_print == True:
            print("Details of the Similarity Score:\n")
            print("Similarity score on Repeating Unit = ", S_repeat_unit, ", Weight for Repeating Unit = ", Repeat_Unit_Weight)
            print("Similarity score on Graph = ", S_graph, ", Weight for Graph = ", Graph_Weight)
            print("Similarity score on End Group = ", "None", ", Weight for End Group = ", End_Group_Weight)
            print("Similarity score Combined in " + Mean_Function + " mean = ", S_combined)
            print("\n")

        return S_combined 

    else: 
        S_end_group = Similarity_Score_EMD(query_smiles_list = query.end_group_smiles_list, 
                     query_smiles_level_list = query.end_group_smiles_level_list, 
                     target_smiles_list = target.end_group_smiles_list, 
                     target_smiles_level_list = target.end_group_smiles_level_list,
                     #level_weight = level_weight,
                     #level_ratio = level_ratio,	
                     embedding_function = embedding_function,
                     similarity_score_function = similarity_score_function,
                     restrain_emd = restrain_emd)
            
        S_combined = Combined_Similarity_Score(Repeat_Unit_Similarity_Score=S_repeat_unit ,
                                       Repeat_Unit_Weight=Repeat_Unit_Weight,
                                       Graph_Similarity_Score=S_graph,
                                       Graph_Weight=Graph_Weight,
                                       End_Group_Similarity_Score = S_end_group,
                                       End_Group_Weight = End_Group_Weight,
                                       Mean_Function = Mean_Function)
        
        if details_print == True:
            print("Details of the Similarity Score:\n")
            print("Similarity score on Repeating Unit = ", S_repeat_unit, ", Weight for Repeating Unit = ", Repeat_Unit_Weight)
            print("Similarity score on Graph = ", S_graph, ", Weight for Graph = ", Graph_Weight)
            print("Similarity score on End Group = ", S_end_group, ", Weight for End Group = ", End_Group_Weight )
            print("Similarity score Combined in " + Mean_Function + " mean = ", S_combined)
            print("\n")
            return S_combined, S_repeat_unit, S_graph, S_end_group

        return S_combined   
    

class Polymer:
    def __init__(self, 
               repeat_unit_smiles_list=None, 
               repeat_unit_smiles_level_list=None, 
               end_group_smiles_list=None, 
               end_group_smiles_level_list=None,
               graph_representation=None):
    
        self.repeat_unit_smiles_list = repeat_unit_smiles_list
        self.repeat_unit_smiles_level_list = repeat_unit_smiles_level_list
        self.end_group_smiles_list = end_group_smiles_list
        self.end_group_smiles_level_list = end_group_smiles_level_list
        self.graph_representation = graph_representation


