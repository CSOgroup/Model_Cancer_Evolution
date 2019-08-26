
import numpy as np
import random
from ete2 import Tree
import time
random.seed(time.time())

OUTPUT_PATH="my_output_path" #define the output path

#evolutionnary parameters
mut_rate=0.0001
fitness=0.0001
fitness_d=0. #change to non zero value to have effect of deleterious passenger mutations

#simulation parameters
initsize=1 #initial number of cell in the first clone
maxcell=500000000 #final size f the tumor
maxgene=20000 #total number of genes
maxdriver=250 #change this value to increase the driver mutation rate
max_success=50 # number of successful simulations

nb_success=0

while nb_success<max_success:
    seed=random.randint(0,100000000)

    Tumor=Tree()
    founder_clone=Tumor.get_tree_root()
    founder_clone.name="founder"
    founder_clone.add_features(cell=1,driver=1,passenger=0,mutation="0\t") # number of cells in the node, number of driver in each cell, number of passengers in each cells, list of mutated genes

    tot_population=1
    generation=1

    while tot_population >0 and tot_population<maxcell:
        tot_population=0

        G=[]
        for node in Tumor.traverse(strategy="postorder"): # traverse tree from the leaf, trick to add new child without disturbing the loop
            nb_driver=int(node.driver) # number of drivers per cell in the node
            nb_passenger=int(node.passenger) # number of passengers per cell in the node                            
            d=0.5*np.power((1.-fitness),nb_driver)/(np.power((1.-fitness_d),nb_passenger)) # probability to die  
            
            max_replicate=int(node.cell)
            nb_replicate=np.random.binomial(max_replicate,1-d,1) # number of cells that replicates with a success probabilitu 1-d given max_replicate can replicate

            if nb_replicate == 0: 
                node.cell=0
                G.append(node) #store the node for further deleting with chidren connected to the next possible parent

            else:
                nb_new_clones=np.random.binomial(nb_replicate,mut_rate,1) # among nb_replicate, nb_new_clones appear, each with a probability mut_rate
                node.cell=int(2*nb_replicate-nb_new_clones) # add to the current node, the number of cells that do not aquire a new mutations
                tot_population=tot_population+int(node.cell)
                mut_list=node.mutation # list of mutation

                for i in range(nb_new_clones):
                    num_gene=np.random.randint(1,maxgene) # gene that is mutated
                    nb_passenger=int(node.passenger) # number of passenger per cell in the node
                    new_clone=node.add_child() # add a new clone
                    new_clone.name=str(generation)
                    generation=generation+1
                    if num_gene <= maxdriver : # the gene is a driver
                        new_clone.add_features(cell=1,driver=nb_driver+1,passenger=nb_passenger,mutation=mut_list+"%d\t"%num_gene)
                    else: # the mutated gene is a passenger
                        new_clone.add_features(cell=1,driver=nb_driver,passenger=nb_passenger+1,mutation=mut_list+"%d\t"%num_gene)

        for node in G: # delete nodes with no replicate
            node.delete(prevent_nondicotomic=False)

    totcell=0
    for node in Tumor.traverse("postorder"):
        totcell=totcell+int(node.cell)

    if totcell > 1:
        nb_success+=1

        #save the whole tree
        tree_fn="%s/Tree_Clonal_evolution_%d.dat"%(OUTPUT_PATH,maxdriver,seed)
        for node in Tumor.traverse(strategy="postorder"):                                                                   
            node.write(features=[], outfile=tree_fn, format=2, is_leaf_fn=None, format_root_node=True)

        #save the size of each clone and its set of mutations
        clone_fn="%s/Clonal_evolution_%d.dat"%(OUTPUT_PATH,maxdriver,seed)
        clone_f=open(clone_fn,"w")
        clone_f.write("#mut_rate\t%f\tfitness\t%f\tfitness_d\t%f\n"%(mut_rate,fitness,fitness_d)) #write evolutionnary parameters on the first line
        for node in Tumor.traverse("postorder"):
            clone_f.write("%d\t%s\n"%(int(node.cell),node.mutation))
        clone_f.close()  
