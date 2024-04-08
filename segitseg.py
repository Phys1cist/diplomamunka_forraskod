import multiprocessing
from functools import partial
import modules
import newick
import numpy as np
import random as rd
from collections import deque


"""
Newick-fa inicializálása a Hofmeister et. al 2019 cikk alapján.
"""
newick_tree=''' ((Törzs:28,((((13.1:29):13,13.2:41):28,13.3:70):41,13.5:80):60):1,((((14.2:35):6,14.3:41):36,14.4:40):32,14.5:72):150):1'''
Node = newick.loads(newick_tree)[0]

def simulate(n,
            moran_steps_in_a_year,   
            starter_list, 
            node,
            mutation_number_expected_value, 
            rds, 
            rds_prop, 
            rds_bud, 
            last_node_name, 
            last_node_children_number, 
            brain_list, 
            results, 
            end_of_branch,
            number_of_sims):
    tree_13_res = np.zeros(15)
    tree_14_res = np.zeros(15)

    for i in range(1):
        sim = modules.DFS(n,
                          moran_steps_in_a_year, 
                          starter_list, 
                          node, 
                          mutation_number_expected_value, 
                          rds, 
                          rds_prop, 
                          rds_bud, 
                          last_node_name, 
                          last_node_children_number, 
                          brain_list, 
                          results, 
                          end_of_branch)
        sim.pop()
        tree_14_mut = [rd.choice(sim[0]), rd.choice(sim[1]), rd.choice(sim[2]), rd.choice(sim[3])]
        tree_13_mut = [rd.choice(sim[4]), rd.choice(sim[5]), rd.choice(sim[6]), rd.choice(sim[7])]
        tree_14_res += np.array(modules.fun(tree_14_mut))
        tree_13_res += np.array(modules.fun(tree_13_mut))

    avg_tree_14 = tree_14_res/number_of_sims
    avg_tree_13 = tree_13_res/number_of_sims

    return avg_tree_14, avg_tree_13

"""
Newick-fa inicializálása a Hofmeister et. al 2019 cikk alapján.
"""
newick_tree=''' ((Törzs:28,((((13.1:29):13,13.2:41):28,13.3:70):41,13.5:80):60):1,((((14.2:35):6,14.3:41):36,14.4:40):32,14.5:72):150):1'''
Node = newick.loads(newick_tree)[0]

d, p, b = modules.randoms(n = 7, moran_steps = 1, num_of_sim = 10)

randoms1 = deque([[x] for x in d])
randoms2 = deque(p)
randoms3 = deque(b)




# Define a function to run simulations in parallel
def run_simulations():
    # Define the range of values for the 6 parameters you want to vary
    n = [7]
    moran_steps_in_a_year = [1]


    # Create a list of all parameter combinations
    all_param_combinations = [(p1, p2) for p1 in n
                              for p2 in moran_steps_in_a_year]

    # Create a partial function with fixed parameters
    partial_simulate = partial(simulate,
                               starter_list = [[[]]*7,[], 1],
                               node = Node,
                                mutation_number_expected_value = 1,
                               rds = randoms1,
                                rds_prop = randoms2, 
                                rds_bud = randoms3,
                                last_node_name = [None], 
                                last_node_children_number = [1], 
                                brain_list = [], 
                                results = [], 
                                end_of_branch = [],
                                number_of_sims = 10)

    # Create a multiprocessing pool
    num_processes = multiprocessing.cpu_count()  # Use all available CPU cores
    pool = multiprocessing.Pool(processes=num_processes)

    # Run simulations in parallel
    pool.starmap(partial_simulate, all_param_combinations)

    # Close the pool
    pool.close()
    pool.join()




if __name__ == "__main__":
    run_simulations()