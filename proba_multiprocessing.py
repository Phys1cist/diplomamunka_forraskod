import multiprocessing
import timeit
import modules
import numpy as np
import newick
import random as rd

"""
Newick-fa inicializálása a Hofmeister et. al 2019 cikk alapján.
"""
newick_tree=''' ((Törzs:28,((((13.1:29):13,13.2:41):28,13.3:70):41,13.5:80):60):1,((((14.2:35):6,14.3:41):36,14.4:40):32,14.5:72):150):1'''
Node = newick.loads(newick_tree)[0]

def execute(n, 
            starter_list, 
            node, 
            moran_steps_in_a_year, 
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
                          starter_list, 
                          node, 
                          moran_steps_in_a_year, 
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


def run_parallel(data):
    pool = multiprocessing.Pool()
    results = pool.starmap(execute, data)
    pool.close()
    pool.join()
    print(len(results))


def run_sync(data):
    results = list(map(lambda l: execute(*l), data))
    print(len(results))


def main(is_parallel, data):
    if is_parallel:
        runtime = timeit.timeit(
            lambda:run_parallel(data),
            number=1
        )
    else:
        runtime = timeit.timeit(
            lambda:run_sync(data),
            number=1
        )
    print(runtime)


