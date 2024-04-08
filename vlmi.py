import multiprocessing
import modules
import timeit
import numpy as np
import random as rd
import matplotlib.pyplot as plt
import newick
from copy import deepcopy
import math
from itertools import chain
from collections import deque
from itertools import product


"""
Newick-fa inicializálása a Hofmeister et. al 2019 cikk alapján.
"""
newick_tree=''' ((Törzs:28,((((13.1:29):13,13.2:41):28,13.3:70):41,13.5:80):60):1,((((14.2:35):6,14.3:41):36,14.4:40):32,14.5:72):150):1'''
# példa futtatás:
Node = newick.loads(newick_tree)[0]


final_tree_14 = np.zeros((1,1,15))
final_tree_13 = np.zeros((1,1,15))


def simulate(k,j):
    d, p, b = modules.randoms(n = k, moran_steps = j, num_of_sim = 1)
    print(len(d),len(p),len(b))
    print(k,j)
    randoms1 = deque([[x] for x in d])
    randoms2 = deque(p)
    randoms3 = deque(b)

    tree_13_res = np.zeros(15)
    tree_14_res = np.zeros(15)
    for _ in range(1):
        sim = modules.DFS(node = Node,
                n = j,
                last_node_name = [None],
                last_node_children_number = [1],
                brain_list = [], 
                starter_list = [[[]]*k,[], 1], 
                results = [], 
                end_of_branch = [], 
                moran_steps_in_a_year = j, 
                mutation_number_expected_value = 1, 
                rds = randoms1, 
                rds_prop = randoms2, 
                rds_bud = randoms3)
        sim.pop()
        tree_14_mut = [rd.choice(sim[0]), rd.choice(sim[1]), rd.choice(sim[2]), rd.choice(sim[3])]
        tree_13_mut = [rd.choice(sim[4]), rd.choice(sim[5]), rd.choice(sim[6]), rd.choice(sim[7])]
        tree_14_res += np.array(modules.fun(tree_14_mut))
        tree_13_res += np.array(modules.fun(tree_13_mut))

    avg_tree_14 = tree_14_res/1
    avg_tree_13 = tree_13_res/1
    
    w = j-2
    u = j-1

    final_tree_14[u][int(w)] = avg_tree_14
    final_tree_13[u][int(w)] = avg_tree_13

    return final_tree_14, final_tree_13


def run_parallel(data):
    pool = multiprocessing.Pool()
    results = pool.starmap(simulate, data)
    pool.close()
    pool.join()
    print(len(results))


def run_sync(data):
    results = list(map(lambda l: simulate(*l), data))
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

def say_hello(): print("hello")
