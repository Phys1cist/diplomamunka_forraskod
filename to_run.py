import module
import module
from collections import deque
import newick
import random as rd
import numpy as np
import math
import time


"""
Newick-fa inicializálása a Hofmeister et. al 2019 cikk alapján.
"""
newick_tree=''' ((Törzs:28,((((13.1:29):13,13.2:41):28,13.3:70):41,13.5:80):60):1,((((14.2:35):6,14.3:41):36,14.4:40):32,14.5:72):150):1'''
# példa futtatás:

Node = newick.loads(newick_tree)[0]

final_tree_14 = np.zeros((4,7,15))
final_tree_13 = np.zeros((4,7,15))

M = [1,2,3,4]
N = [2,4,6,8,16,32,64]


def run(M,N):

    for k in M:
        for j in N:
            d, p, b = module.randoms(n = j, moran_steps = k, num_of_sim = 800)
            randoms1 = deque([[x] for x in d])
            randoms2 = deque(p)
            randoms3 = deque(b)

            tree_13_res = np.zeros(15)
            tree_14_res = np.zeros(15)
            for i in range(800):
                sim = module.DFS(node = Node,
                        n = j,
                        last_node_name = [None],
                        last_node_children_number = [1],
                        brain_list = [], 
                        starter_list = [[[]]*j,[], 1], 
                        results = [], 
                        end_of_branch = [], 
                        moran_steps_in_a_year = k, 
                        mutation_number_expected_value = 1, 
                        rds = randoms1, 
                        rds_prop = randoms2, 
                        rds_bud = randoms3)
                sim.pop()
                tree_14_mut = [rd.choice(sim[0]), rd.choice(sim[1]), rd.choice(sim[2]), rd.choice(sim[3])]
                tree_13_mut = [rd.choice(sim[4]), rd.choice(sim[5]), rd.choice(sim[6]), rd.choice(sim[7])]
                tree_14_res += np.array(module.fun(tree_14_mut))
                tree_13_res += np.array(module.fun(tree_13_mut))

            avg_tree_14 = tree_14_res/800
            avg_tree_13 = tree_13_res/800
            
            w = int(math.log(j,2))-1
            u = k-1

            final_tree_14[u][int(w)] = avg_tree_14
            final_tree_13[u][int(w)] = avg_tree_13
    return final_tree_14, final_tree_13

#start_time = time.time()

res14, res13 = run(M,N)
with open("tree_14_res.txt", "w") as file:
    # Write all elements of the list to the file
    file.write("\n".join(map(str, res14)))
with open("tree_13_res.txt", "w") as file:
    # Write all elements of the list to the file
    file.write("\n".join(map(str, res13)))
#elapsed_time = time.time() - start_time
#print(elapsed_time)


