import modules
from collections import deque
import newick
import random as rd
import numpy as np
"""
Newick-fa inicializálása a Hofmeister et. al 2019 cikk alapján.
"""
newick_tree=''' ((Törzs:28,((((13.1:29):13,13.2:41):28,13.3:70):41,13.5:80):60):1,((((14.2:35):6,14.3:41):36,14.4:40):32,14.5:72):150):1'''

d, p, b = modules.randoms(n = 7, moran_steps = 1, num_of_sim = 1000)

randoms1 = deque([[x] for x in d])
randoms2 = deque(p)
randoms3 = deque(b)

# példa futtatás:
tree_13_res = np.zeros(15)
tree_14_res = np.zeros(15)

for i in range(100):
    sim = modules.DFS(node = newick.loads(newick_tree)[0],
                    n = 7,
                    last_node_name = [None], 
                    last_node_children_number = [1], 
                    brain_list = [], starter_list = [[[]]*7,[], 1], 
                    results = [], 
                    end_of_branch = [], 
                    moran_steps_in_a_year = 1, 
                    mutation_number_expected_value = 1, 
                    rds = randoms1, 
                    rds_prop = randoms2, 
                    rds_bud = randoms3)
    sim.pop()
    tree_14_mut = [rd.choice(sim[0]), rd.choice(sim[1]), rd.choice(sim[2]), rd.choice(sim[3])]
    tree_13_mut = [rd.choice(sim[4]), rd.choice(sim[5]), rd.choice(sim[6]), rd.choice(sim[7])]
    tree_14_res += np.array(modules.fun(tree_14_mut))
    tree_13_res += np.array(modules.fun(tree_13_mut))

avg_tree_14 = tree_14_res/100
avg_tree_13 = tree_13_res/100

print(avg_tree_14, avg_tree_13)