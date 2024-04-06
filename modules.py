import numpy as np
import random as rd
import matplotlib.pyplot as plt
import newick
from copy import deepcopy
import math
from itertools import chain
from collections import deque



#ezeket választjuk ki osztódásra: bemenetek: l-kezdeti list, k-hányat választunk ki
def choose_from_population(l, k):
    """
    Kiválaszt a populációból k sejtet
    params:
    l: list, bemeneti populáció
    k: int, hányat sejtet válasszunk ki osztódásra
    """
    chosen_samples=rd.sample(range(len(l)),k)
    return chosen_samples

#sejtek elpusztítása
def kill(duplicated_list,for_death):
    """
    A megadott indexű elemeket kitörli a populációból
    params:
    duplicated_list: list, azon populáció, amelyben már osztódtak a sejtek
    for_death: list, azokat az indexeket tartalmazza, amiket ki kell hagyni
    """

    result2=[]
    for i in range(len(duplicated_list)):
        if i not in for_death:    
            result2.append(duplicated_list[i])
    return result2

def duplication_and_mutation_100(culture, for_duplication, m):
    """
    A populációból a kiválasztott sejteket osztatja és  biztosan ad mutációkat is ad az új sejtekhez
    params:
    culture: list, a populációt tartalmazza
    for_duplication: list, azon sejtek indexei, amelyeket osztódnak
    m: int, kezdeti mutáció száma. Érdemes 1-el indítani


    """
    
    after_duplication=[]
    for i in range(len(culture)):
        mutates=[]
        if i in for_duplication:
            mutates.append(deepcopy(culture[i]))
            mutates.append(deepcopy(culture[i]))  
            for j in range(2):
                mutates[j].append(m)   
                m = m+1
                after_duplication.append(deepcopy(mutates[j]))   
        else:
            after_duplication.append(deepcopy(culture[i]))
    return after_duplication,m

#sejtek osztódása és mutációk gyűjtése
def duplication_and_mutation(culture, for_duplication, m, mutation_number_expected_value):
    """
    A populációból a kiválasztott sejteket osztatja és mutációkat is ad az új sejtekhez
    params:
    culture: list, a populációt tartalmazza
    for_duplication: list, azon sejtek indexei, amelyeket osztódnak
    m: int, kezdeti mutáció száma. Érdemes 1-el indítani
    mutation_number_expected_value: int, a poisson eloszlás a várható értéke

    """

    after_duplication=[]
    for i in range(len(culture)):
        mutates=[]
        if i in for_duplication:
            mutates.append(deepcopy(culture[i]))
            mutates.append(deepcopy(culture[i]))  
            for j in range(2):
                mut_num = np.random.poisson(mutation_number_expected_value)
                if mut_num != 0:
                    for l in range(mut_num):
                        mutates[j].append(m)   
                        m = m+1
                    after_duplication.append(deepcopy(mutates[j]))
                else:
                    after_duplication.append(deepcopy(mutates[j]))    
        else:
            after_duplication.append(deepcopy(culture[i]))
    return after_duplication,m



def decider(n):
    """
    Ez a fv két számot ad vissza, hányszor kell duplázni a sejteket ágképzésnél és hány sejtet kell még 1-szer osztatni.
    params: 
    n: int
    """

    k = 0
    while 2**(k+1)<=n:
        k += 1
    r = n-pow(2,k)
    return k,r

def creating_buds_random(after_moran, n, rds_bud:deque):
    """
    Ez a fv random kiveszi az osztódás utáni populáció egyik tagját.
    params:
    after_moran: list, a Moran lépés utáni listát kapja meg
    n: int, populációméret
    """
    #print(after_moran)
    a = rds_bud.popleft()
    #print(a)
    leolt = kill(after_moran, [a]) 
    #print(leolt)
    bud = after_moran[a] #ide ugyanazt a számot beírni, mint egy sorral feljebb.
    return (leolt, bud)


#ez pedig az osztódás rész
def propagation(n,single_cell, starter_mutation, mutation_number_expected_value, rds_prop:deque): 
    """
    Az ágképzést valósítja meg. 1 sejtből csinál n sejtet.
    params:
    n: int, populáció mérete
    single_cell: list, 1 sejt, amit megkap az ágból
    starter_mutation: int, kezdeti mutációs szám
    mutation_number_expected_value: int, a poisson eloszlás a várható értéke
    """

    bemenet=deepcopy(single_cell)
    mutate_m=starter_mutation
    for i in range(decider(n)[0]):
        #print('Culture:', bemenet)
        #print(bemenet)
        valasztott=choose_from_population(bemenet,len(bemenet))
        #print(valasztott)
        kettozott_es_mutalt, mutate_m=duplication_and_mutation_100(bemenet,valasztott, mutate_m)
        #print(kettozott_es_mutalt)
        bemenet=deepcopy(kettozott_es_mutalt)
        #print('Culture:', bemenet)
    valasztott = rds_prop.popleft()
    
    #print(valasztott)
    kettozott_es_mutalt, mutate_m=duplication_and_mutation_100(bemenet, valasztott, mutate_m)
    #print(kettozott_es_mutalt)
    bemenet=deepcopy(kettozott_es_mutalt)
    #print(bemenet)
    last_mutation=mutate_m
    return(bemenet, last_mutation)


    
def simulate(to, n, starter_culture, starter_mutation, moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud): 
    """
    két elágazás közt végigmegy a szimuláció. Ha mellékág kezdődik, akkor először felsokszorozza a sejteket, majd megy tovább
    params:
    to: int, évek száma.
    n: int, sejtpopuláció mérete
    starter_culture: list, kezdeti populáció. 
    starter_mutation: int, kezdeti mutációs szám
    moran_steps_in_a_year: int, 1 év alatt hány lépést tegyen a szimuláció, azaz hány rügyet rakjon le. Összes lépés = to*moran_steps_in_a_year
    mutation_number_expected_value: int, a poisson eloszlás a várható értéke
    rds: list, tartalmazza a random számokat a fa létrehozására
    """
    
    mutate_m = starter_mutation
    culture = deepcopy(starter_culture)
    branching_buds = []
    #print('Kezdeti populáció:', culture)
    lim = to
    if len(culture) == 1:
        culture, mutate_m=propagation(n,culture,mutate_m, mutation_number_expected_value, rds_prop)
    for i in range(moran_steps_in_a_year * lim):
        #print('\n','Populáció a(z)', i+1,'. Moran lépés kezdetén', culture)
        k = 1

        choices = rds.popleft()
        #print('Sorszám osztódásra:',choices)
        duplicated_and_mutated, mutate_m = duplication_and_mutation_100(culture, choices, mutate_m)
        #print('Megkettőződtek és mutáltak',duplicated_and_mutated)

        #after_creating_buds = deepcopy(creating_buds_random(duplicated_and_mutated, n)[0])
        #print('Rügyképzés után:',after_creating_buds)  
        #if i == moran_steps_in_a_year * lim - 1:
        #    branching_buds=[deepcopy(creating_buds_random(duplicated_and_mutated, n)[1])]
        v1 , v2 = creating_buds_random(duplicated_and_mutated, n, rds_bud)
        after_creating_buds = deepcopy(v1)
        #print('Rügyképzés után:',after_creating_buds)  
        if i == moran_steps_in_a_year * lim - 1:
            branching_buds=[deepcopy(v2)]
        culture=deepcopy(after_creating_buds)
    #print('Kimeneti populáció', culture)
    #print('Rügy:', branching_buds)
    next_mutation = mutate_m
    return culture, branching_buds,next_mutation


def DFS(node, n, last_node_name, last_node_children_number,  brain_list, starter_list, results, end_of_branch, moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud):
    """
    Ez a fv az all in all, végigmegy az egész fán depth-first search algoritmussal, megcsinál minden Moran-lépést, ágképzést.
    params:
    node: a Newick-könyvtárral beolvasott fa-struktúra. 
    n: int, populációméret.
    last_node_name: list, az utolsó elágazási pont neve. Rekurzív a DFS fv, így [None]-t kell megadni kezdeti értéknek.
    last_node_children_number: int, az utolsó nód után hány ág van. Kezdeti érték legyen [1]
    brain_list: list, mindig a fő ág szimulációjának eredményét tartalmazza. Pl. végigmegy egy mellékágon, ami után már nincs több elágazás, visszaugrásnál
    kiveszi a lista utolsó elemét és azzal megy tovább.
    starter_list: list, kezdeti populáció. [[]]*n legyen
    results: list, a sim fv eredményeit adja vissza, [] kezdetben
    end_of_branch: list, az ág végén lévő populációt adja vissza (n sejt a mutációkkal)
    moran_steps_in_a_year: int, 1 év alatt hány lépést tegyen a szimuláció, azaz hány rügyet rakjon le. Összes lépés = to*moran_steps_in_a_year
    mutation_number_expected_value: int, a poisson eloszlás a várható értéke

    """
    lengths = []
    last_node_name.append(node.name)
    a = node.length
    children = node.descendants[::-1]
    last_node_children_number.append(len(children))

    if len(children) == 2 and last_node_children_number[-2] == 2:
        res = simulate(int(a), n , starter_list[1], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]
        results.append(res)
        brain_list.append(res[0])

    elif len(children) == 2 and last_node_children_number[-2] == 1:
        res = simulate(int(a), n , starter_list[0], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]
        results.append(res)
        brain_list.append(res[0])

    elif len(children) == 1 and last_node_children_number[-2] == 2:
        res = simulate(int(a), n , starter_list[1], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]

    elif len(children) == 1 and last_node_children_number[-2] == 1:
        res = simulate(int(a), n , starter_list[0], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]

    elif len(children) == 2 and last_node_children_number[-2] == 0:
        res = simulate(int(a), n , brain_list[-1], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]
        results.append(res)
        if len(brain_list) != 0:
            brain_list.pop() 
        brain_list.append(res[0])


    elif len(children) == 1 and last_node_children_number[-2] == 0: 
        res = simulate(int(a), n , brain_list[-1], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]
        if len(brain_list) != 0:
            brain_list.pop() 
    
    elif len(children) == 0 and last_node_children_number[-2] == 0:
        res = simulate(int(a), n , brain_list[-1], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]
        results.append(res)
        end_of_branch.append(res[0])
        if len(brain_list) != 0:
            brain_list.pop() 
    
    elif len(children) == 0 and last_node_children_number[-2] == 1:
        res = simulate(int(a), n , starter_list[0], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]
        results.append(res)
        end_of_branch.append(res[0])
    
    elif len(children) == 0 and last_node_children_number[-2] == 2: 
        res = simulate(int(a), n , starter_list[1], starter_list[2], moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)
        starter_list[0] = res[0]
        starter_list[1] = res[1]
        starter_list[2] = res[2]
        results.append(res)
        end_of_branch.append(res[0])
    
    lengths.append(a)
    for descendants in node.descendants[::-1]:

        kisebb_fa=DFS(descendants, n, last_node_name, last_node_children_number, brain_list, starter_list, results, end_of_branch, moran_steps_in_a_year, mutation_number_expected_value, rds, rds_prop, rds_bud)

    return end_of_branch

def flatten_chain(matrix):
    return list(chain.from_iterable(matrix))

def remove(my_list, my_set):
    for element in my_set:
        if element in my_list:
            my_list.remove(element)
    return my_list

def fun(nested_list):

    list1 = nested_list[0]
    list2 = nested_list[1]
    list3 = nested_list[2]
    list4 = nested_list[3]

    set1 = set(list1)
    set2 = set(list2)
    set3 = set(list3)
    set4 = set(list4)

    found_in_all = set1 & set2 & set3 & set4
    num_found_in_all = len(found_in_all)
    set1.difference_update(found_in_all)
    set2.difference_update(found_in_all)
    set3.difference_update(found_in_all)
    set4.difference_update(found_in_all)

    found_in_123 = set1 & set2 & set3
    num_found_in_123 = len(found_in_123)
    set1.difference_update(found_in_123)
    set2.difference_update(found_in_123)
    set3.difference_update(found_in_123)

    found_in_124 = set1 & set2 & set4
    num_found_in_124 = len(found_in_124)
    set1.difference_update(found_in_124)
    set2.difference_update(found_in_124)
    set4.difference_update(found_in_124)

    found_in_134 = set1 & set3 & set4
    num_found_in_134 = len(found_in_134)
    set1.difference_update(found_in_134)
    set3.difference_update(found_in_134)
    set4.difference_update(found_in_134)

    found_in_234 = set2 & set3 & set4
    num_found_in_234 = len(found_in_234)
    set2.difference_update(found_in_234)
    set3.difference_update(found_in_234)
    set4.difference_update(found_in_234)

    found_in_12 = set1 & set2
    num_found_in_12 = len(found_in_12)
    set1.difference_update(found_in_12)
    set2.difference_update(found_in_12)

    found_in_13 = set1 & set3
    num_found_in_13 = len(found_in_13)
    set1.difference_update(found_in_13)
    set3.difference_update(found_in_13)

    found_in_14 = set1 & set4
    num_found_in_14 = len(found_in_14)
    set1.difference_update(found_in_14)
    set4.difference_update(found_in_14)

    found_in_23 = set2 & set3
    num_found_in_23 = len(found_in_23)
    set2.difference_update(found_in_23)
    set3.difference_update(found_in_23)

    found_in_24 = set2 & set4
    num_found_in_24 = len(found_in_24)
    set2.difference_update(found_in_24)
    set4.difference_update(found_in_24)

    found_in_34 = set3 & set4
    num_found_in_34 = len(found_in_34)
    set3.difference_update(found_in_34)
    set4.difference_update(found_in_34)

    found_in_1 = set1
    num_found_in_1 = len(found_in_1)
    found_in_2 = set2
    num_found_in_2 = len(found_in_2)
    found_in_3 = set3
    num_found_in_3 = len(found_in_3)
    found_in_4 = set4
    num_found_in_4 = len(found_in_4)

    return num_found_in_1, num_found_in_2, num_found_in_3, num_found_in_4, num_found_in_12, num_found_in_13, num_found_in_14, num_found_in_23, num_found_in_24, num_found_in_34, num_found_in_123, num_found_in_124, num_found_in_134, num_found_in_234, num_found_in_all



def pairing(l_14, l_13):
    """
    A szimuláció eredményeit párosítja össze a tree_13-ra és a tree_14-re.
    params:
    l_14: tree_14-es fa szimulációjának eredménye
    l_13: tree_13-es fa szimulációjának eredménye
    """

    y = [ 0, #l_14[0], l_14[1], l_14[2], l_14[3],
            #l_13[0], l_13[1], l_13[2], l_13[3],
            l_14[4], l_14[5], l_14[6], l_14[7], l_14[8], l_14[9],
            l_13[4], l_13[5], l_13[6], l_13[7], l_13[8], l_13[9],
            l_14[10], l_14[11], l_14[12], l_14[13],
            l_13[10], l_13[11], l_13[12], l_13[13],
            #l_14[14],
            #l_13[14]
            ]
    return y

def plotting(tree_14, tree_13,):

    """
    A plottolást valósítja meg. bar plotot készít, amelyen szerepel a szimuláció eredménye és a mért adatok is.

    """
    fig, ax = plt.subplots(1, figsize=(12, 5))

    x_values = ['[1,2,3,4,5,6,7,8]',
                #'[1,0,0,0,0,0,0,0]', '[0,1,0,0,0,0,0,0]', '[0,0,1,0,0,0,0,0]', '[0,0,0,1,0,0,0,0]',
                #'[0,0,0,0,1,0,0,0]', '[0,0,0,0,0,1,0,0]', '[0,0,0,0,0,0,1,0]', '[0,0,0,0,0,0,0,1]',
                '[1,1,0,0,0,0,0,0]', '[1,0,1,0,0,0,0,0]', '[1,0,0,1,0,0,0,0]', '[0,1,1,0,0,0,0,0]', '[0,1,0,1,0,0,0,0]', '[0,0,1,1,0,0,0,0]',
                '[0,0,0,0,1,1,0,0]', '[0,0,0,0,1,0,1,0]', '[0,0,0,0,1,0,0,1]', '[0,0,0,0,0,1,1,0]', '[0,0,0,0,0,1,0,1]', '[0,0,0,0,0,0,1,1]',
                '[1,1,1,0,0,0,0,0]', '[1,1,0,1,0,0,0,0]', '[1,0,1,1,0,0,0,0]', '[0,1,1,1,0,0,0,0]',
                '[0,0,0,0,1,1,1,0]', '[0,0,0,0,1,1,0,1]', '[0,0,0,0,1,0,1,1]', '[0,0,0,0,0,1,1,1]',
                #'[1,1,1,1,0,0,0,0]',
                #'[0,0,0,0,1,1,1,1]'
                ]
    width = 0.45
    x_positions = np.arange(len(x_values))
    #y_values = [0, 12, 47, 30, 31, 25, 21, 34, 22, 0, 0, 1, 1, 1, 8, 0, 0, 0, 1, 0, 4, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0 ]
    y_values = [0, 0, 0, 1, 1, 1, 8, 0, 0, 0, 1, 0, 4, 0, 0, 0, 3, 0, 0, 0, 3 ]
    y_values_sim = pairing(tree_14, tree_13)
    #y_values_sim = pairing(np.array(tree_14), np.array(tree_13))

    for i in range(1,21):
        ax.text(x_positions[i] - width/2, y_values[i] + 0.3, f'{y_values[i]:.3f}', ha='center', va='bottom', fontsize=10, color='black', rotation=90)
        ax.text(x_positions[i] + width/2, y_values_sim[i] + 0.3, f'{y_values_sim[i]:.3f}', ha='center', va='bottom', fontsize=10, color='black', rotation=90)  

    bar1 = ax.bar(x_positions - width/2, y_values, width, color = 'green', label = 'Original data')
    bar2 = ax.bar(x_positions + width/2, y_values_sim, width, color = 'yellow', label = 'Simulation')
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_values, rotation = 90)
    for y_value in range(1, math.ceil(max(y_values))+1):
        ax.axhline(y = y_value, color = 'gray', linestyle = '--', linewidth=0.5)
    ax.set_title(f'Tuti kap mutációt, Merisztéma populáció: 7 sejt, Moran lépések száma: 1, random kivétel')
    ax.set_ylabel('SNPs')
    ax.legend()
    #plt.savefig('Random_kivétel_5_sejt', format = 'png')
    plt.show()


def loglikelihood(mu, f, k):
    L = np.log(mu)*np.sum(k) + np.sum(k*np.log(f)) - mu*np.sum(f)
    return L

def randoms(n, moran_steps, num_of_sim):
    prop = []
    v = int(np.power(2, np.floor(math.log(n, 2))))
    dupl = rd.choices(range(n), k = 804 * moran_steps * num_of_sim)
    for i in range(8 * num_of_sim):
        prop.append(rd.sample(range(v), n-v))
    bud = rd.choices(range(n), k = 804 * moran_steps * num_of_sim)

    return dupl, prop, bud
