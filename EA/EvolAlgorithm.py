# -*- coding: utf-8 -*-
from copy import copy
from math import sqrt

from Popul import Popul
from Indiv import Indiv
import numpy as np
from collections import defaultdict
from random import shuffle, randint, sample
import itertools



class EvolAlgorithm:
    
    def __init__ (self, popsize, numits, noffspring,blocks,mat):
        self.popsize = popsize
        self.numits = numits
        self.noffspring = noffspring
        self.blocks = blocks
        self.matdist = mat

    def generate_indvs(self):
        """
        generate the individuals for the first generation using the previously built blocks
        """
        res = []
        for i in range(self.popsize):
            x = self.blocks.copy()
            shuffle(x)
            indv = Indiv(self.matdist, list(itertools.chain.from_iterable(x)))
            res.append(indv)
        return res
        
    def initPopul(self):
        """
        create the first generation of individuals
        """
        indivs=self.generate_indvs()
        self.popul= Popul(self.popsize,indivs)
        
    def iteration(self,mode):
        """
        iterative process depending on the mode for each generation and further selection, recombination
        and reinsertion
        param mode: mode 1 will perform a random mutation in the population being the default one;
        mode 2 will be activated when the best fitness of the population remains the same for the last 100 iterations.
        This will perform a type 3 mutation in the 4th to the 10th individual with the best fitness and
        random mutations in the worst 70% of the population; mode 3 is activated when the best fitness remains the
        same for the previous 50 iterations and performs random mutations in the worst 50% of the population;
        mode 4 is activated when the best fitness remains the same for the previous 250 iterations. This mode will
        perform a type 2 mutation in the individual with best fitness and random mutations in the worst 50% of the
        population
        """
        if mode == 1:
            self.popul.random_mutations()
        elif mode==2:
            best_indexes2 = self.popul.getRanking()[4:10]
            for i in best_indexes2:
                self.popul.getIndiv(i[0]).mutation(3)
            self.popul.updateRanking()
            self.popul.random_mutations(0.7)

        elif mode==3:
            self.popul.random_mutations(0.5)

        elif mode==4:
            self.popul.random_mutations(0.5)
            self.popul.getIndiv(self.popul.getRanking()[0][0]).mutation(3)
            self.popul.updateRanking()


        parents = self.popul.selection(self.noffspring)
        offspring = self.popul.recombination(parents, self.noffspring)
        self.popul.reinsertion(offspring)

    def run(self):
        """
        Will start the iteration process generation a population and retrieving the best fitness for each iteration.
        It will further select the modes for the iteration process and retrieve at the end the best sequence and
        respective fitness value
        """
        self.initPopul()
        self.bestsol = []
        self.bestfit = self.popul.getFitnesses()[0]
        l=1
        previous=self.popul.bestFitness()
        for i in range(self.numits+1):
            if l%50==0 and l%100!=0 and l%250!=0:
                self.iteration(3)
                print("---------Random mutations in the worst 50 % ------------")
            elif l%100==0:
                self.iteration(2)
                print("---------Random mutations in the worst 70 % and 3-swap mutation in the best 4 to 10------------")
            elif l%250==0:
                self.iteration(4)
                print("---------Random mutations in the worst 50 % and 2-swap mutation in the best ------------")
                l=1
            else:
                self.iteration(1)

            bs, bf = self.popul.bestSolution()

            if bf!=previous:
                l=1
                previous=bf
            else:
                l+=1

            if bf < self.bestfit:
                self.bestfit = bf
                self.bestsol = bs
            print("Iteration:", i, " ", "Best: ", self.popul.bestFitness())
        print("Best solution: ", self.bestsol.getGenes())
        print("Best fitness: ", self.bestfit)
        #self.bestsol, self.bestfit = self.popul.bestSolution()

def parser(file):
    """
    creates a dictionary with the node id as keys and the corresponding coordinates as values
    param file: text file with the node id and respective coordinates
    """
    with open(file) as f:
        lines=f.readlines()
        res={}
        go=False
        for line in lines:
            if line.strip("\n") == "EOF":
                go=False
            if go:
                line=line.strip("\n")
                l=line.split(' ')
                res[int(l[0])]=(float(l[1]), float(l[2]))
            if line.strip("\n")=="NODE_COORD_SECTION":
                go=True
        return res

def distmat(dic):
    """
    creates a distance matrix calculating the distance between each node
    param dic: dictionary with the node id as keys and respective coordinates as values
    """
    res= np.zeros((len(dic.keys()),(len(dic.keys()))))
    for i in dic.keys():
        for j in range(i+1,len(dic.keys())+1):
            d=dist(i,j,dic)
            res[i-1][j-1]=d
            res[j-1][i-1]= d
    return res

def dist(i,j,dic):
    """
    calculates the euclidean distance between each node
    param i: line of the matrix
    param j: column of the matrix
    param dic: dictionary with the node id as keys and respective coordinates as values
    """
    return int(sqrt((dic[i][0]-dic[j][0])**2+(dic[i][1]-dic[j][1])**2))

def generate_blocks(mat,perc):
    """
    generates blocks of close nodes further merging the blocks with common nodes between them. Returns a list of blocks
    param mat: distance matrix
    param perc: percentage of blocks to be kept for individual generation
    """
    res = []
    for i in range(len(mat)):
        min=None
        for j in range(i+1,len(mat)):
            if min==None:
                min = (i+1,j+1,mat[i][j])
            elif min[2] > mat[i][j]:
                min = (i+1,j+1,mat[i][j])
        if min is not None:
            res.append(min)
    n = len(res)
    for i in range(n):
        for j in range(0, n - i - 1):
            if res[j][2] > res[j + 1][2]:
                res[j], res[j + 1] = res[j + 1], res[j]
    n_blocks = int(len(res)*perc)
    filtered = res[:n_blocks]
    blocks = [ [r[0],r[1]] for r in filtered ]
    merged = merge_common(blocks)
    merged = list(merged)
    merged_blocks = list(itertools.chain.from_iterable(merged))
    for g in [i for i in range(1,len(mat[0])+1)]:
        if g not in merged_blocks:
            merged.append([g])
    return merged

def merge_common(lists):
    """
    merge function to  merge all sublist having common elements
    param lists: lists of blocks
    """
    neigh = defaultdict(set)
    visited = set()
    for each in lists:
        for item in each:
            neigh[item].update(each)

    def comp(node, neigh=neigh, visited=visited, vis=visited.add):
        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            vis(node)
            nodes |= neigh[node] - visited
            yield node

    for node in neigh:
        if node not in visited:
            yield sorted(comp(node))





if __name__=="__main__":
    dic = parser("qa194.tsp")
    mat = distmat(dic)
    blocks=generate_blocks(mat,0.86)
    ea = EvolAlgorithm(120, 1000, 50,blocks,mat)
    ea.run()



