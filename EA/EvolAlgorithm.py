# -*- coding: utf-8 -*-
from copy import copy, deepcopy
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

    def generate_indvs(self,block):
        """
        generate the individuals for the first generation using the previously built blocks
        """
        res = []
        for i in range(self.popsize):
            x = block.copy()
            shuffle(x)
            indv = Indiv(self.matdist, list(itertools.chain.from_iterable(x)))
            res.append(indv)
        return res
        
    def initPopul(self):
        """
        create the first generation of individuals in each island
        """
        self.populs=[]
        for block in self.blocks:
            indivs = self.generate_indvs(block)
            self.populs.append(Popul(self.popsize,indivs))
        
    def iteration(self,mode):
        """
        iterative process depending on the mode for each generation and further selection, recombination
        and reinsertion
        param mode: mode 1 will perform a random mutation in populations being the default one;
        mode 2 will be activated when the best fitness of the populations remains the same for the last 100 iterations.
        This will perform a random mutations in the 5th to the 9th individual with the best fitness and
        random mutations in the worst 70% individuals of the population, also in the 2nd to the 14th individual with
        the best fitness the mutation type 4 is performed;
         mode 3 is activated when the best fitness remains the
        same for the previous 50 iterations and performs random mutations in the worst 50% of the population;
        mode 4 is activated when the best fitness remains the same for the previous 150 iterations. This mode will
        perform a type 3 mutation in the individual with best fitness, random mutations in the worst 50% of the
        population and will promote migrations between islands
        """
        for i in range(len(self.populs)):
            if mode == 1:
                self.populs[i].random_mutations()
            elif mode==2:
                best_indexes2 = self.populs[i].getRanking()[4:10]
                for j in best_indexes2:
                    self.populs[i].getIndiv(j[0]).random_mutations_indiv()
                self.populs[i].updateRanking()
                self.populs[i].random_mutations(0.7)
                best_indexes2 = self.populs[i].getRanking()[1:15]
                for j in best_indexes2:
                    self.populs[i].getIndiv(j[0]).mutation(4)
                self.populs[i].updateRanking()

            elif mode==3:
                self.populs[i].random_mutations(0.5)

            elif mode==4 and i>0:
                n=randint(0,1)
                self.populs[i].random_mutations(0.5)
                self.populs[i].updateRanking()

                self.populs[i].migration(self.populs[i-1])
                self.populs[i-1].migration(self.populs[i])



            parents = self.populs[i].selection(self.noffspring)
            shuffle(parents)
            offspring = self.populs[i].recombination(parents, self.noffspring)
            self.populs[i].reinsertion(offspring)

    def run(self):
        """
        Will start the iteration process generation a population and retrieving the best fitness for each iteration.
        It will further select the modes for the iteration process and retrieve at the end the best sequence and
        respective fitness value
        """
        self.initPopul()
        self.bestsol = []
        self.bestfit = self.populs[0].getFitnesses()[0]
        l=1
        bs, bf = self.populs[0].getIndiv(self.populs[0].getRanking()[0][0]), self.populs[0].getIndiv(
            self.populs[0].getRanking()[0][0]).getFitness()

        previous=bf
        for i in range(self.numits+1):
            if l%50==0 and l%100!=0 and l%150!=0:
                self.iteration(3)
                print("---------Random mutations in the worst 50 % ------------")
            elif l%100==0:
                self.iteration(2)
                print("---------Random mutations in the worst 70 % and type 3 mutation in the best 3 to 9 + block swap------------")
            elif l%150==0:
                self.iteration(4)
                print("---------Random mutations in the worst 50 % and type 3 mutation in the best and migration between islands------------")
                l=1
            else:
                self.iteration(1)


            for popul in self.populs:
                newfit = popul.getIndiv(
                    popul.getRanking()[0][0]).getFitness()
                if bf > newfit:
                    bs, bf = popul.getIndiv(popul.getRanking()[0][0]), popul.getIndiv(popul.getRanking()[0][0]).getFitness()


            if bf!=previous:
                l=1
                previous=bf
            else:
                l+=1

            if bf < self.bestfit:
                self.bestfit = bf
                self.bestsol = deepcopy(bs)

            # if i % 100 == 0:
            #     print(bs.getGenes())
            #     print(self.bestsol.getGenes())
            print("Iteration:", i, " ", "Best: ", self.bestfit)
        print("Best solution: ", self.bestsol.getGenes())
        print("Best fitness: ", self.bestfit)

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
    blocks=[]
    for i in range(3):
        blocks.append(generate_blocks(mat, 0.86+(i/100)))
    ea = EvolAlgorithm(50, 10000, 26,blocks,mat)
    ea.run()



