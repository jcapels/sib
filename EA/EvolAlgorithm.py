# -*- coding: utf-8 -*-

from Popul import Popul
import numpy as np
from collections import defaultdict

class EvolAlgorithm:
    
    def __init__ (self, popsize, numits, noffspring, indsize):
        self.popsize = popsize
        self.numits = numits
        self.noffspring = noffspring
        self.indsize = indsize
        
    def initPopul(self, indsize):
        self.popul = Popul(self.popsize, indsize)
        
    def iteration(self):
        parents = self.popul.selection(self.noffspring)
        offspring = self.popul.recombination(parents, self.noffspring)
        self.evaluate(offspring)
        self.popul.reinsertion(offspring)
    
    def evaluate(self, indivs):
        for i in range(len(indivs)):
            ind = indivs[i]
            fit = 0.0
            for x in ind.getGenes():
                if x == 1: fit += 1.0
            ind.setFitness(fit)
        return None
    
    def run(self):
        self.initPopul(self.indsize)
        self.evaluate(self.popul.indivs)
        self.bestsol = []
        self.bestfit = 0.0
        for i in range(self.numits+1):            
            self.iteration()
            bs, bf = self.popul.bestSolution()
            if bf > self.bestfit:
                self.bestfit = bf
                self.bestsol = bs
            print("Iteration:", i, " ", "Best: ", self.popul.bestFitness())
        #self.bestsol, self.bestfit = self.popul.bestSolution()

def parser(file):
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
    res= np.zeros((len(dic.keys()),(len(dic.keys()))))
    for i in dic.keys():
        for j in range(i+1,len(dic.keys())+1):
            d=dist(i,j,dic)
            res[i-1][j-1]=d
            res[j-1][i-1]= d
    return res

def dist(i,j,dic):
    return (dic[i][0]-dic[j][0])**2+(dic[i][1]-dic[j][1])**2

def generate_blocks(mat,perc):
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
    return list(merged)


# merge function to  merge all sublist having common elements.
def merge_common(lists):
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

def test():
    ea = EvolAlgorithm(100, 4000, 50, 1000)
    ea.run()

if __name__=="__main__":
    dic = parser("qa194.tsp")
    mat = distmat(dic)
    print(generate_blocks(mat,0.3))




