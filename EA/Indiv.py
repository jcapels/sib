# -*- coding: utf-8 -*-
import unittest
from random import randint, random, shuffle
from math import sqrt
import numpy as np

class Indiv:
    
    def __init__(self, distMat, genes = []):
        self.distMat = distMat
        if genes:
            self.genes = genes
            self.calculate_fitness()


    def calculate_fitness(self):
        res=0
        i=0
        while i<len(self.genes)-1:
            res+=self.dist(self.genes[i],self.genes[i+1])
            i+=1
        res+=self.dist(self.genes[-1],self.genes[0])
        self.setFitness(res)

    def setFitness(self, fit):
        self.fitness = fit

    def getFitness(self):
        return self.fitness
    
    def getGenes(self):
        return self.genes

    def initRandom(self, size):
        self.genes = []
        for i in range(size):
            self.genes.append(randint(0, self.maxValue))

    def dist(self, i, j):
        return int(self.distMat[i-1][j-1])

    def mutation(self, type=2):
        if type==2:
            s = len(self.genes)
            pos1=0
            pos2=0
            # while pos1 == pos2 or pos1 - pos2 in [-1, 1, 193, -193]:
            while pos1==pos2:
                pos1 = randint(0, s - 1)
                pos2 = randint(0, s - 1)
            gene1 = self.genes[pos1]
            gene2 = self.genes[pos2]

            if pos1-pos2 in [-1,1,193,-193]:
                self.genes[pos1], self.genes[pos2] = gene2, gene1
                self.calculate_fitness()
            else:
                self.fitness-=(self.dist(self.genes[pos1-1],self.genes[pos1]) \
                          + self.dist(self.genes[pos1],self.genes[(pos1+1) % (s)]) \
                            + self.dist(self.genes[pos2],self.genes[(pos2+1) % (s)]) \
                          + self.dist(self.genes[pos2-1],self.genes[pos2]))

                self.fitness += self.dist(self.genes[pos1 - 1], self.genes[pos2]) \
                            + self.dist(self.genes[pos2], self.genes[(pos1+1) % (s)]) \
                            + self.dist(self.genes[pos1], self.genes[(pos2+1) % (s)]) \
                            + self.dist(self.genes[pos2 - 1],self.genes[pos1])
                self.genes[pos1], self.genes[pos2] = gene2, gene1

        elif type == 3:
            s = len(self.genes)
            pos1 = 0
            pos2 = 0
            pos3 = 0
            # while pos1 == pos2 or pos2 == pos3 or pos1 == pos3 or pos1 - pos2 in [1, -1, 193, -193] or pos1 - pos3 in [1, -1, 193, -193] or pos2 - pos3 in [1, -1, 193, -193]:
            while pos1 == pos2 or pos2 == pos3 or pos1==pos3:
                pos1 = randint(0, s - 1)
                pos2 = randint(0, s - 1)
                pos3 = randint(0, s - 1)

            if pos1-pos2 in [1,-1,193,-193] or pos1-pos3 in [1,-1,193,-193] or pos2-pos3 in [1,-1,193,-193]:
                self.genes[pos1], self.genes[pos2], self.genes[pos3] = self.genes[pos2], self.genes[pos3], self.genes[
                    pos1]
                self.calculate_fitness()

            else:
                value1 = (self.dist(self.genes[pos1 - 1], self.genes[pos1]) \
                           + self.dist(self.genes[pos1],self.genes[(pos1+1) % (s)]) \
                           + self.dist(self.genes[pos2], self.genes[(pos2+1) % (s)]) \
                            + self.dist(self.genes[pos2 - 1], self.genes[pos2]) \
                            + self.dist(self.genes[pos3 - 1], self.genes[pos3])  \
                            + self.dist(self.genes[pos3], self.genes[(pos3+1) % (s)]))

                self.fitness-=value1

                value2 = self.dist(self.genes[pos1 - 1], self.genes[pos2]) \
                            + self.dist(self.genes[pos2], self.genes[(pos1+1) % (s)]) \
                            + self.dist(self.genes[pos3 - 1], self.genes[pos1]) \
                            + self.dist(self.genes[pos1], self.genes[(pos3+1) % (s)]) \
                            + self.dist(self.genes[pos2 - 1], self.genes[pos3]) \
                            + self.dist(self.genes[pos3], self.genes[(pos2+1) % (s)])

                self.fitness += value2

                self.genes[pos1], self.genes[pos2], self.genes[pos3] = self.genes[pos2], self.genes[pos3], self.genes[
                    pos1]


    def crossover(self, indiv2):
        return self.one_pt_crossover(indiv2)

    def one_pt_crossover(self, indiv2):
        child1 = self.crossover_pmx((self.genes.copy(),indiv2.getGenes().copy()))
        child2 = self.crossover_pmx((indiv2.getGenes(),self.genes))
        return Indiv(self.distMat, child1), Indiv(self.distMat, child2)

    def crossover_pmx(self,parents):
        pos1 = 0
        pos2 = 0
        size = len(parents[0])
        while pos1 == pos2:
            pos1 = randint(0, size)
            pos2 = randint(0, size)
        if pos1>pos2:
            pos1,pos2 = pos2,pos1
        fragment1 = parents[0][pos1:pos2]
        fragment2 = parents[1][pos1:pos2]
        i_lst = []
        j_lst = []
        child = [0]*len(parents[0])
        child[pos1:pos2] = fragment1
        for i in range(len(fragment2)):
            if fragment2[i] not in fragment1:
                i_lst.append(fragment2[i])
                j_lst.append(fragment1[i])
        for i in range(len(i_lst)):
            pos = parents[1].index(j_lst[i])
            if child[pos] == 0:
                child[pos] = i_lst[i]
            else:
                k = child[pos]
                ind = parents[1].index(k)
                while child[ind] != 0:
                    k = child[ind]
                    ind = parents[1].index(k)
                child[ind] = i_lst[i]
        for c in range(len(child)):
            if child[c] == 0:
                child[c] = parents[1][c]
        return child

    def random_mutations_indiv(self):
        n = randint(1, 3)
        for i in range(n):
            m = randint(2, 3)
            self.mutation(m)


def distmat(dic):
    res= np.zeros((len(dic.keys()),(len(dic.keys()))))
    for i in dic.keys():
        for j in range(i+1,len(dic.keys())+1):
            d=dist(i,j,dic)
            res[i-1][j-1]=d
            res[j-1][i-1]= d
    return res

# def dist(i,j,dic):
#     return sqrt((dic[i][0]-dic[j][0])**2+(dic[i][1]-dic[j][1])**2)
#
# def parser(file):
#     with open(file) as f:
#         lines=f.readlines()
#         res={}
#         go=False
#         for line in lines:
#             if line.strip("\n") == "EOF":
#                 go=False
#             if go:
#                 line=line.strip("\n")
#                 l=line.split(' ')
#                 res[int(l[0])]=(float(l[1]), float(l[2]))
#             if line.strip("\n")=="NODE_COORD_SECTION":
#                 go=True
#         return res
# if __name__=="__main__":
#     genes = [122, 78, 96, 84, 83, 45, 12, 10, 9, 22, 163, 21, 17, 11, 14, 171, 148, 126, 125, 140, 169, 179, 173, 174, 138, 38, 109, 110, 6, 8, 156, 73, 155, 76, 129, 112, 69, 7, 4, 75, 1, 3, 5, 15, 19, 31, 41, 27, 39, 47, 65, 43, 53, 120, 55, 116, 42, 30, 32, 35, 44, 46, 40, 34, 117, 48, 52, 56, 36, 61, 119, 67, 82, 191, 37, 68, 79, 81, 111, 20, 51, 85, 103, 98, 130, 134, 137, 139, 181, 188, 193, 189, 16, 192, 194, 182, 186, 187, 190, 185, 133, 115, 105, 91, 72, 2, 102, 92, 97, 152, 184, 128, 124, 123, 29, 28, 141, 24, 26, 74, 106, 118, 131, 147, 175, 183, 66, 54, 165, 60, 153, 157, 154, 62, 150, 158, 136, 143, 135, 49, 107, 108, 100, 64, 57, 177, 70, 13, 88, 93, 95, 18, 121, 160, 166, 87, 162, 114, 33, 101, 58, 127, 149, 172, 176, 164, 50, 161, 86, 146, 145, 142, 113, 80, 71, 151, 77, 25, 23, 99, 94, 89, 59, 144, 63, 90, 132, 104, 167, 170, 178, 180, 168, 159]
#
#     dic = parser("qa194.tsp")
#     mat = distmat(dic)
#     new = Indiv(mat,genes)
#     print(new.getFitness())