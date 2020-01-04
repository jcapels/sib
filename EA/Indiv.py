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
        """
        calculates the fitness of the individual through the distance between the genes
        """
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

    def dist(self, i, j):
        """
        retrieves the distance in the distance matrix corresponding to the cell ij
        param i: line number
        param j: column number
        """
        return self.distMat[i-1][j-1]

    def mutation(self, type=2):
        """
        performs a swap mutation depending on the specified type, and further calculates the new fitness.
        It takes into account special cases where the swapped genes were neighbors
        param type: it can be either 2 or 3 or 4. In type 2 it will do a swap of two random genes. Type 3 will swap two
        consecutive genes. Type 4 will move a block of genes
        """
        if type==2:
            s = len(self.genes)
            pos1=0
            pos2=0
            while pos1==pos2:
                pos1 = randint(0, s - 1)
                pos2 = randint(0, s - 1)
            gene1 = self.genes[pos1]
            gene2 = self.genes[pos2]

            if pos1-pos2 in [-1,193]:
                self.fitness -= self.dist(self.genes[pos1 - 1], self.genes[pos1]) \
                             + self.dist(self.genes[pos2], self.genes[(pos2 + 1) % (s)])

                self.fitness += self.dist(self.genes[pos1 - 1], self.genes[pos2]) \
                            + self.dist(self.genes[pos1], self.genes[(pos2 + 1) % (s)])

                self.genes[pos1], self.genes[pos2] = gene2, gene1

            elif pos1-pos2 in [1,-193]:
                self.fitness -= self.dist(self.genes[pos2 - 1], self.genes[pos2]) \
                                + self.dist(self.genes[pos1], self.genes[(pos1 + 1) % (s)])

                self.fitness += self.dist(self.genes[pos2 - 1], self.genes[pos1]) \
                                + self.dist(self.genes[pos2], self.genes[(pos1 + 1) % (s)])

                self.genes[pos1], self.genes[pos2] = gene2, gene1
                #self.calculate_fitness()
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
            # while pos1 == pos2 or pos1 - pos2 in [-1, 1, 193, -193]:
            pos1 = randint(0, s - 1)
            pos2 = (pos1 + 1) % s
            gene1 = self.genes[pos1]
            gene2 = self.genes[pos2]

            if pos1 - pos2 in [-1, 193]:
                self.fitness -= self.dist(self.genes[pos1 - 1], self.genes[pos1]) \
                                + self.dist(self.genes[pos2], self.genes[(pos2 + 1) % (s)])

                self.fitness += self.dist(self.genes[pos1 - 1], self.genes[pos2]) \
                                + self.dist(self.genes[pos1], self.genes[(pos2 + 1) % (s)])

                self.genes[pos1], self.genes[pos2] = gene2, gene1

            elif pos1 - pos2 in [1, -193]:
                self.fitness -= self.dist(self.genes[pos2 - 1], self.genes[pos2]) \
                                + self.dist(self.genes[pos1], self.genes[(pos1 + 1) % (s)])

                self.fitness += self.dist(self.genes[pos2 - 1], self.genes[pos1]) \
                                + self.dist(self.genes[pos2], self.genes[(pos1 + 1) % (s)])

                self.genes[pos1], self.genes[pos2] = gene2, gene1

        elif type==4:
            s = len(self.genes)
            pos1 = 0
            pos2 = 0
            while pos1 == pos2:
                pos1 = randint(0, s - 1)
                pos2 = randint(0, s - 1)

            if pos1>pos2:
                pos1,pos2=pos2,pos1
            if pos1 > pos2:
                pos1, pos2 = pos2, pos1
            fragment1 = self.genes[pos1:pos2]
            fragment2 = self.genes[:pos1] + self.genes[pos2:]
            pos3 = randint(0,s - 1)
            child = [0] * s
            for i in range(len(fragment1)):
                child[(pos3+i)%s]=fragment1[i]
            j=0
            for elem in fragment2:
                go=True
                while go and j<s:
                    if child[j]==0:
                        child[j] = elem
                        go=False
                    j+=1
            self.genes = child
            self.calculate_fitness()


    def crossover(self, indiv2):
        return self.one_pt_crossover(indiv2)

    def one_pt_crossover(self, indiv2):
        """
        performs a partially mapped crossover returning two new individuals
        param indiv2: genes from a second individual for recombination
        """
        child1 = self.crossover_pmx((self.genes.copy(),indiv2.getGenes().copy()))
        child2 = self.crossover_pmx((indiv2.getGenes().copy(),self.genes.copy()))
        return Indiv(self.distMat, child1), Indiv(self.distMat, child2)

    def crossover_pmx(self,parents):
        """
        performs a partially mapped crossover by choosing the middle fragment of the sequence of one parent and
        completing the rest of the sequence with the genes from the other parent.
        For the second child the middle fragment of the other parent is used
        param parents: list with both parents
        """
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
        """
        performs random mutations (random type) n times
        """
        n = randint(1, 3)
        for i in range(n):
            m = randint(2, 3)
            self.mutation(m)

def distmat(dic):
    res = np.zeros((len(dic.keys()), (len(dic.keys()))))
    for i in dic.keys():
        for j in range(i + 1, len(dic.keys()) + 1):
            d = dist(i, j, dic)
            res[i - 1][j - 1] = d
            res[j - 1][i - 1] = d
    return res


def dist(i, j, dic):
    return int(sqrt((dic[i][0] - dic[j][0]) ** 2 + (dic[i][1] - dic[j][1]) ** 2))


def parser(file):
    with open(file) as f:
        lines = f.readlines()
        res = {}
        go = False
        for line in lines:
            if line.strip("\n") == "EOF":
                go = False
            if go:
                line = line.strip("\n")
                l = line.split(' ')
                res[int(l[0])] = (float(l[1]), float(l[2]))
            if line.strip("\n") == "NODE_COORD_SECTION":
                go = True
        return res


if __name__ == "__main__":
    genes = [138, 139, 144, 141, 147, 152, 153, 150, 157, 154, 126, 125, 127, 132, 134, 137, 140, 145, 142, 146, 149, 156, 161, 163, 164, 169, 176, 172, 174, 173, 175, 177, 181, 186, 179, 182, 194, 190, 187, 183, 184, 189, 192, 191, 188, 122, 103, 17, 14, 13, 11, 7, 5, 9, 10, 12, 15, 19, 22, 18, 21, 24, 26, 28, 29, 27, 31, 32, 30, 35, 38, 34, 39, 47, 40, 43, 41, 44, 46, 48, 53, 52, 54, 42, 50, 49, 55, 56, 58, 51, 37, 45, 33, 60, 57, 64, 66, 61, 67, 73, 68, 70, 77, 84, 81, 79, 83, 88, 92, 95, 93, 96, 97, 100, 105, 106, 118, 107, 108, 110, 112, 115, 116, 117, 121, 120, 123, 124, 128, 129, 131, 136, 135, 133, 143, 148, 155, 151, 158, 159, 165, 162, 160, 166, 167, 168, 170, 171, 185, 193, 180, 178, 130, 98, 86, 85, 6, 4, 8, 16, 76, 75, 78, 74, 72, 25, 23, 71, 82, 80, 87, 102, 91, 69, 3, 2, 1, 20, 65, 63, 36, 59, 62, 89, 90, 94, 99, 101, 104, 111, 114, 109, 113, 119]


    dic = parser("qa194.tsp")
    mat = distmat(dic)
    new = Indiv(mat, genes)
    print(new.getFitness())