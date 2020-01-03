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
        return int(self.distMat[i-1][j-1])

    def mutation(self, type=2):
        """
        performs a swap mutation of 2 or 3 genes, depending on the specified type, and further
        calculates the new fitness. It takes into account special cases where the swapped genes were neighbors
        param type: it can be either 2 or 3 depending on the number of swapped genes
        """
        if type==2:
            s = len(self.genes)
            pos1=0
            pos2=0
            #while pos1 == pos2 or pos1 - pos2 in [-1, 1, 193, -193]:
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
            pos2 = 0
            pos3 = 0
            #while pos1 == pos2 or pos2 == pos3 or pos1 == pos3 or pos1 - pos2 in [1, -1, 193, -193] or pos1 - pos3 in [1, -1, 193, -193] or pos2 - pos3 in [1, -1, 193, -193]:
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
        """
        performs a partially mapped crossover returning two new individuals
        param indiv2: genes from a second individual for recombination
        """
        child1 = self.crossover_pmx((self.genes.copy(),indiv2.getGenes().copy()))
        child2 = self.crossover_pmx((indiv2.getGenes().copy(),self.genes.copy()))
        return Indiv(self.distMat, child1), Indiv(self.distMat, child2)

    def crossover_pmx(self,parents):
        """
        performs a partially mapped crossover by choosing the middle fragment of the genes of one parent and completing
        the rest of the sequence with the genes from the other parent. For the second child the middle fragment of the
        other parent is used
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
        performs swap mutations (2 or 3 genes) n times
        """
        n = randint(1, 3)
        for i in range(n):
            m = randint(2, 3)
            self.mutation(m)

