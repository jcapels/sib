# -*- coding: utf-8 -*-

from random import randint, random, shuffle

class Indiv:
    
    def __init__(self, size, maxVal, genes = []):
        if genes:
            self.genes = genes
        else:
            self.initRandom(size)
        self.maxValue = maxVal
        self.fitness = None
        
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

    def mutation(self):
        #alterar
        s = len(self.genes)
        pos = randint(0, s - 1)
        self.genes[pos] = randint(0, self.maxValue)

    def crossover(self, indiv2):
        return self.one_pt_crossover(indiv2)

    def one_pt_crossover(self, indiv2):
        offsp1 = []
        offsp2 = []
        s = len(self.genes)
        pos = randint(0,s-1)
        for i in range(pos):
            offsp1.append(self.genes[i])
            offsp2.append(indiv2.genes[i])
        for i in range(pos, s):
            offsp2.append(self.genes[i])
            offsp1.append(indiv2.genes[i])
        return Indiv(s, offsp1), Indiv(s, offsp2)
    


