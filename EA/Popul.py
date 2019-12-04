# -*- coding: utf-8 -*-

from Indiv import Indiv
from random import random
import numpy as np

class Popul:

    def __init__(self, popsize, indsize, maxValue, indivs = []):
        self.popsize = popsize
        self.indsize = indsize
        self.maxValue = maxValue
        if indivs: self.indivs = indivs
        else: self.initRandomPop()

    def getIndiv(self, index):
        return self.indivs[index]

    def initRandomPop(self):
        self.indivs = []
        for i in range(self.popsize):
            indiv_i = Indiv(self.indsize, self.maxValue, [])
            self.indivs.append(indiv_i)

    def getFitnesses(self):
        fitnesses = []
        for ind in self.indivs:
            fitnesses.append(ind.getFitness())
        return fitnesses
        
    def bestFitness(self):
        return max(self.getFitnesses())

    def bestSolution(self):
        fitnesses = self.getFitnesses()
        bestf = fitnesses[0]
        bestsol = 0
        for i in range(1,len(fitnesses)):
            if fitnesses[i] > bestf:
                bestf = fitnesses[i]
                bestsol = i
        return self.getIndiv(bestsol), bestf
    
    def selection(self, n):
        res = []
        fitnesses = list(self.linscaling(self.getFitnesses()))
        #fitnesses = list(self.getFitnesses())
        for i in range(n):
            sel = self.roulette(fitnesses)
            fitnesses[sel] = 0.0
            res.append(sel)
        return res
    
    def roulette(self, f):
        tot = sum(f)
        val = random()
        acum = 0.0
        ind = 0
        while acum < val:
            acum += (f[ind] / tot)
            ind += 1
        return ind-1
    
    def linscaling(self, fitnesses):
        mx = max(fitnesses)
        mn = min(fitnesses)
        res = []
        for f in fitnesses:
            val = (f-mn)/(mx-mn)
            res.append(val)
        return res
    
    def recombination(self, parents, noffspring):  #numero de descendentes igual ao número de pais
        offspring = []
        new_inds = 0
        while new_inds < noffspring:
            parent1 = self.indivs[parents[new_inds]]
            parent2 = self.indivs[parents[new_inds+1]]
            offsp1, offsp2 = parent1.crossover(parent2)
            offsp1.mutation()
            offsp2.mutation()
            offspring.append(offsp1)
            offspring.append(offsp2)
            new_inds += 2
        return offspring
       
    def reinsertion(self, offspring):   #formação de uma população com descendentes da recombinação e da seleção
        tokeep = self.selection(self.popsize-len(offspring))
        ind_offsp = 0
        for i in range(self.popsize):
            if i not in tokeep:
                self.indivs[i] = offspring[ind_offsp]
                ind_offsp += 1
        

    

















