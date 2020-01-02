# -*- coding: utf-8 -*-

from Indiv import Indiv
from random import sample, random, randint, uniform
import numpy as np

class Popul:

    def __init__(self, popsize, indivs = []):
        self.popsize = popsize
        self.indivs = indivs
        self.ranking = self.getFitnessesAndIndivsSel()


    def addIndivs(self,indiv):
        self.indivs.extend(indiv)
        self.updateRanking()
        inds=[]
        for i in range(len(indiv)):
            ind=self.ranking[-i-1][0]
            inds.append(self.indivs[ind])
        for individuo in inds:
            self.indivs.remove(individuo)
        self.updateRanking()


    def getIndiv(self, index):
        return self.indivs[index]

    def getFitnesses(self):
        fitnesses = []
        for ind in self.indivs:
            fitnesses.append(ind.getFitness())
        return fitnesses

    def getFitnessesAndIndivsSel(self):
        fitnessesAndIndvs = []
        for ind in range(len(self.indivs)):
            fitnessesAndIndvs.append((ind,self.indivs[ind].getFitness()))
        size = len(fitnessesAndIndvs)
        for i in range(size):
            for j in range(0, size - i - 1):
                if fitnessesAndIndvs[j][1] > fitnessesAndIndvs[j + 1][1]:
                    fitnessesAndIndvs[j], fitnessesAndIndvs[j + 1] = fitnessesAndIndvs[j + 1], fitnessesAndIndvs[j]
        return fitnessesAndIndvs

    def bestFitness(self):
        return min(self.getFitnesses())

    def bestSolution(self):
        fitnesses = self.getFitnesses()
        bestf = fitnesses[0]
        bestsol = 0
        for i in range(1,len(fitnesses)):
            if fitnesses[i] < bestf:
                bestf = fitnesses[i]
                bestsol = i
        return self.getIndiv(bestsol), bestf
    
    def selection(self, n ,type="t"):
        res = []
        if type=="t":
            fitnessesAndIndvs=self.ranking[:n]
            for i in range(n):
                res.append(fitnessesAndIndvs[i][0])
        elif type=="roul":
            fitnesses = list(self.linscaling(self.getFitnesses()))
            #fitnesses = list(self.getFitnesses())
            for i in range(n):
                sel = self.roulette(fitnesses)
                fitnesses[sel] = 0.0
                res.append(sel)
        return res

    def getRanking(self):
        return self.ranking

    def updateRanking(self):
        self.ranking=self.getFitnessesAndIndivsSel()

    def roulette(self, f):
        tot = sum(f)
        val = uniform(0.8,0.9)
        acum = 1
        ind = 0
        while acum > val:
            acum -= (f[ind] / tot)
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
            offsp1.random_mutations_indiv()
            offsp2.random_mutations_indiv()
            offspring.append(offsp1)
            offspring.append(offsp2)
            new_inds += 2
        return offspring

    def random_mutations(self,perc=0.25,direction="last"):
        inds = []
        n = int(len(self.indivs)*perc)
        if direction=="last":
            for ind in self.ranking[len(self.indivs)-n:]:
                inds.append(ind[0])
        else:
            for ind in self.ranking[:n]:
                inds.append(ind[0])
        n = randint(1,10)
        for ind in inds:
            for i in range(n):
                m = randint(2, 3)
                self.getIndiv(ind).mutation(m)
        self.ranking = self.getFitnessesAndIndivsSel()


    def reinsertion(self, offspring):
        tokeep = self.selection(self.popsize-len(offspring))
        ind_offsp = 0
        for i in range(self.popsize):
            if i not in tokeep:
                self.indivs[i] = offspring[ind_offsp]
                ind_offsp += 1
        self.ranking=self.getFitnessesAndIndivsSel()
        

    

















