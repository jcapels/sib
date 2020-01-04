# -*- coding: utf-8 -*-
from copy import copy

from Indiv import Indiv
from random import sample, random, randint, uniform
import numpy as np

class Popul:

    def __init__(self, popsize, indivs = []):
        self.popsize = popsize
        self.indivs = indivs
        self.ranking = self.getFitnessesAndIndivsSel()

    def getIndiv(self, index):
        return self.indivs[index]

    def getFitnesses(self):
        fitnesses = []
        for ind in self.indivs:
            fitnesses.append(ind.getFitness())
        return fitnesses

    def getFitnessesAndIndivsSel(self):
        """
        returns the individual and the respective fitness and further sorts by best fitness
        """
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
        """
        returns the best fitness
        """
        return min(self.getFitnesses())

    def bestSolution(self):
        """
        returns a tuple with the gene sequence corresponding to the best fitness and the respective fitness
        """
        fitnesses = self.getFitnesses()
        bestf = fitnesses[0]
        bestsol = 0
        for i in range(1,len(fitnesses)):
            if fitnesses[i] < bestf:
                bestf = fitnesses[i]
                bestsol = i
        return self.getIndiv(bestsol), bestf
    
    def selection(self, n):
        """
        returns the n best individuals of a population
        param n: number of individuals to be selected in the ranking
        """
        res = []
        fitnessesAndIndvs=self.ranking[:n]
        for i in range(n):
            res.append(fitnessesAndIndvs[i][0])
        return res

    def getRanking(self):
        return self.ranking

    def updateRanking(self):
        self.ranking=self.getFitnessesAndIndivsSel()
    
    def recombination(self, parents, noffspring):  #numero de descendentes igual ao número de pais
        """
        recombination of the population, starting with the selection of the parents and further random
        mutations in the offspring
        param parents: list of individuals for parents
        param noffspring: list of existing offspring
        """
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

    def random_mutations(self,perc=0.25):
        """
        random mutations in the worst percentage of the population. It performs mutations type 2 or 3
        a random number of times
        param perc: percentage of the population chosen to perform random mutations
        """
        inds = []
        n = int(len(self.indivs)*perc)
        for ind in self.ranking[len(self.indivs)-n:]:
            inds.append(ind[0])
        n = randint(1,10)
        for ind in inds:
            for i in range(n):
                m = randint(2, 3)
                self.getIndiv(ind).mutation(m)
        self.ranking = self.getFitnessesAndIndivsSel()


    def reinsertion(self, offspring):
        """
        insertion of the offspring into the population
        param offspring: list with all the offspring of a generation
        """
        tokeep = self.selection(self.popsize-len(offspring))
        ind_offsp = 0
        for i in range(self.popsize):
            if i not in tokeep:
                self.indivs[i] = offspring[ind_offsp]
                ind_offsp += 1
        self.ranking=self.getFitnessesAndIndivsSel()
        
    def migration(self,popul):
        """
        perfoms the migration of the best 15 individuals between islands
        param popul: the other island
        return:
        """
        best = popul.selection(15)
        best_indvs = []
        for i in best:
            best_indvs.append(copy(popul.getIndiv(i)))
        tokeep = self.selection(self.popsize - 15)
        ind_offsp = 0
        for i in range(self.popsize):
            if i not in tokeep:
                self.indivs[i] = best_indvs[ind_offsp]
                ind_offsp += 1
        self.updateRanking()

















