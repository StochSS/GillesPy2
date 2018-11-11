from gillespy2.core import GillesPySolver
import random
import math
import numpy as np
import sys


class TauLeapingSolver(GillesPySolver):
    name = "TauLeaping"
    """ TODO
    """
    # initializes the various variables needed for tau leaping algorithm
    @classmethod
    def initialize(self):
        self.curr_state['vol'] = self.model.volume
        self.results['time'] = []
        self.currentTime = 0
        self.outputTime = 0      
        self.nextTime = 0
        self.listOfCriticalSpecies = []
        self.listOfNonCriticalSpecies = []
        self.failedLeaps = 0

        for parameter in self.model.listOfParameters:
            self.curr_state[parameter] = self.model.listOfParameters[parameter].value        
        
        for species in self.model.listOfSpecies:
            # Initialize species populations
            self.curr_state[species] = self.model.listOfSpecies[species].initial_value    
            self.results[species]=[]
            self.listOfAffectedReactions[species] = [] 
            self.listOfNonCriticalSpecies.append(species) # default all species to non-critical
        
        for reaction in self.model.listOfReactions:
            # build a list of which reactions are affected by each species
            for reactant in self.model.listOfReactions[reaction].reactants:
                self.listOfAffectedReactions[str(reactant)].append(reaction)
            self.isCritical[reaction] = False # default every reaction to be critical 
            self.previousReactionCounts[reaction] = []


    # update which reactions are critical
    @classmethod
    def updateCriticalLists(self):
        # from critical list to non-critical list
        toRemove = []
        for criticalSpecies in self.listOfCriticalSpecies:
            if self.curr_state[criticalSpecies] > self.criticalThreshold:
                self.listOfNonCriticalSpecies.append(criticalSpecies)
                # modify the tag for every reaction affected by this species
                for affectedReaction in self.listOfAffectedReactions[criticalSpecies]:
                    changeTag = True
                    # check if other reactants for this reaction are critical
                    for affectedSpecies in self.model.listOfReactions[affectedReaction].reactants:
                        if affectedSpecies in self.listOfCriticalSpecies:
                            changeTag = False
                            break
                    if changeTag:
                        self.isCritical[affectedReaction] = False
                # add criticalSpecies to list to be removed
                toRemove.append(criticalSpecies)
        # remove former critical species
        for species in toRemove:
            self.listOfCriticalSpecies.remove(species)
        toRemove = []
        # from non-critical list to critical list
        for nonCriticalSpecies in self.listOfNonCriticalSpecies:
            if self.curr_state[nonCriticalSpecies] <= self.criticalThreshold:
                self.listOfCriticalSpecies.append(nonCriticalSpecies)
                # modify tag for each reaction affected by this species
                for affectedReaction in self.listOfAffectedReactions[nonCriticalSpecies]:
                    self.isCritical[affectedReaction] = True
                # add nonCriticalSpecies to list to be removed
                toRemove.append(nonCriticalSpecies)
        # remove former non-critical species
        for species in toRemove:
            self.listOfNonCriticalSpecies.remove(species)


    # selects tau based on the algorithm found in Cao et al 2006
    @classmethod
    def selectTau(self):
        # for critical reactions
        self.criticalPropensitySum = 0
        for reaction in self.model.listOfReactions:
            if self.isCritical[reaction]:
                self.criticalPropensitySum += self.propensities[reaction]
        if self.criticalPropensitySum != 0.0:
            self.criticalStepsize = -1*math.log(random.random())/self.criticalPropensitySum
        else:
            self.criticalStepsize = -1
       
        # for noncritical reactions
        # generate mu and sigmaSquared
        mu = {}
        sigmaSquared = {}
        for species in self.model.listOfSpecies:
            mu[species] = 0
            sigmaSquared[species] = 0
        for reaction in self.model.listOfReactions:
            if not self.isCritical[reaction]:
                for product in self.model.listOfReactions[reaction].products:
                    mu[str(product)] += self.model.listOfReactions[reaction].products.get(product, 0) * self.propensities[reaction]
                    sigmaSquared[str(product)] += self.model.listOfReactions[reaction].products.get(product, 0)
                for reactant in self.model.listOfReactions[reaction].reactants:
                    mu[str(reactant)] -= self.model.listOfReactions[reaction].reactants.get(reactant, 0) * self.propensities[reaction]
                    sigmaSquared[str(reactant)] -= self.model.listOfReactions[reaction].reactants.get(reactant, 0)
                for species in sigmaSquared:
                    sigmaSquared[species] = sigmaSquared[species] * sigmaSquared[species] * self.propensities[reaction]
        self.noncriticalStepsize = sys.maxsize
        g = 3.0 # should improve
        for species in self.model.listOfSpecies:
            numerator = self.epsilon * self.curr_state[species]/g
            numerator = max(numerator, 1.0)
            if mu[species] != 0:
                temp1 = numerator/abs(mu[species])
            else:
                temp1 = sys.maxsize
            if sigmaSquared[species] != 0:
                temp2 = numerator*numerator/sigmaSquared[species]
            else:
                temp2 = sys.maxsize
            self.noncriticalStepsize = min(self.noncriticalStepsize, temp1)
            self.noncriticalStepsize = min(self.noncriticalStepsize, temp2)

    # selects non-critical reaction to fire and handles poisson distribution
    @classmethod
    def selectReactions(self):
        self.criticalReaction = None
        if self.runCritical:
            # handle both critical and noncritical reactions
            # for critical reactions
            self.criticalReaction = None
            # generate uniform random number between 0 and criticalPropensitySum
            reactionNum = random.uniform(0, self.criticalPropensitySum)
            jsum = 0
            for reaction in self.model.listOfReactions:
                # for critical reactions
                if self.isCritical[reaction]:
                    self.previousReactionCounts[reaction] = 0
                    if jsum < reactionNum:
                       self.criticalReaction = reaction
                       jsum += self.propensities[reaction]
                # for noncritical reactions
                else:
                    self.previousReactionCounts[reaction] = np.random.poisson(self.propensities[reaction]*self.tau)
        else:
            # handle only noncritical reactions
            for reaction in self.model.listOfReactions:
                if self.isCritical[reaction]:
                    self.previousReactionCounts[reaction] = 0
                else:
                    self.previousReactionCounts[reaction] = np.random.poisson(self.propensities[reaction]*self.tau)
        #reactionsLastLeap = norm of self.previousReactionCounts -- for after selectTau is working

    # simmulates the reaction firings
    @classmethod
    def fireReactions(self):
        for species in self.model.listOfSpecies:
            self.populationChange[species] = 0
        # calculate the population changes
        for reaction in self.model.listOfReactions:
            for reactant in self.model.listOfReactions[reaction].reactants:
                self.populationChange[str(reactant)] -= self.model.listOfReactions[reaction].reactants[reactant]*self.previousReactionCounts[reaction]
            for product in self.model.listOfReactions[reaction].products:
                self.populationChange[str(product)] += self.model.listOfReactions[reaction].products[product]*self.previousReactionCounts[reaction]
        # append changes to curr_state
        for species in self.model.listOfSpecies:
            self.curr_state[str(species)] += self.populationChange[str(species)]

        # if we need to handle a critical reaction
        if self.criticalReaction is not None:
            for reactant in self.model.listOfReactions[self.criticalReaction].reactants:
                self.curr_state[str(reactant)] -= self.model.listOfReactions[self.criticalReaction].reactants[reactant]
            for product in self.model.listOfReactions[self.criticalReaction].products:
                self.curr_state[str(product)] += self.model.listOfReactions[self.criticalReaction].products[product]

        # check for negative populations
        self.negativePopulation = False
        self.leapFailed = False
        for species in self.model.listOfSpecies:
            if self.curr_state[str(species)] < 0:
                self.negativePopulation = True
                break
        if self.negativePopulation:
            # we need to undo what we just did
            for species in self.model.listOfSpecies:
                self.curr_state[str(species)] -= self.populationChange[str(species)]
            if self.criticalReaction is not None:
                # reverse self.criticalReaction firing. Should be criticalRollBack()
                for reactant in self.model.listOfReactions[self.criticalReaction].reactants:
                    self.curr_state[str(reactant)] += self.model.listOfReactions[self.criticalReaction].reactants[reactant]
                for product in self.model.listOfReactions[self.criticalReaction].products:
                    self.curr_state[str(product)] -= self.model.listOfReactions[self.criticalReaction].products[product]
            self.leapFailed = True
        else:
            self.leapFailed = False

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            debug=False, profile=False, show_labels=False, **kwargs):
        self.model = model
        self.simulation_data = []
        self.tau = 0
        self.curr_state = {}
        self.propensities = {}
        self.results = {}
        self.listOfAffectedReactions = {}
        self.isCritical = {}    # keyed by reaction
        self.criticalThreshold = 5  # threshold for when a species becomes critical
        self.epsilon = 0.03 # bound on relative change per time step
        self.previousReactionCounts = {}
        self.populationChange = {}
        self.criticalStepsize = 0
        self.nonCriticalStepsize = 0
        self.SSAThreshold = .005

        for traj_num in range(number_of_trajectories):
            # initialize everything
            self.initialize()

            # run the algorithm
            while(self.currentTime < t):
                while self.currentTime >= self.outputTime and self.outputTime < t:
                    self.results['time'].append(self.outputTime)
                    for species in self.model.listOfSpecies:
                        self.results[species].append(self.curr_state[species])
                    self.outputTime += increment

                # update the lsits of critical species
                self.updateCriticalLists()
                
                # update propensities
                for reaction in self.model.listOfReactions: 
                    self.propensities[reaction] = eval(self.model.listOfReactions[reaction].propensity_function, self.curr_state)
                
                # select the tau
                self.selectTau()
      
                # if we don't need to worry about critical reactions
                if self.noncriticalStepsize < self.criticalStepsize or self.criticalStepsize == -1:
                    self.tau = self.noncriticalStepsize
                    self.runCritical = False
                # otherwise, we do need to worry about critical reactions
                else:
                    self.tau = self.criticalStepsize
                    self.runCritical = True

                # avoid leaping past next output time
                if self.currentTime + self.tau > self.outputTime: # and currentTime <= t:
                    self.tau = self.outputTime - self.currentTime
                    self.nextTime = self.outputTime
                else:
                    self.nextTime = self.currentTime + self.tau

                # if we need to do SSA
                doSSA = False
                if doSSA:
                    # do SSA
                    1+1
                    #for species in model.listOfSpecies:
                        #model.listOfSpecies[species].initial_value = self.curr_state[str(species)]
                else:
                    # do tau leaping

                    self.selectReactions()
                    self.fireReactions()

                    if self.leapFailed:
                        self.failedLeaps += 1
                    else:
                        self.currentTime = self.nextTime
                        self.failedLeaps = 0

                    # record output
                    while self.currentTime >= self.outputTime:
                        self.results['time'].append(self.outputTime)
                        for species in self.model.listOfSpecies:
                            self.results[species].append(self.curr_state[species])
                        self.outputTime += increment
                        if self.outputTime >= t:
                            break

            return self.results

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
       # else:
            #TODO: need to account for 'show_labels'




