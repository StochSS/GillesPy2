import gillespy2
from .gillespySolver import GillesPySolver
import random
import math
import numpy as np
import sys

# selects Tau based on algorithm in Cao et al. 2006
# def selectTau(self, propensities):
#    hello = updateMuAndSigmaSquared(propensities)
#    return 0.05

# def updateMuAndSigmaSquared(self, propensities):
#    for reaction in model.listOfReactions:
#        mu[reaction] = propensities[reaction]*(model.listOfReactions[reaction].products - model.listOfReactions[reaction].reactants)
        #mu = propensities[reaction]
        #hi = model.listOfReactions[reaction].products - model.listOfReactions[reaction].reactants

#    return mu

class TauLeapingSolver(GillesPySolver):
    """ TODO
    """

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        self.simulation_data = []
 
        curr_state = {}
        propensities = {}
        results = {}
        listOfAffectedReactions = {}
        isCritical = {}    # keyed by reaction
        criticalThreshold = 25  # threshold for when a species becomes critical
        previousReactionCounts = {}
        populationChange = {}

        for traj_num in range(number_of_trajectories):
            # put into own function (initialize())?
            curr_state['vol'] = model.volume
            results['time'] = []
            currentTime = 0
            outputTime = 0      
            nextTime = 0
            vj = []
            listOfCriticalSpecies = []
            listOfNonCriticalSpecies = []
            failedLeaps = 0

            for species in model.listOfSpecies:
                # Initialize species populations
                curr_state[species] = model.listOfSpecies[species].initial_value    
                results[species]=[]
                listOfAffectedReactions[species] = [] 
                listOfCriticalSpecies.append(species)

            for reaction in model.listOfReactions:
                # build a list of which reactions are affected by each species
                for reactant in model.listOfReactions[reaction].reactants:
                    listOfAffectedReactions[str(reactant)].append(reaction)
                # default every reaction to be critical
                isCritical[reaction] = True
                previousReactionCounts[reaction] = []

            for parameter in model.listOfParameters:
                curr_state[parameter] = model.listOfParameters[parameter].value        
        
            # run the algorithm
            while(currentTime < t):
                while currentTime >= outputTime and outputTime < t:
                    results['time'].append(outputTime)
                    for species in model.listOfSpecies:
                        results[species].append(curr_state[species])
                    outputTime += increment

                # update tag lists --> determine which reactions are critical
                # should be own function
                # from critical list to non-critical list
                for criticalSpecies in listOfCriticalSpecies:
                    if curr_state[criticalSpecies] > criticalThreshold:
                        listOfNonCriticalSpecies.append(criticalSpecies)
                        # modify the tag for every reaction affected by this species
                        for affectedReaction in listOfAffectedReactions[criticalSpecies]:
                            changeTag = True
                            # check if other reactants for this reaction are critical
                            for affectedSpecies in model.listOfReactions[affectedReaction].reactants:
                                if affectedSpecies in listOfCriticalSpecies:
                                    changeTag = False
                                    break
                            if changeTag:
                                isCritical[affectedReaction] = False
                        # erase criticalSpecies from critical list
                        listOfCriticalSpecies.remove(criticalSpecies)
                # from non-critical list to critical list
                for nonCriticalSpecies in listOfNonCriticalSpecies:
                    if curr_state[nonCriticalSpecies] <= criticalThreshold:
                        listOfCriticalSpecies.append(nonCriticalSpecies)
                        # modify tag for each reaction affected by this species
                        for affectedReaction in listOfAffectedReactions[nonCriticalSpecies]:
                            isCritical[affectedReaction] = True
                        listOfNonCriticalSpecies.remove(nonCriticalSpecies)

                # evaluate propensities
                for reaction in model.listOfReactions: 
                    propensities[reaction] = eval(model.listOfReactions[reaction].propensity_function, curr_state)
   
                # select tau. SHOULD BE OWN FUNCTION
                # for critical reactions
                criticalStepsize = 0
                criticalPropensitySum = 0
                for reaction in model.listOfReactions:
                    if isCritical[reaction]:
                        criticalPropensitySum += propensities[reaction]
                if criticalPropensitySum != 0:
                    criticalStepsize = -1*math.log(random.random())/criticalPropensitySum
                else:
                    criticalStepSize = -1
                # for noncritical reactions
                noncriticalStepsize = sys.maxsize
                # generate mu and sigmaSquared
#                muSum = 0
#                for reaction in model.listOfReactions:
#                    for product in model.listOfReactions[reaction].products:
#                        muSum += model.listOfReactions[reaction].products[product]*propensities[reaction]
#                    for reactant in model.listOfReactions[reaction].reactants:
#                        muSum -= model.listOfReactions[reaction].reactants[reactant]*propensities[reaction]
                # if we don't need to worry about critical reactions
                if noncriticalStepsize < criticalStepsize or criticalStepsize == -1:
                    tau = noncriticalStepsize
                    runCritical = False
                # otherwise, we do need to worry about critical reactions
                else:
                    tau = criticalStepsize
                    runCritical = True
                tau = 0.5
                
                # avoid leaping past next output time
                if currentTime + tau > outputTime and currentTime <= t:
                    tau = outputTime - currentTime
                    nextTime = outputTime
                    outputTime += increment    
                else:
                    nextTime = currentTime + tau
                
                # select reactions --> should be function selectReactions(tau, runCritical)
                if runCritical:
                    # handle both critical and noncritical reactions
                    # for critical reactions
                    criticalReaction = None
                    # generate uniform random number between 0 and criticalPropensitySum
                    reactionNum = random.uniform(0, criticalPropensitySum)
                    jsum = 0
                    for reaction in model.listOfReactions:
                        # for critical reactions
                        if isCritical[reaction]:
                            previousReactionCounts[reaction] = 0
                            if jsum < reactionNum:
                               criticalReaction = reaction
                               jsum += propensities[reaction]
                        # for noncritical reactions
                        else:
                            previousReactionCounts[reaction] = np.random.poisson(propensities[reaction]*tau)
                else:
                    # handle only noncritical reactions
                    for reaction in model.listOfReactions:
                        if isCritical[reaction]:
                            previousReactionCounts[reaction] = 0
                        else:
                            previousReactionCounts[reaction] = np.random.poisson(propensities[reaction]*tau)
                #reactionsLastLeap = norm of previousReactionCounts -- for after selectTau is working

                # should be fireReactions(previousReaction)
                for species in model.listOfSpecies:
                    populationChange[species] = 0
                # calculate the population changes
                for reaction in model.listOfReactions:
                    for reactant in model.listOfReactions[reaction].reactants:
                        populationChange[str(reactant)] -= model.listOfReactions[reaction].reactants[reactant]*previousReactionCounts[reaction]
                    for product in model.listOfReactions[reaction].products:
                        populationChange[str(product)] += model.listOfReactions[reaction].products[product]*previousReactionCounts[reaction]
                # append changes to curr_state
                for species in model.listOfSpecies:
                    curr_state[str(species)] += populationChange[str(species)]

                # if we need to handle a critical reaction. Should be function fireCritical
                if criticalReaction is not None:
                    for reactant in model.listOfReactions[criticalReaction].reactants:
                        curr_state[str(reactant)] -= model.listOfReactions[criticalReaction].reactants[reactant]
                    for product in model.listOfReactions[criticalReaction].products:
                        curr_state[str(product)] += model.listOfReactions[criticalReaction].products[product]

                # check for negative populations
                negativePopulation = False
                leapFailed = False
                for species in model.listOfSpecies:
                    if curr_state[str(species)] < 0:
                        negativePopulation = True
                        break
                if negativePopulation:
                    # we need to undo what we just did
                    for species in model.listOfSpecies:
                        curr_state[str(species)] -= populationChange[str(species)]
                    if criticalReaction is not None:
                        # reverse criticalReaction firing. Should be criticalRollBack()
                        for reactant in model.listOfReactions[criticalReaction].reactants:
                            curr_state[str(reactant)] += model.listOfReactions[criticalReaction].reactants[reactant]
                        for product in model.listOfReactions[criticalReaction].products:
                            curr_state[str(product)] -= model.listOfReactions[criticalReaction].products[product]
                    leapFailed = True
                else:
                    # calculateAllPropensities?
                    leapFailed = False

                if leapFailed:
                    failedLeaps += 1
                else:
                    currentTime = nextTime
                    failedLeaps = 0

                # record output. (Is this in the right place?)
                results['time'].append(outputTime)
                for species in model.listOfSpecies:
                    results[species].append(curr_state[species])

            return results

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
       # else:
            #TODO: need to account for 'show_labels'




