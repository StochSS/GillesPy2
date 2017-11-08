import gillespy2
from .gillespySolver import GillesPySolver
import random
import math
import numpy as np
<<<<<<< HEAD
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
=======
>>>>>>> 810e9cdd3fd84752cb700dab30d2a3572129f787

class TauLeapingSolver(GillesPySolver):
    """ TODO
    """
<<<<<<< HEAD

=======
>>>>>>> 810e9cdd3fd84752cb700dab30d2a3572129f787
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        self.simulation_data = []
<<<<<<< HEAD
 
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
=======
        listOfAffectedSpecies = {}
        listOfAffectedReactions = {}
        criticalThreshold = 20 # switch this to be a parameter later
        epsilon = 0.5 # switch this to be a parameter later

        # initialize tagList -- a dict of booleans keyed by reactions to indicate if a reaction is critical
        tagList = [False]*len(model.listOfReactions)
        tagList = zip(model.listOfReactions, tagList)
        tagList = dict(tagList)

        # create dict of species affected by a given reaction
        for reaction in model.listOfReactions:
            for species in model.listOfReactions[reaction].reactants:
                affectedSpecies[reaction].append(species)

        # create dict of reactions affected by a given species
        for species in model.listOfSpecies:
            for reaction in model.listOfReactions:
                if species in model.listOfReactions[reaction].reactants:
                    affectedReactions[species].append(reaction)

        # create stoichiometry vector and square it
        for reaction in model.listOfReactions:
            stoichiometry[reaction] = {}
            for species in model.listOfSpecies:
                if species in model.listOfReactions[reaction].products:
                    stoichiometry[reaction].append(species: model.listOfReactions[reaction].products[species])
                else:
                    stoichiometry[reaction].append(species: 0)
                if species in model.listOfReactions[reaction].reactants:
                    stoichiometry[reaction][species] -= model.listOfReactions[reaction].reactants[epcies]
        squaredVj = np.square(stoichiometry)

        # run the algorithm
        return simmulate(self)


    # updates tagLists
    def updateTagLists(self):
        
        # from critical species to noncritical species
        for species in listOfCriticalSpecies:
            if curr_state[species] > criticalThreshold:
                listOfNonCriticalSpecies.append(species)
                for reaction in affectedReactions[species]
                    changeTag = True
                    for affectedSpecies in listOfAffectedSpecies[reaction]:
                        if curr_state[affectedSpecies] <= criticalThreshold:
                            changeTag = False
                            break
                    if changeTag:
                        tagList[reaction] = False
            del listOfCriticalSpecies[species]

        # from non-critical species to critical species
        for species in listOfNonCriticalSpecies:
            if curr_state[species] <= crticialThreshold:
                listOfCriticalSpecies.append(species)
                for affectedReaction in listOfAffectedReactions[species]:
                    tagList[affectedReaction] = True
                del listOfNonCriticalSpecies[species]
                    
    # selects tau
    def selectTau(self, nonCriticalStepsize, criticalStepsize):
        # for critical reactions
        criticalPropensitySum = 0
        for reaction in tagList:
            if tagList[reaction] == true:
                criticalPropensitySum += propensities[reaction]
        if criticalPropensitySum != 0:
            criticalStepsize = np.random.exponential(1.0/criticalPropensitySum) # possible error
        else:
            criticalStepsize = -1

        # for noncritical reactions
        # update mu and sigmaSquared
        mu = np.dot(propensities, stoichiometry)
        sigmaSquared = np.dot(propensities, squaredVj)

        nonCriticalStepsize = sys.maxsize
        for species in listOfNonCriticalSpecies:
            numerator = epsilon*curr_state[species]/g
            numerator = max(numerator, 1.0)
            temp1 = numerator/abs(mu[species])
            temp2 = numerator*numerator/sigmaSquared[species]
            nonCriticalStepsize = min(nonCriticalStepsize, temp1)
            nonCriticalStepsize = min(nonCriticalStepsize, temp2)

    # selects the reactions
    def selectReactions(self, leapsize, runCritical):
        if runCritical: # handle both critical and noncritical reactions
            rand = 0
            while rand == 0:
                rand = random.uniform(0, criticalPropensitySum)
            jsum = 0

            for reaction in model.listOfReactions:
                # for critical
                if tagList[reaction] == true:
                    previousReactionCounts[reaction] = 0
                    if (jsum < rand):
                        previousReaction = reaction
                        jsum += propensity[reaction]
                # for noncritical
                else:
                    previousReactionCounts[reaction] = np.random.poisson(leapSize*propensities[reaction])
            return previousReaction
                else:
                    #for noncritical
                    previousReactionCounts[reaction] = np.random.poisson(leapSize*propensities[reaction])
        else: #handle only noncritical reactions
            for reaction in model.listOfReactions:
                if tagList[reaction] == true:
                    previousReactionCounters[reaction] = 0
                else
                    previousReactionCounts[reaction] = np.random.poisson(leapSize*propensities[reaction])
            return None


    # fires the reaction
    def fireReactions(self, reaction):
        populationChange = np.dot(previousReactionCounts, stoich)
        curr_state += populationChange
        if reaction != None:
            critical_fireReaction(reaction)
    
        # check for negative populations
        negativePopulation = False
        for species in model.listOfSpecies:
            if curr_state[species] < 0.0:
                negativePopulation = True
                break
        if negativePopulation:
            curr_state -= populationChange
            if reaction != None
                # the heck is this?
                critical_rollBack(reaction)
            return False
        else
            calculateAllPropensities()
            return True

    # calculates the propensities (should eventually inherrit from SSA_Solver, correct?)
    def (calculateAllPropensities(self):
        prop_sum = 0
        for reaction in model.listOfReactions:
            propensity[reaction] = eval(model.listOfReactions[reaction].propensity_function, curr_state)
            prop_sum += propensities[reaction]
        


    # fires a critical reaction
    def critical_fireReaction(self, reaction):
        if (reaction == None):
            return False
        else:
            # check this!
            for species in stoichiometry[reaction]:
                curr_state[species] += stoichiometry[reaction][species]
            curr_state += stoichiometry[reaction]
            return True


    # runs the tau-leaping algorithm
    def simmulate(self):
        
        curr_state = {}
        propensities = {}
        results = {}
        poissonValues = {}
        expectedChanges = {}
        criticalStepsize = 0
        nonCriticalStepsize = 0

        for traj_num in range(number_of_trajectories):
            # put this into a different function?
            for species in model.listOfSpecies:   #Initialize Species population
                curr_state[species] = model.listOfSpecies[species].initial_value    
                results[species]=[]
            listOfCriticalSpecies = {}
            listOfNonCriticalSpecies = {}
            curr_state['vol'] = model.volume
            results['time'] = []
            for time in range (0, t, increment):
                results['time'].append(time)
            curr_time = 0
            currentTimeInterval = 0
            reactionsLastLeap = sys.maxsize

            for p in model.listOfParameters:
                curr_state[p] = model.listOfParameters[p].value        
            
            previousReactionCounts = {}

            # do the algorithm 
            while curr_time < t
                while (curr_time >= results['time'][currentTimeInterval]) and (currentTimeInterval < len(results['time']):
                    for species in model.listOfSpecies:
                        results[species].append(curr_state[species])
                        currentTimeInterval += 1
                
                if false: # if reactionsLastLeap < threshold) {do the normal SSA}
                    # do SSA

                else:
                    # do tau-leaping
                    updateTagLists()
                    selectTau(nonCriticalStepsize, criticalStepsize)
                    if (nonCriticalStepsize < criticalStepsize) or (criticalStepsize == -1):
                        tau = nonCriticalStepsize
                        runCritical = False
                    else:
                        tau = criticalStepsize
                        runCritical = True
                    if (curr_time + tau > results['time'][currentTimeInterval]):
                        tau = results['time'][currentTimeInterval] - curr_time
                        nextTime = curr_time + tau
                        runCritical = False
                    else:
                        nextTime = curr_time + tau
                    criticalReaction = selectReactions(tau, runCritical)
                    reactionsLastLeap = np.linalg.norm(previousReactionCounts)
                    if runCritical:
                        reactionsLastLeap += reactionsLastLeap
                        totalReactionsDuringLeaps += reactionsLastLeap
                    if fireReactions(criticalReaction):
                        curr_time = nextTime
                        failedLeaps = 0
                        totalLeapsTaken += totalLeapsTaken
                    else
                        failedLeaps += failedLeaps
                        # put in some errors homie

        while (curr_time >= results['time'][currentTimeInterval]):
            for species in model.listOfSpecies:
                results[species].append(curr_state[species])
            currentTimeInterval += currentTimeInterval
            if currentTimeInterval >= len(results['time']
                break

        return results
>>>>>>> 810e9cdd3fd84752cb700dab30d2a3572129f787

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
       # else:
            #TODO: need to account for 'show_labels'
<<<<<<< HEAD




=======
>>>>>>> 810e9cdd3fd84752cb700dab30d2a3572129f787
