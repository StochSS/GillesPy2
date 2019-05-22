'''
This method initailizes the state for tau-leaping selections to be made.
Based on Cao, Y.; Gillespie, D. T.; Petzold, L. R. (2006). "Efficient step size selection for the tau-leaping simulation method" (PDF). The Journal of Chemical Physics. 124 (4): 044109. Bibcode:2006JChPh.124d4109C. doi:10.1063/1.2159468. PMID 16460151
'''
def initialize(model, epsilon):

    HOR = {} # Highest Order Reaction of species
    reactants = set()  # a list of all species in the model which act as reactants
    mu_i = {}   # mu_i for each species
    sigma_i = {}  # sigma_i squared for each species
    g_i = {}    # Relative species error allowance denominator
    epsilon_i = {} # Relative error allowance of species
    critical_threshold = 10  # Reactant Population to be considered critical

    #initialize Highest Order Reactions
    for s in model.listOfSpecies:
        HOR[s] = 0

    #Create list of all reactants
    for r in model.listOfReactions:
        #Calculate this reaction's order
        reaction_order = sum(model.listOfReactions[r].reactants.values())
        for reactant, count in model.listOfReactions[r].reactants.items():
            # Build reactant list
            reactants.add(reactant)
            # Initialize mu and sigma for each reactant
            mu_i[reactant] = 0
            sigma_i[reactant] = 0
            # if this reaction's order is higher than previous, set HOR
            if reaction_order > HOR[reactant.name]:
                HOR[reactant.name] = reaction_order
                if count == 2 and reaction_order == 2: g_i[reactant] = lambda x: 2+(1/(x-1))
                elif count == 2 and reaction_order == 3: g_i[reactant] = lambda x: (3/2)*(2+(1/(x-1)))
                elif count == 3: g_i[reactant] = lambda x: (3+(1/(x-1))+(2/(x-2)))
                else: 
                    g_i[reactant] = HOR[reactant.name]
                    epsilon_i[reactant] = epsilon / g_i[reactant]

    # Return components for tau selection
    return HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold

'''
Tau Selection method based on Cao, Y.; Gillespie, D. T.; Petzold, L. R. (2006). "Efficient step size selection for the tau-leaping simulation method" (PDF). The Journal of Chemical Physics. 124 (4): 044109. Bibcode:2006JChPh.124d4109C. doi:10.1063/1.2159468. PMID 16460151
'''
def select(*tau_args):
    
    HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold, model, propensities, curr_state, curr_time, save_time, debug = tau_args
    tau_step = None
    crit_taus = {}
    critical_reactions = []
    critical = False
    critical_tau = 1000 #arbitrarily large value
    non_critical_tau = 0

    #Determine if there are any critical reactions
    for rxn in model.listOfReactions:
        for r, v in model.listOfReactions[rxn].reactants.items():
            if curr_state[r.name] / v < critical_threshold and propensities[rxn] > 0:
                critical_reactions.append(rxn)
                critical = True
                print('CRITICAL: ', r)

    # If a critical reaction is present, estimate tau for a single firing of each
    # critical reaction with propensity > 0
    if critical:
        for rxn in model.listOfReactions:
            if propensities[rxn] > 0:
                crit_taus[rxn] = 1/propensities[rxn]
        critical_tau = min(crit_taus.values())


    for r in g_i:
        #If a reactant's HOR requires >1 of that reactant, evaluate lambda at curr_state
        if callable(g_i[r]):
            g_i[r] = g_i[r](curr_state[r.name])
            epsilon_i[r] = self.epsilon / g_i[r]

    if debug:
        print("curr_state = {", end='')
        for i, s in enumerate(model.listOfSpecies):
            print("'{0}' : {1}, ".format(s, curr_state[s]), end='')
        print("}")

    tau_i = {}  # estimated tau for non-critical reactions
    non_critical_tau = None
    mu_i = {species: 0 for species in model.listOfSpecies.values()}
    sigma_i = {species: 0 for species in model.listOfSpecies.values()}

    for r in model.listOfReactions:
        #For non-critical Reactions
        if not r in critical_reactions:
            #Calculate abs mean and standard deviation for each reactant
            for reactant in model.listOfReactions[r].reactants:
                if reactant not in mu_i:
                    mu_i[reactant] = 0
                if reactant not in sigma_i:
                    sigma_i[reactant] = 0
                mu_i[reactant] += model.listOfReactions[r].reactants[reactant] * propensities[
                    r]  # Cao, Gillespie, Petzold 29a
                sigma_i[reactant] += model.listOfReactions[r].reactants[reactant] ** 2 * propensities[
                    r]  # Cao, Gillespie, Petzold 29b
    for r in reactants:
        calculated_max = epsilon_i[r] * curr_state[r.name]
        #print('calculated max: ', calculated_max)
        max_pop_change_mean = max(calculated_max, 1)
        max_pop_change_sd = max(calculated_max, 1) ** 2
        if mu_i[r] > 0:
            # Cao, Gillespie, Petzold 33
            tau_i[r] = min(
                    max_pop_change_mean / mu_i[r], 
                    max_pop_change_sd / sigma_i[r])

    non_critical_tau = min(tau_i.values())
    tau = non_critical_tau if critical_tau is None else min(non_critical_tau, critical_tau)
    tau_step = min(max(tau, 1e-10), save_time - curr_time)
    return tau_step
