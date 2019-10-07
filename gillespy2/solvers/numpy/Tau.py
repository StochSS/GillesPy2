'''
This Python module contains the initialization and selection methods for the Tau-Leaping method described in Cao, Y.; Gillespie, D. T.; Petzold, L. R. (2006). "Efficient step size selection for the tau-leaping simulation method" (PDF). The Journal of Chemical Physics. 124 (4): 044109. Bibcode:2006JChPh.124d4109C. doi:10.1063/1.2159468. PMID 16460151.

This module is for use in the basic_tau_leaping_solver and basic_tau_hybrid solver only.
'''

def initialize(model, epsilon):
    '''
    This method initailizes the state for tau-leaping selections to be made.
    Based on Cao, Y.; Gillespie, D. T.; Petzold, L. R. (2006). "Efficient step size selection for the tau-leaping simulation method" (PDF). The Journal of Chemical Physics. 124 (4): 044109. Bibcode:2006JChPh.124d4109C. doi:10.1063/1.2159468. PMID 16460151
    '''

    HOR = {} # Highest Order Reaction of species
    reactants = set()  # a list of all species in the model which act as reactants
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
            # if this reaction's order is higher than previous, set HOR
            if reaction_order > HOR[str(reactant)]:
                HOR[str(reactant)] = reaction_order
                if count == 2 and reaction_order == 2: g_i[str(reactant)] = lambda x: 2+(1/(x-1))
                elif count == 2 and reaction_order == 3: g_i[str(reactant)] = lambda x: (3/2)*(2+(1/(x-1)))
                elif count == 3: g_i[str(reactant)] = lambda x: (3+(1/(x-1))+(2/(x-2)))
                else: 
                    g_i[str(reactant)] = HOR[str(reactant)]
                    epsilon_i[str(reactant)] = epsilon / g_i[str(reactant)]

    # Return components for tau selection
    return HOR, reactants, g_i, epsilon_i, critical_threshold

def select(*tau_args):
    '''
    Tau Selection method based on Cao, Y.; Gillespie, D. T.; Petzold, L. R. (2006). "Efficient step size selection for the tau-leaping simulation method" (PDF). The Journal of Chemical Physics. 124 (4): 044109. Bibcode:2006JChPh.124d4109C. doi:10.1063/1.2159468. PMID 16460151
    '''
    
    HOR, reactants, g_i, epsilon_i, epsilon, critical_threshold, \
    model, propensities, curr_state, curr_time, save_time = tau_args
    tau_step = 0
    crit_taus = {} # Estimated time to single-firing of critical reactions
    critical_reactions = [] # List of critical reactions at this step
    critical = False # system-wide flag, true when any reaction is critical
    critical_tau = 0 # holds the smallest tau time for critical reactions
    non_critical_tau = 0 # holds the smallest tau time for non-critical reactions
    tau = 0

    #Determine if there are any critical reactions
    for rxn in model.listOfReactions:
        for r, v in model.listOfReactions[rxn].reactants.items():
            if curr_state[r.name] / v < critical_threshold and propensities[rxn] > 0:
                critical_reactions.append(rxn)
                critical = True
        non_critical_reactions = [rxn for rxn in model.listOfReactions if not rxn in critical_reactions]

    # If a critical reaction is present, estimate tau for a single firing of each
    # critical reaction with propensity > 0, and take the smallest tau
    if critical:
        for rxn in model.listOfReactions:
            if propensities[rxn] > 0:
                crit_taus[rxn] = 1/propensities[rxn]
        critical_tau = min(crit_taus.values())


    # If a reactant's HOR requires >1 of that reactant, evaluate lambda at curr_state
    for r in g_i:
        if callable(g_i[str(r)]):
            g_i[str(r)] = g_i[str(r)](curr_state[str(r)])
            epsilon_i[str(r)] = epsilon / g_i[str(r)]

    tau_i = {}  # estimated tau for non-critical reactions
    mu_i = {str(species): 0 for species in model.listOfSpecies.values()}
    sigma_i = {str(species): 0 for species in model.listOfSpecies.values()}

    for r in non_critical_reactions:
        #Calculate abs mean and standard deviation for each reactant
        for reactant in model.listOfReactions[r].reactants:
            mu_i[str(reactant)] -= model.listOfReactions[r].reactants[reactant] * propensities[
                r]  # Cao, Gillespie, Petzold 32a
            sigma_i[str(reactant)] += model.listOfReactions[r].reactants[reactant] ** 2 * propensities[
                r]  # Cao, Gillespie, Petzold 32b
        #Calculate abs mean and standard deviation for each product
        for product in model.listOfReactions[r].products:
            mu_i[str(product)] += model.listOfReactions[r].products[product] * propensities[
                r]  # Cao, Gillespie, Petzold 32a
            sigma_i[str(product)] += model.listOfReactions[r].products[product] ** 2 * propensities[
                r]  # Cao, Gillespie, Petzold 32b

    for spec in model.listOfSpecies.values():
        if mu_i[str(spec)] != 0 and str(spec) in epsilon_i:
            calculated_max = epsilon_i[str(spec)] * curr_state[spec.name]
            max_pop_change_mean = max(calculated_max, 1)
            max_pop_change_sd = max(calculated_max, 1) ** 2
            # Cao, Gillespie, Petzold 33
            tau_i[str(spec)] = min(
                    abs(max_pop_change_mean / mu_i[str(spec)]), 
                    max_pop_change_sd / sigma_i[str(spec)])

    if len(tau_i) > 0: non_critical_tau = min(tau_i.values())

    # If all reactions are non-critical, use non-critical tau.
    if not critical:
        tau = non_critical_tau
    # If all rxns are critical, use critical tau.
    elif len(tau_i) == 0:
        tau = critical_tau
    # If there are both critical and non-critical reactions,
    # take the shortest tau between critical and non-critical.
    else:
        tau = min(non_critical_tau, critical_tau)
    # If selected tau exceeds save time, integrate to save time
    if tau > 0: tau = max(tau, 1e-10) # set minimum to prevent integration errors
    if tau > 0:
        tau_step = min(tau, save_time - curr_time)
    else:
        tau_step = save_time - curr_time
    return tau_step, mu_i, sigma_i
