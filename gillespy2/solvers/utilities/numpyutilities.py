import ast
import gillespy2


def species_parse(model,custom_prop_fun):
    parsed_species = []
    class SpeciesParser(ast.NodeTransformer):
        def visit_Name(self, node):
            if isinstance(model.get_element(node.id),gillespy2.core.Species):
                parsed_species.append(model.get_element(node.id))


    expr = custom_prop_fun
    expr = ast.parse(expr,mode='eval')
    expr = SpeciesParser().visit(expr)
    return parsed_species


def dependency_grapher(model,reactions):
    dependent_rxns = {}
    for i in reactions:
        cust_spec = []
        if model.listOfReactions[i].type == 'customized':
            cust_spec = (species_parse(model, model.listOfReactions[i].propensity_function))

        for j in reactions:

            if i not in dependent_rxns:
                dependent_rxns[i] = {'dependencies': []}
            if j not in dependent_rxns:
                dependent_rxns[j] = {'dependencies': []}
            if i == j:
                continue

            reactantsI = list(model.listOfReactions[i].reactants.keys())
            reactantsJ = list(model.listOfReactions[j].reactants.keys())

            if j not in dependent_rxns[i]['dependencies']:
                if any(elem in reactantsI for elem in reactantsJ):
                    if i not in dependent_rxns[j]['dependencies']:
                        dependent_rxns[j]['dependencies'].append(i)
                    dependent_rxns[i]['dependencies'].append(j)

            if i not in dependent_rxns[j]['dependencies']:
                if any(elem in list(model.listOfReactions[i].products.keys()) for elem in
                       list(model.listOfReactions[j].reactants.keys())):
                    dependent_rxns[j]['dependencies'].append(i)

            if cust_spec:
                if any(elem in cust_spec for elem in list(model.listOfReactions[j].reactants)) or any\
                            (elem in cust_spec for elem in list(model.listOfReactions[j].products)):
                    dependent_rxns[i]['dependencies'].append(j)

    return dependent_rxns



