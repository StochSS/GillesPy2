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
        for j in reactions:

            if i not in dependent_rxns:
                dependent_rxns[i] = {'dependencies': []}
            if j not in dependent_rxns:
                dependent_rxns[j] = {'dependencies': []}
            if i == j:
                continue

            reactantsI = list(model.listOfReactions[i].reactants.keys())
            reactantsJ = list(model.listOfReactions[j].reactants.keys())
            productsI = list(model.listOfReactions[i].products.keys())

            if model.listOfReactions[i].type == 'customized':
                reactantsI.extend(species_parse(model,model.listOfReactions[i].propensity_function))

            if model.listOfReactions[j].type == 'customized':
                reactantsJ.extend(species_parse(model,model.listOfReactions[j].propensity_function))

            if j not in dependent_rxns[i]['dependencies']:
                if any(elem in reactantsI for elem in reactantsJ):
                    if i not in dependent_rxns[j]['dependencies']:
                        dependent_rxns[j]['dependencies'].append(i)
                    dependent_rxns[i]['dependencies'].append(j)

            if i not in dependent_rxns[j]['dependencies']:
                if any(elem in list(model.listOfReactions[i].products.keys()) for elem in
                       list(model.listOfReactions[j].reactants.keys())):
                    dependent_rxns[j]['dependencies'].append(i)
    return dependent_rxns