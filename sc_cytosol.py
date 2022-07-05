from Data.yeastGEM_io import read_yeast_model, write_yeast_model
from Model.model_utils import summarize_properties
from cobra import Model, Reaction, Metabolite
import optlang

# --- step 1: load scaffold model
model = read_yeast_model() # loads existing SBML yeastGEM model 
# check content of scaffold
# summarize_properties(model)



# --- step 2: adjust model (add new and adjust knocked-down/deleted or overexpressed 
#                       reactions and corresponding metabolites)

# How to add new reaction: Create object and add respective metabolites
# How to 'knock-out' an existing reaction: set upper and lower bound to 0 (= no flux can pass through)
# How to overexpress a reaction: ?

reactions_list = list()

# new metabolites:
aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')

# biological reaction: 1.0 GPP -> 1.0 aPinene
react_gpp_alphapinene = Reaction('r_gpp_apinene')
react_gpp_alphapinene.name = 'GPP -> aPinene'
react_gpp_alphapinene.subsystem = 'Cytoplasm'
react_gpp_alphapinene.lower_bound = 0.  # This is the default
react_gpp_alphapinene.upper_bound = 1000.  # This is the default
react_gpp_alphapinene.add_metabolites({
    model.metabolites.get_by_id("s_0745"): -1.0, # GPP (C10H17O7P2)
    model.metabolites.get_by_id("s_0633"): 1.0, # diphosphate/ppi (HO7P2)
    aPinene: 1.0
})
react_gpp_alphapinene.gene_reaction_rule = '(AlphapineneSnythase)'
#TODO: add associated gene(s)!

# reaction that consumes aPinene
react_alphapinen_con = Reaction('r_apinene_con')
react_alphapinen_con.name = 'aPinene -> '
react_alphapinen_con.subsystem = 'Cytoplasm'
react_alphapinen_con.lower_bound = 0.  # This is the default
react_alphapinen_con.upper_bound = 1000.  # This is the default
react_alphapinen_con.add_metabolites({
    aPinene: -1.0
})

reactions_list.append(react_gpp_alphapinene)
reactions_list.append(react_alphapinen_con)

model.add_reactions(reactions_list)


# --- step 3: set model objective (reaction(s) whose fluxes shall be maximized)
# model.objective = "r_apinene_con"
# mehrere objectives definieren:
model.objective = {model.reactions.get_by_id("r_apinene_con"): 0.5, model.reactions.get_by_id("r_2111"): 0.5}

print(f"\nModel objective expression: {model.objective.expression}")
print(f"Model objective direction: {model.objective.direction}\n")

# adjust lower bounds to simulate overexpression of underlying genes
model.reactions.get_by_id("r_gpp_apinene").lower_bound = 2


# --- step 4: run solver
# declare which solver to use (best one: Gurobipy)
model.solver = 'gurobi'

solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")

# Get any value in the solution
solution.fluxes.get('r_apinene_con')

# save solution object
#solution.fluxes.to_csv("out.csv", index=False)

# summarize solution
print(model.summary())





# last step: save model
#write_yeast_model(model)   # saving, writes SBML file