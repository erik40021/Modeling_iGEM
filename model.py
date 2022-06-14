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
# dmapp = Metabolite('DMAPP', formula='C5H9O7P2', name='Dimethylallylpyrophosphat', compartment='c')
# ------> Already existing metabolite s_1376
# TODO: correct formula of aPinene
aPinene = Metabolite('aPinene', formula='C10H17O7P2', name='Alphapinene', compartment='c')

# biological reaction: 1.0 IPP + 1.0 DMAPP -> 1.0 GPP
# react_ipp_gpp = Reaction('r_ipp_mFPS144_gpp')
# react_ipp_gpp.name = 'mFPS144 IPP -> GPP '
# react_ipp_gpp.subsystem = 'Cytoplasm'
# react_ipp_gpp.lower_bound = 0.  # This is the default
# react_ipp_gpp.upper_bound = 1000. # This is the default
# react_ipp_gpp.add_metabolites({
#     model.metabolites.get_by_id("s_0943"): -1.0, # IPP id: s_0943
#     model.metabolites.get_by_id("s_0745"): 1.0,  # GPP id: s_0745
#     dmapp: -1.0
# }) ---------> Already existing reaction r_0355!


# biological reaction: 1.0 Mevalonate -> 1.0 DMAPP (Simplified!!)
#TODO: correct reaction (in reality: three different reactions!) AND search for possible existing reactions
react_mevalonate_dmapp = Reaction('r_mevalonate_dmapp') #mevalonate id: s_0028
react_mevalonate_dmapp.name = 'Mevalonate -> DMAPP'
react_mevalonate_dmapp.subsystem = 'Cytoplasm'
react_mevalonate_dmapp.lower_bound = 0.  # This is the default
react_mevalonate_dmapp.upper_bound = 1000.  # This is the default
react_mevalonate_dmapp.add_metabolites({
    model.metabolites.get_by_id("s_0028"): -1.0,
    model.metabolites.get_by_id("s_1376"): 1.0
})
# TODO: correct stochiometry!

# biological reaction: 1.0 GPP -> 1.0 aPinene
react_gpp_alphapinene = Reaction('r_gpp_apinene')
react_gpp_alphapinene.name = 'GPP -> aPinene'
react_gpp_alphapinene.subsystem = 'Cytoplasm'
react_gpp_alphapinene.lower_bound = 0.  # This is the default
react_gpp_alphapinene.upper_bound = 1000.  # This is the default
react_gpp_alphapinene.add_metabolites({
    model.metabolites.get_by_id("s_0745"): -1.0,
    aPinene: 1.0
})
react_gpp_alphapinene.gene_reaction_rule = '(AlphapineneSnythase)'
#TODO: add associated gene(s)!


reactions_list.append(react_mevalonate_dmapp)
reactions_list.append(react_gpp_alphapinene)
#reactions_list.append(react_gpp_fpp)

model.add_reactions(reactions_list)

# model.reactions.get_by_id("r_0462").bounds = (0.0, 10.0) # knock out reaction
# alternative: use reaction.knock_out()



# --- step 3: set model objective (reaction(s) whose fluxes shall be maximized)
# TODO
#model.objective = "r_gpp_apinene" #vorläufig, am Ende: Reaktion, die zu Verbenon führt

# mehrere objectives definieren:
# model.objective = {model.reactions.get_by_id("r_ipp_mFPS144_gpp"): 0.5, model.reactions.get_by_id("r_2111"): 0.5}

print(f"\nModel objective expression: {model.objective.expression}")
print(f"Model objective direction: {model.objective.direction}\n")




# --- step 4: run solver
# declare which solver to use (best one: Gurobipy)
model.solver = 'gurobi'

solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")

# save solution object
#solution.fluxes.to_csv("out.csv", index=False)


# 
#x = model.slim_optimize() # if only one reaction flux (default: objective value) is needed
#print(x)

# summarize solution
print(model.summary())





# last step: save model
#write_yeast_model(model)   # saving, writes SBML file