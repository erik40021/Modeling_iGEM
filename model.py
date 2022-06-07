from Data.yeastGEM_io import read_yeast_model, write_yeast_model
from Model.model_utils import summarize_properties
from cobra import Model, Reaction, Metabolite

# --- step 1: load scaffold model
model = read_yeast_model() # loads existing SBML yeastGEM model 
# check content of scaffold
summarize_properties(model)



# --- step 2: adjust model (add new and adjust knocked-down/deleted or overexpressed 
#                       reactions and corresponding metabolites and genes)

reactions_list = list()
react_ipp_gpp = Reaction('r_ipp_mFPS144_gpp')
# biological reaction: 1.0 IPP + 1.0 DMAPP -> 1.0 GPP
react_ipp_gpp.name = 'mFPS144 IPP -> GPP '
react_ipp_gpp.subsystem = 'Cytoplasm'
react_ipp_gpp.lower_bound = 0.  # This is the default
react_ipp_gpp.upper_bound = 1000.  # This is the default
# summarize_properties(model) # check content again after adding reactions
DMAPP = Metabolite(
    'DMAPP',
    formula='C5H9O7P2',
    name='Dimethylallylpyrophosphat',
    compartment='c')
react_ipp_gpp.add_metabolites({
    model.metabolites.get_by_id("s_0943"): -1.0, # IPP id: s_0943
    model.metabolites.get_by_id("s_0745"): 1.0,  # GPP id: s_0745
    DMAPP: -1.0
})

#TODO: correct reaction (in reality: three different reactions!)
react_mevalonate_dmapp = Reaction('r_mevalonate_dmapp') #mevalonate id: s_0028
# biological reaction: 1.0 Mevalonate -> 1.0 DMAPP (Simplified!!)
react_mevalonate_dmapp.name = '_ Mevalonate -> DMAPP'
react_mevalonate_dmapp.subsystem = 'Cytoplasm'
react_mevalonate_dmapp.lower_bound = 0.  # This is the default
react_mevalonate_dmapp.upper_bound = 1000.  # This is the default

react_mevalonate_dmapp.add_metabolites({
    model.metabolites.get_by_id("s_0028"): -1.0,
    DMAPP: 1.0
})
reactions_list.append(react_ipp_gpp)
reactions_list.append(react_mevalonate_dmapp)


model.add_reactions(reactions_list)






# --- step 3: set model objective (reaction(s) whose fluxes shall be maximized)
# TODO
model.objective = "r_ipp_mFPS144_gpp" #vorläufig, am Ende: Reaktion, die zu Verbenon führt

print(f"\nModel objective expression: {model.objective.expression}")
print(f"Model objective direction: {model.objective.direction}\n")




# --- step 4: run solver
solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")

# 
#x = model.slim_optimize() # if only one reaction flux (default: objective value) is needed
#print(x)

# summarize solution
print(model.summary())





# last step: save model
write_yeast_model(model)   # saving, writes SBML file