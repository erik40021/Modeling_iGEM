from Data.yeastGEM_io import read_yeast_model
from Model.model_utils import summarize_properties
from cobra import Model, Reaction, Metabolite

# --- step 1: load scaffold model
model = read_yeast_model() # loads existing SBML yeastGEM model 
# check content of scaffold
summarize_properties(model)

# --- step 2: adjust model (add new and adjust knocked-down/deleted or overexpressed 
#                       reactions and corresponding metabolites and genes)
# TODO
r1 = Reaction("Test_reaction1")
r1.name = "test"
# summarize_properties(model) # check content again after adding reactions

# --- step 3: set model objective (reaction(s) of which fluxes shall be maximized)
# TODO
# model.objective = "r_0001"

print(f"\nModel objective expression: {model.objective.expression}")
print(f"Model objective direction: {model.objective.direction}\n")

# --- step 4: run solver
solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")

# 
x = model.slim_optimize() # if only one reaction flux (default: objective value) is needed
print(x)

# summarize solution
print(model.summary())

# example in-detail summary:
# model.metabolites.nadh_c.summary()
# model.metabolites.atp_c.summary()


# last step: save model
#io.write_yeast_model(model)   # saving, writes SBML file