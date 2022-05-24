from Data.yeastGEM_io import read_yeast_model
from Model.model_utils import summarize_properties
from cobra import Model, Reaction, Metabolite

# step 1: load scaffold model
model = read_yeast_model() # loads existing SBML yeastGEM model 


# step 2: check model content
summarize_properties(model)

# step 3: adjust model
#TODO
r1 = Reaction("Test_reaction1")
r1.name = "test"
summarize_properties(model)


# last step: save model
#io.write_yeast_model(model)   # saving, writes SBML file