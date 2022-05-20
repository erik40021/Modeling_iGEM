from Data.yeastGEM_io import read_yeast_model
from Model.model_utils import summarize_properties

# step 1: load scaffold model
model = read_yeast_model() # loads existing SBML yeastGEM model 


# step 2: check model content
summarize_properties(model)

# step 3: adjust model
#TODO




# last step: save model
#io.write_yeast_model(model)   # saving, writes SBML file