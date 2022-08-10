from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from cobra.io import read_sbml_model
# --- Step 1: Create reference model and model with heterologous reactions for alpha-pienen synthesis in peroxisomes ---

# create reference and engineered model
model = read_sbml_model("Data/iYLI647.xml")

# initilize metabolite objects for heterologous reactions in peroxisomes


# initilize reaction objects for heterologous reactions


# Add metabolites and reaction stoichiometry to reaction objects


# Add reactions to model
reactionlist=[
    
]

for reaction in reactionlist:
    model.add_reaction(reaction)
