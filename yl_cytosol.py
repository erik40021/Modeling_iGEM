from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from cobra.io import read_sbml_model
# --- Step 1: Create reference model and model with heterologous reactions for alpha-pienen synthesis in cytosol ---

# create reference and engineered model
model = read_sbml_model("Data/iYLI647.xml")

# initilize metabolite objects for heterologous reactions in cytosol
aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')
PPi = model.metabolites.get_by_any("ppi_c")
GPP = model.metabolites.get_by_id("grdp_c")

# initilize reaction objects for heterologous reactions

APS_reaction = Reaction(id="APinen_s", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", subsystem="Cytoplasm", upper_bound=1000.0)
Extract_APinene =Reaction(id="APinen_ex", name="extract Alpha-Pinene", subsystem="Cytoplasm", upper_bound=1000.0)
# Add metabolites and reaction stoichiometry to reaction objects

APS_reaction.add_metabolites({
    GPP: -1.0,
    aPinene: 1.0,
    PPi: 1.0
})

Extract_APinene.add_metabolites({
    aPinene: -1.0
})
# Add reactions to model
reactionlist=[
APS_reaction,
Extract_APinene
]

for reaction in reactionlist:
    model.add_reaction(reaction)

model.objective= "Extract_APinene"

solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")
