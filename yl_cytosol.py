from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from cobra.io import read_sbml_model

#constants:
OVEREXPRESSION_LOWER_BOUND = 0.2
KNOCK_DOWN_HIGHER_BOUND = 500
APINENE_OBJECTIVE_COEFFICIENT = 0.2
GROWTH_OBJECTIVE_COEFFICIENT = 0.8

# --- Step 1: Create reference model and model with heterologous reactions for alpha-pienen synthesis in cytosol ---

# create reference and engineered model
model = read_sbml_model("Data/iYLI647.xml")
reactionlist =[]
def smallProcess():
    # initilize metabolite objects for heterologous reactions in cytosol
    aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')
    GPP = model.metabolites.get_by_id("grdp_c")
    FPP = model.metabolites.get_by_id("frdp_c")
    IPP=model.metabolites.get_by_id("ipdp_c")
    DMAPP=model.metabolites.get_by_id("dmpp_c")
    Diphosphate = model.metabolites.get_by_id("ppi_c")

    # initilize reaction objects for heterologous reactions
    GPPS_reaction =Reaction(id="GPP_s", name= "1.0 IPP+1.0 DMAPP--> 1.0 GPP + 1.0 Diphosphate",subsystem="Cytoplasm",upper_bound=1000.0)
    mFPPS_reaction =Reaction(id="FPP_s")
    APS_reaction = Reaction(id="aPinen_s", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", subsystem="Cytoplasm", upper_bound=1000.0)
    Extract_aPinene =Reaction(id="aPinene_ex", name="extract Alpha-Pinene", subsystem="Cytoplasm", upper_bound=1000.0)
    # Add metabolites and reaction stoichiometry to reaction objects

    GPPS_reaction.add_metabolites({
        IPP: -1.0,
        DMAPP: -1.0,
        Diphosphate: 1.0,
        GPP: 1.0
    })

    APS_reaction.add_metabolites({
        GPP: -1.0,
        aPinene: 1.0,
        Diphosphate: 1.0
    })

    Extract_aPinene.add_metabolites({
        aPinene: -1.0
    })
    # Add reactions to model
    reactionlist.append(GPPS_reaction)
    reactionlist.append(APS_reaction)
    reactionlist.append(Extract_aPinene)


smallProcess()

for reaction in reactionlist:
    model.add_reaction(reaction)

#model.objective= "aPinene_ex"
model.objective={model.reactions.get_by_id("Biomass_Climit"):0.2, model.reactions.get_by_id("aPinene_ex"):0.8}

solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")
