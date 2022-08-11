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

aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')
NPP = Metabolite('NPP', formula='C5H9O7P2', name='NPP', compartment='c' )
GPP = model.metabolites.get_by_id("grdp_c")
FPP = model.metabolites.get_by_id("frdp_c")
IPP = model.metabolites.get_by_id("ipdp_c")
DMAPP = model.metabolites.get_by_id("dmpp_c")
Diphosphate = model.metabolites.get_by_id("ppi_c")


def AP_Syn():
    APS_GPP_reaction = Reaction(id="aPinen_gpp_s", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
    APS_NPP_reaction = Reaction(id="aPinen_fpp_s", name="1.0 NPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
    Extract_aPinene =Reaction(id="aPinene_ex", name="extract Alpha-Pinene", subsystem="Cytoplasm", upper_bound=1000.0)
    APS_GPP_reaction.add_metabolites({
        GPP: -1.0,
        aPinene: 1.0,
        Diphosphate: 1.0
    })
    APS_NPP_reaction.add_metabolites({
        NPP: -1.0,
        aPinene: 1.0,
        Diphosphate: 1.0
    })


    Extract_aPinene.add_metabolites({
        aPinene: -1.0
    })

    reactionlist.append(APS_GPP_reaction)
    reactionlist.append(APS_NPP_reaction)
    reactionlist.append(Extract_aPinene)


def GPP_Pathway():
    GPPS_reaction = Reaction(id="GPP_s", name= "1.0 IPP+1.0 DMAPP--> 1.0 GPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound= OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
    mFPPS_gpp_reaction = Reaction(id="FPP_GPP_s",name="1.0 IPP+1.0 DMAPP--> 1.0 GPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND,upper_bound=1000.0)
    mFPPS_fpp_reaction = Reaction(id="FPP_FPP_s",name="1.0 GPP+1.0 DMAPP--> 1.0 FPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND,upper_bound=1000.0)
    
    GPPS_reaction.add_metabolites({
        IPP: -1.0,
        DMAPP: -1.0,
        Diphosphate: 1.0,
        GPP: 1.0
    })

    mFPPS_gpp_reaction.add_metabolites({
        IPP: -1.0,
        DMAPP: -1.0,
        Diphosphate: 1.0,
        GPP: 1.0
    })
    mFPPS_gpp_reaction.gene_reaction_rule='mFPPs'
    
    mFPPS_fpp_reaction.add_metabolites({
        GPP: -1.0,
        DMAPP: -1.0,
        Diphosphate: 1.0,
        FPP: 1.0
    })
    mFPPS_fpp_reaction.gene_reaction_rule='mFPPs'

    reactionlist.append(GPPS_reaction)
    reactionlist.append(mFPPS_gpp_reaction)
    reactionlist.append(mFPPS_fpp_reaction)

def NPP_Pathway():
    SINPPS_reaction = Reaction(id="SINPP_s", name= "1.0 IPP+1.0 DMAPP--> 1.0 NPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound= OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
    SINPPS_reaction.add_metabolites({
        IPP: -1.0,
        DMAPP: -1.0,
        Diphosphate: 1.0,
        NPP: 1.0
    })
    reactionlist.append(SINPPS_reaction)

def erg20_knockdown():
    model.reactions.get_by_id("DMATT").higher_bound=KNOCK_DOWN_HIGHER_BOUND
    model.reactions.get_by_id("GRTT").higher_bound=KNOCK_DOWN_HIGHER_BOUND


#activate functions to build model
GPP_Pathway()
NPP_Pathway()
AP_Syn()
erg20_knockdown() #erg20 knockdown 
#model.genes.YALI0E05753g.knock_out() #erg20 knockout

for reaction in reactionlist:
    model.add_reaction(reaction)

#objective function
model.objective={model.reactions.get_by_id("Biomass_Climit"):GROWTH_OBJECTIVE_COEFFICIENT, model.reactions.get_by_id("aPinene_ex"):APINENE_OBJECTIVE_COEFFICIENT}

solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")
