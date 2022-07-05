from cobra import Model, Reaction, Metabolite, read_sbml_model
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
# --- Step 1: Create reference model and model with heterologous reactions for alpha-pienen synthesis in peroxisomes ---

# create reference and engineered model
model = read_sbml_model("Data/iYLl647.xml")

# initilize metabolite objects for heterologous reactions in peroxisomes
Acetyl_CoA = model.metabolites.get_by_id("M_accoa_x")
Acetoacetyl_CoA = Metabolite(id="M_aacoa_x", formula="C25H36N7O18P3S", name="acetoacetyl-CoA[p]", charge=-4, compartment="x") 
CoA = model.metabolites.get_by_id("M_coa_x")
HMG_CoA = Metabolite(id="M_hmgcoa_x", formula="C27H39N7O20P3S", name="3-hydroxy-3-methylglutaryl-CoA[p]", charge=-5, compartment="x")
NADPH = model.metabolites.get_by_id("M_nadph_x")
NADP = model.metabolites.get_by_id("M_nadp_x")
H_p = model.metabolites.get_by_id("M_h_x")
Mevalonate = Metabolite(id="M_mev_x", formula="C6H11O4", name="(R)-mevalonate[p]", charge=-1, compartment="x")
Mevalonate_P = Metabolite(id="M_5pmev_x", formula="C6H10O7P", name="(R)-5-phosphomevalonic acid[p]", charge=-3, compartment="x")
Mevalonate_PP = Metabolite(id="M_5dpmev_x", formula="C6H10O10P2", name="(R)-5-diphosphomevalonic acid[p]", charge=-4, compartment="x")
ATP = model.metabolites.get_by_id("M_atp_x")
ADP = model.metabolites.get_by_id("M_adp_x")
IPP = Metabolite(id="M_ipdp_x", formula="C5H9O7P2", name="isopentyl diphosphate[p]", charge=-3, compartment="x")
DMAPP = Metabolite(id="M_dmpp_x", formula="C5H9O7P2", name="dimethylallyl diphosphate[p]", charge=-3, compartment="x")
GPP = Metabolite(id="M_grdp_x", formula="C10H17O7P2", name="geranyl diphosphate[p]", charge=-3, compartment="x")
Alpha_pinene = Metabolite(id="M_0000_x", formula="C6H16", name="(+)-alpha-pinene[p]", charge=0, compartment="x")#metabolite not excisting in any other compartment
Phosphate = model.metabolites.get_by_id("M_pi_x")
CO2 = model.metabolites.get_by_id("M_co2_x")
Diphosphate = model.metabolites.get_by_id("M_ppi_x")

# initilize reaction objects for heterologous reactions
Erg10_reaction = Reaction(id="MVA1", name="2.0 Acetyl-CoA --> 1.0 Acetoacetyl-CoA + 1.0 CoA", upper_bound=1000.0)
Erg13_reaction = Reaction(id="MVA2", name="1.0 Acetyl-CoA + 1.0 Acetoacetyl-CoA --> 1.0 HMG-CoA + 1.0 CoA", upper_bound=1000.0)
tHMGR_reaction = Reaction(id="MVA3", name="1.0 HMG-CoA + 2.0 NADPH + 2 h_p --> 1.0 Mevalonate + 2.0 NADP + 1.0 CoA", upper_bound=1000.0)

Erg12_reaction = Reaction(id="MVA4", name="1.0 Mevalonate + 1.0 ATP --> 1.0 Mevalonat-P + 1.0 ADP + 1.0 h_p", upper_bound=1000.0)
Erg8_reaction = Reaction(id="MVA5", name="1.0 Mevalonate-P + 1.0 ATP --> 1.0 Mevalonat-PP + 1.0 ADP", upper_bound=1000.0)
Erg19_reaction = Reaction(id="MVA6", name="1.0 Mevalonate-PP + 1.0 ATP --> 1.0 IPP + 1.0 ADP + 1.0 PO4 + 1.0 CO2", lower_bound=None, upper_bound=1000.0) #irreversible reactions, so lower_bound = None
Idi1_reaction = Reaction(id="MVA7", name="1.0 IPP <--> 1.0 DMAPP", upper_bound=1000.0)

mFPS144_reaction = Reaction(id="MVA8", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0)
APS_reaction = Reaction(id="MVA9", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0)

# Add metabolites and reaction stoichiometry to reaction objects
Erg10_reaction.add_metabolites({
    Acetyl_CoA: -2.0,
    Acetoacetyl_CoA: 1.0,
    CoA: 1.0
})

Erg13_reaction.add_metabolites({
    Acetyl_CoA: -1.0,
    Acetoacetyl_CoA: -1.0,
    HMG_CoA: 1.0,
    CoA: 1.0
})

tHMGR_reaction.add_metabolites({
    HMG_CoA: -1.0,
    NADPH: -2.0,
    H_p: -2.0,
    Mevalonate: 1.0,
    NADP: 2.0,
    CoA: 1.0
})

Erg12_reaction.add_metabolites({
    Mevalonate: -1.0,
    ATP: -1.0,
    Mevalonate_P: 1.0,
    ADP: 1.0,
    H_p: 1.0
})

Erg8_reaction.add_metabolites({
    Mevalonate_P: -1.0,
    ATP: -1.0,
    Mevalonate_PP: 1.0,
    ADP: 1.0
})

Erg19_reaction.add_metabolites({
    Mevalonate_PP: -1.0,
    ATP: -1.0,
    IPP: 1.0,
    ADP: 1.0,
    Phosphate: 1.0,
    CO2: 1.0
})

Idi1_reaction.add_metabolites({
    IPP: -1.0,
    DMAPP: 1.0
}, reversibly=True)

mFPS144_reaction.add_metabolites({
    IPP: -1.0,
    DMAPP: -1.0,
    GPP: 1.0,
    Diphosphate: 1.0
})

APS_reaction.add_metabolites({
    GPP: -1.0,
    Alpha_pinene: 1.0,
    Diphosphate: 1.0
})

# Add reactions to model
reactionlist = [
    Erg10_reaction,
    Erg13_reaction,
    tHMGR_reaction,
    Erg12_reaction,
    Erg8_reaction,
    Erg19_reaction,
    Idi1_reaction,
    mFPS144_reaction,
    APS_reaction
]

for reaction in reactionlist:
    model.add_reaction(reaction)

