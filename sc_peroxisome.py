from Data.yeastGEM_io import read_yeast_model, write_yeast_model
from Model.model_utils import summarize_properties
from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
# --- Step 1: Create reference model and model with heterologous reactions for alpha-pienen synthesis in peroxisomes ---

# create reference and engineered model
model = read_yeast_model()
reference = read_yeast_model()

# initilize metabolite objects for heterologous reactions in peroxisomes
Acetyl_CoA = reference.metabolites.get_by_id("s_0378")
Acetoacetyl_CoA = Metabolite(id="s_4270", formula="C25H36N7O18P3S", name="acetoacetyl-CoA[p]", charge=-4, compartment="p") #the original model's last metabolite has the id s_4269
CoA = reference.metabolites.get_by_id("s_0534")
HMG_CoA = Metabolite(id="s_4271", formula="C27H39N7O20P3S", name="3-hydroxy-3-methylglutaryl-CoA[p]", charge=-5, compartment="p")
NADPH = reference.metabolites.get_by_id("s_1215")
NADP = reference.metabolites.get_by_id("s_1211")
H_p = reference.metabolites.get_by_id("s_0801")
Mevalonate = Metabolite(id="s_4272", formula="C6H11O4", name="(R)-mevalonate[p]", charge=-1, compartment="p")
Mevalonate_P = Metabolite(id="s_4273", formula="C6H10O7P", name="(R)-5-phosphomevalonic acid[p]", charge=-3, compartment="p")
Mevalonate_PP = Metabolite(id="s_4274", formula="C6H10O10P2", name="(R)-5-diphosphomevalonic acid[p]", charge=-4, compartment="p")
ATP = reference.metabolites.get_by_id("s_0439")
ADP = reference.metabolites.get_by_id("s_0399")
IPP = Metabolite(id="s_4275", formula="C5H9O7P2", name="isopentyl diphosphate[p]", charge=-3, compartment="p")
DMAPP = Metabolite(id="s_4276", formula="C5H9O7P2", name="dimethylallyl diphosphate[p]", charge=-3, compartment="p")
GPP = Metabolite(id="s_4277", formula="C10H17O7P2", name="geranyl diphosphate[p]", charge=-3, compartment="p")
Alpha_pinene = Metabolite(id="s_4278", formula="C6H16", name="(+)-alpha-pinene[p]", charge=0, compartment="p")
Phosphate = Metabolite(id="s_4279", formula="HO4P", name="phosphate[p]", charge=-2, compartment="p")
CO2 = reference.metabolites.get_by_id("s_0462")
Diphosphate = reference.metabolites.get_by_id("s_0638")

# initilize reaction objects for heterologous reactions
Erg10_reaction = Reaction(id="MVA1", name="2.0 Acetyl-CoA --> 1.0 Acetoacetyl-CoA + 1.0 CoA", upper_bound=1000.0)
Erg13_reaction = Reaction(id="MVA2", name="1.0 Acetyl-CoA + 1.0 Acetoacetyl-CoA --> 1.0 HMG-CoA + 1.0 CoA", upper_bound=1000.0)
tHMGR_reaction = Reaction(id="MVA3", name="1.0 HMG-CoA + 2.0 NADPH + 2 h_p --> 1.0 Mevalonate + 2.0 NADP + 1.0 CoA", upper_bound=1000.0)

Erg12_reaction = Reaction(id="MVA4", name="1.0 Mevalonate + 1.0 ATP --> 1.0 Mevalonat-P + 1.0 ADP + 1.0 h_p", upper_bound=1000.0)
Erg8_reaction = Reaction(id="MVA5", name="1.0 Mevalonate-P + 1.0 ATP --> 1.0 Mevalonat-PP + 1.0 ADP", upper_bound=1000.0)
Erg19_reaction = Reaction(id="MVA6", name="1.0 Mevalonate-PP + 1.0 ATP --> 1.0 IPP + 1.0 ADP + 1.0 PO4 + 1.0 CO2", lower_bound=None, upper_bound=1000.0) #irreversible reactions, so lower_bound = None
Idi1_reaction = Reaction(id="MVA7", name="1.0 IPP <--> 1.0 DMAPP", upper_bound=1000.0)

AgGPPS2_GPP_reaction = Reaction(id="MVA8", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0)
AgGPPS2_FPP_reaction = Reaction(id="MVA9", name="", upper_bound=1000.0) #TODO
mFPS144_reaction = Reaction(id="MVA10", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0)
APS_npp_reaction = Reaction(id="MVA11", name="1.0 FPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0)
APS_gpp_reaction = Reaction(id="MVA12", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0)
#TODO add missing reactions (SiNPPS1)

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

APS_gpp_reaction.add_metabolites({
    GPP: -1.0,
    Alpha_pinene: 1.0,
    Diphosphate: 1.0
})
#TODO add missing metabolites to reactions

# Add reactions to model
reactionlist = [
    Erg10_reaction,
    Erg13_reaction,
    tHMGR_reaction,
    Erg12_reaction,
    Erg8_reaction,
    Erg19_reaction,
    Idi1_reaction,
    mFPS144_reaction, #7
    AgGPPS2_GPP_reaction,
    AgGPPS2_FPP_reaction,
    APS_gpp_reaction,
    APS_npp_reaction
]

for reaction in reactionlist:
    model.add_reaction(reaction)

def simulate_all_changes():
    #TODO
    return

def solve_and_get_reaction_fluxes(model, change_biological="", change_technical=""):
    solution = model.optimize()
    print("\n--- Change made to model: " + change_biological + "/" + change_technical + " ---")
    l = [solution.fluxes.get('r_2111')] + reactionlist
    for i in range(0,len(l)):
        if l[i] == None:
            l[i] = 0.0
    ap_prod = l[2] + l[5]
    gpp_prod = l[3] + l[6] + l[7]
    print("aPinene production flux: " + str(ap_prod) + " of which")
    print(" - "+str(l[2])+" by aps with gpp"); print(" - "+str(l[5])+" by aps with npp")
    print("fpp production flux: " + str(l[10]))
    print("biomass flux: " + str(l[0])) #biomass
    print("gpp production flux: "+ str(gpp_prod) + " of which")
    print(" - "+str(l[9])+" by AgGPPS2"); print(" - "+str(l[8])+" by mFPS144")
    print("npp production flux: " + str(l[4]))
    return l + [change_biological]

# --- Step 2: FBA of both models ---
referenceState = reference.optimize()
engineeredState = model.optimize()

# --- Step 3: Compare fluxes of reference and model with heterologous reactions ---
df = pd.concat([referenceState.fluxes, engineeredState.fluxes], axis=1)
df.columns = ["Fluxes_ref", "Fluxes_engineered"]

#calculate change in Flux and clean data
df["Fold_change"] = df["Fluxes_ref"] / df["Fluxes_engineered"]
df["Fold_change"] = df["Fold_change"].replace([np.inf, -np.inf], np.nan)
df["Fold_change"] = df["Fold_change"].fillna(1)
df["Fold_change"] = df["Fold_change"].apply(lambda x: round(x, 2))
print(df[df["Fold_change"] != 1])

print("max: {}".format(np.max(df["Fold_change"])))
print("min: {}".format(np.min(df["Fold_change"])))
print("Recations with changes in flux: {}".format(len(df[(df["Fold_change"] != 1)])))

plt.hist(df["Fold_change"], bins=5)
plt.show()


