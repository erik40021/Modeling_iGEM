from Media.medium_analysis import NUTRIENTS_SC, medium_objectivevalue_xlsx, run_medium_test
from Utils.yeastGEM_io import read_yeast_model, write_yeast_model
from Utils.model_utils import plot_fluxes, summarize_properties
from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np



# constants
OVEREXPRESSION_LOWER_BOUND = 0.0 
KNOCK_DOWN_HIGHER_BOUND = 500
APINENE_OBJECTIVE_COEFFICIENT = 1.0
GROWTH_OBJECTIVE_COEFFICIENT = 1 - APINENE_OBJECTIVE_COEFFICIENT

def run_media_analysis(filename):
    model = simple_run()
    model.objective = {model.reactions.get_by_id("r_apinene_con"): APINENE_OBJECTIVE_COEFFICIENT, 
        model.reactions.get_by_id("r_2111"): GROWTH_OBJECTIVE_COEFFICIENT}
    exchange_reactions = model.exchanges #list_import_reactions_SC(model)
    solutions = run_medium_test(model, exchange_reactions, NUTRIENTS_SC) # run analysis
    medium_objectivevalue_xlsx(solutions, filename) #store as excel file


def simple_run():
    # builds and returns model with heterologuous reactions
    model = read_yeast_model()
    Acetyl_CoA = model.metabolites.get_by_id("s_0378")
    Acetoacetyl_CoA = Metabolite(id="s_4270", formula="C25H36N7O18P3S", name="acetoacetyl-CoA[p]", charge=-4, compartment="p") #the original model's last metabolite has the id s_4269
    CoA = model.metabolites.get_by_id("s_0534")
    HMG_CoA = Metabolite(id="s_4271", formula="C27H39N7O20P3S", name="3-hydroxy-3-methylglutaryl-CoA[p]", charge=-5, compartment="p")
    NADPH = model.metabolites.get_by_id("s_1215")
    NADP = model.metabolites.get_by_id("s_1211")
    H_p = model.metabolites.get_by_id("s_0801")
    Mevalonate = Metabolite(id="s_4272", formula="C6H11O4", name="(R)-mevalonate[p]", charge=-1, compartment="p")
    Mevalonate_P = Metabolite(id="s_4273", formula="C6H10O7P", name="(R)-5-phosphomevalonic acid[p]", charge=-3, compartment="p")
    Mevalonate_PP = Metabolite(id="s_4274", formula="C6H10O10P2", name="(R)-5-diphosphomevalonic acid[p]", charge=-4, compartment="p")
    ATP = model.metabolites.get_by_id("s_0439")
    ADP = model.metabolites.get_by_id("s_0399")
    IPP = Metabolite(id="s_4275", formula="C5H9O7P2", name="isopentyl diphosphate[p]", charge=-3, compartment="p")
    DMAPP = Metabolite(id="s_4276", formula="C5H9O7P2", name="dimethylallyl diphosphate[p]", charge=-3, compartment="p")
    GPP = Metabolite(id="s_4277", formula="C10H17O7P2", name="geranyl diphosphate[p]", charge=-3, compartment="p")
    Alpha_pinene = Metabolite(id="s_4278", formula="C10H16", name="(+)-alpha-pinene[p]", charge=0, compartment="p")
    Phosphate = Metabolite(id="s_4279", formula="HO4P", name="phosphate[p]", charge=-2, compartment="p")
    FPP = Metabolite(id="s_4280", formula="C15H25O7P2", name="farnesyl diphosphate[p]", charge=-3, compartment="p")
    NPP = Metabolite(id="s_4281",formula='C10H17O7P2', name='neryl diphosphate', charge=-3, compartment='p')
    CO2 = model.metabolites.get_by_id("s_0462")
    Diphosphate = model.metabolites.get_by_id("s_0638")

    # initialize reaction obchange_info_offsetects for heterologous reactions
    Erg10_reaction = Reaction(id="MVA1", name="2.0 Acetyl-CoA --> 1.0 Acetoacetyl-CoA + 1.0 CoA", upper_bound=1000.0, lower_bound=0)
    Erg13_reaction = Reaction(id="MVA2", name="1.0 Acetyl-CoA + 1.0 Acetoacetyl-CoA --> 1.0 HMG-CoA + 1.0 CoA", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    tHMGR_reaction = Reaction(id="MVA3", name="1.0 HMG-CoA + 2.0 NADPH + 2 h_p --> 1.0 Mevalonate + 2.0 NADP + 1.0 CoA", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    Erg12_reaction = Reaction(id="MVA4", name="1.0 Mevalonate + 1.0 ATP --> 1.0 Mevalonat-P + 1.0 ADP + 1.0 h_p", upper_bound=1000.0)
    Erg8_reaction = Reaction(id="MVA5", name="1.0 Mevalonate-P + 1.0 ATP --> 1.0 Mevalonat-PP + 1.0 ADP", upper_bound=1000.0)
    Erg19_reaction = Reaction(id="MVA6", name="1.0 Mevalonate-PP + 1.0 ATP --> 1.0 IPP + 1.0 ADP + 1.0 PO4 + 1.0 CO2", upper_bound=1000.0, lower_bound=0)
    Idi1_reaction = Reaction(id="MVA7", name="1.0 IPP <--> 1.0 DMAPP", upper_bound=1000.0, lower_bound=-1000) #reversible reaction
    AgGPPS2_reaction = Reaction(id="MVA8", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    mFPS144_GPP_reaction = Reaction(id="MVA9", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    mFPS144_FPP_reaction = Reaction(id="MVA10", name="1.0 IPP + 1.0 DMAPP --> 1.0 FPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    SiNPPS1_reaction = Reaction(id="MVA11", name="1.0 IPP + 1.0 DMAPP --> 1.0 NPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    APS_npp_reaction = Reaction(id="MVA12", name="1.0 NPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    APS_gpp_reaction = Reaction(id="MVA13", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
    APinene_con_reaction = Reaction(id="r_apinene_con", name="1.0 Alpha-Pinene ->", upper_bound=1000.0, lower_bound=0)

    # Add metabolites and reaction stoichiometry to reaction obchange_info_offsetects
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
    AgGPPS2_reaction.add_metabolites({
        IPP: -1.0,
        DMAPP: -1.0,
        GPP: 1.0,
        Diphosphate: 1.0
    })
    mFPS144_GPP_reaction.add_metabolites({
        IPP: -1.0,
        DMAPP: -1.0,
        GPP: 1.0,
        Diphosphate: 1.0
    })
    mFPS144_FPP_reaction.add_metabolites({
        GPP: -1.0,
        IPP: -1.0,
        FPP: 1.0,
        Diphosphate: 1.0
    })
    APS_gpp_reaction.add_metabolites({
        GPP: -1.0,
        Alpha_pinene: 1.0,
        Diphosphate: 1.0
    })
    APS_npp_reaction.add_metabolites({
        NPP: -1.0,
        Alpha_pinene: 1.0,
        Diphosphate: 1.0
    })
    SiNPPS1_reaction.add_metabolites({
        IPP: -1.0,
        DMAPP: -1.0,
        NPP: 1.0,
        Diphosphate: 1.0
    })
    APinene_con_reaction.add_metabolites({
        Alpha_pinene: -1.0
    })
    reactionlist = [Erg10_reaction,Erg13_reaction,tHMGR_reaction,Erg12_reaction,Erg8_reaction,Erg19_reaction,Idi1_reaction,
        mFPS144_GPP_reaction,mFPS144_FPP_reaction,AgGPPS2_reaction,SiNPPS1_reaction,APS_gpp_reaction,APS_npp_reaction,APinene_con_reaction]
    # Add reactions to model
    for reaction in reactionlist:
        model.add_reaction(reaction)
    return model



# ----- Erik's (old) playground (visualisation) -----

# visualisation constants:
FLUX_LABELS = ["biomass", "fpp_prod", "aps_gpp_apinene", "aggpps2_gpp", "npp_prod", "aps_npp_apinene", "mfps144_gpp"]
COLORS = ['seagreen', 'yellow', 'red', 'cornflowerblue', 'sandybrown', 'darkred', 'lightsteelblue']
STACK_INDICES = [(5,6), (2,3)]
CHANGE_INFOS = ['+Erg10', '++Erg13', '++tHMGR', '+Erg12', '+Erg8', '+Erg19', '+Idi1', '++mFPS144', '++AgGPPS2',
    '++SiNPPS1', '++APS & OF']

def simulate_all_changes():
    results = []
    change_info_offset = 0
    for i in range(0, len(reactionlist)-3):
        model.add_reaction(reactionlist[i])
        if i == 7: # special case: no solving/info reasonable if 'secondary' reaction of enzyme is added (because the manipulation's impact comprises all its reactions at once in reality)
            change_info_offset += 1
            continue
        results.append(solve_and_get_reaction_fluxes(model, change_biological=CHANGE_INFOS[i-change_info_offset]))
    # add aps reactions plus consumption manually and adapt objective function (!)
    model.add_reactions([reactionlist[11], reactionlist[12], reactionlist[13]])
    model.objective = {model.reactions.get_by_id("r_apinene_con"): APINENE_OBJECTIVE_COEFFICIENT, 
                    model.reactions.get_by_id("r_2111"): GROWTH_OBJECTIVE_COEFFICIENT}
    results.append(solve_and_get_reaction_fluxes(model, change_biological=CHANGE_INFOS[10]))
    plot_fluxes(results, FLUX_LABELS, COLORS, STACK_INDICES, 0.1, (10,7), save_as='scp_all.png')

def solve_and_get_reaction_fluxes(model, change_biological="", change_technical=""):
    solution = model.optimize()
    print("\n--- Change made to model: " + change_biological + "/" + change_technical + " ---")
    l = [solution.fluxes.get('r_2111'), solution.fluxes.get('MVA10'), solution.fluxes.get('MVA13'), solution.fluxes.get('MVA9'),
        solution.fluxes.get('MVA11'), solution.fluxes.get('MVA12'), solution.fluxes.get('MVA9')]
    for i in range(0,len(l)):
        if l[i] == None:
            l[i] = 0.0
    ap_prod = l[2] + l[5]
    gpp_prod = l[3] + l[6]
    print("aPinene production flux: " + str(ap_prod) + " of which")
    print(" - "+str(l[2])+" by aps with gpp"); print(" - "+str(l[5])+" by aps with npp")
    print("fpp production flux: " + str(l[1]))
    print("biomass flux: " + str(l[0])) #biomass
    print("gpp production flux: "+ str(gpp_prod) + " of which")
    print(" - "+str(l[3])+" by AgGPPS2"); print(" - "+str(l[6])+" by mFPS144")
    print("npp production flux: " + str(l[4]))
    print("-> SOLVER STATUS: "+ solution.status + " <-")
    return l + [change_biological]


# ----- Lasse's (even older) playground: -----

# model = read_yeast_model()
# reference = read_yeast_model() # brauchen wir das noch?

# # initilize metabolite obchange_info_offsetects for heterologous reactions in peroxisomes
# Acetyl_CoA = reference.metabolites.get_by_id("s_0378")
# Acetoacetyl_CoA = Metabolite(id="s_4270", formula="C25H36N7O18P3S", name="acetoacetyl-CoA[p]", charge=-4, compartment="p") #the original model's last metabolite has the id s_4269
# CoA = reference.metabolites.get_by_id("s_0534")
# HMG_CoA = Metabolite(id="s_4271", formula="C27H39N7O20P3S", name="3-hydroxy-3-methylglutaryl-CoA[p]", charge=-5, compartment="p")
# NADPH = reference.metabolites.get_by_id("s_1215")
# NADP = reference.metabolites.get_by_id("s_1211")
# H_p = reference.metabolites.get_by_id("s_0801")
# Mevalonate = Metabolite(id="s_4272", formula="C6H11O4", name="(R)-mevalonate[p]", charge=-1, compartment="p")
# Mevalonate_P = Metabolite(id="s_4273", formula="C6H10O7P", name="(R)-5-phosphomevalonic acid[p]", charge=-3, compartment="p")
# Mevalonate_PP = Metabolite(id="s_4274", formula="C6H10O10P2", name="(R)-5-diphosphomevalonic acid[p]", charge=-4, compartment="p")
# ATP = reference.metabolites.get_by_id("s_0439")
# ADP = reference.metabolites.get_by_id("s_0399")
# IPP = Metabolite(id="s_4275", formula="C5H9O7P2", name="isopentyl diphosphate[p]", charge=-3, compartment="p")
# DMAPP = Metabolite(id="s_4276", formula="C5H9O7P2", name="dimethylallyl diphosphate[p]", charge=-3, compartment="p")
# GPP = Metabolite(id="s_4277", formula="C10H17O7P2", name="geranyl diphosphate[p]", charge=-3, compartment="p")
# Alpha_pinene = Metabolite(id="s_4278", formula="C6H16", name="(+)-alpha-pinene[p]", charge=0, compartment="p")
# Phosphate = Metabolite(id="s_4279", formula="HO4P", name="phosphate[p]", charge=-2, compartment="p")
# FPP = Metabolite(id="s_4280", formula="C15H25O7P2", name="farnesyl diphosphate[p]", charge=-3, compartment="p")
# NPP = Metabolite(id="s_4281",formula='C10H17O7P2', name='neryl diphosphate', charge=-3, compartment='p')
# CO2 = reference.metabolites.get_by_id("s_0462")
# Diphosphate = reference.metabolites.get_by_id("s_0638")

# # initilize reaction obchange_info_offsetects for heterologous reactions
# Erg10_reaction = Reaction(id="MVA1", name="2.0 Acetyl-CoA --> 1.0 Acetoacetyl-CoA + 1.0 CoA", upper_bound=1000.0, lower_bound=0)

# Erg13_reaction = Reaction(id="MVA2", name="1.0 Acetyl-CoA + 1.0 Acetoacetyl-CoA --> 1.0 HMG-CoA + 1.0 CoA", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
# tHMGR_reaction = Reaction(id="MVA3", name="1.0 HMG-CoA + 2.0 NADPH + 2 h_p --> 1.0 Mevalonate + 2.0 NADP + 1.0 CoA", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)

# Erg12_reaction = Reaction(id="MVA4", name="1.0 Mevalonate + 1.0 ATP --> 1.0 Mevalonat-P + 1.0 ADP + 1.0 h_p", upper_bound=1000.0)
# Erg8_reaction = Reaction(id="MVA5", name="1.0 Mevalonate-P + 1.0 ATP --> 1.0 Mevalonat-PP + 1.0 ADP", upper_bound=1000.0)
# Erg19_reaction = Reaction(id="MVA6", name="1.0 Mevalonate-PP + 1.0 ATP --> 1.0 IPP + 1.0 ADP + 1.0 PO4 + 1.0 CO2", upper_bound=1000.0, lower_bound=0)
# Idi1_reaction = Reaction(id="MVA7", name="1.0 IPP <--> 1.0 DMAPP", upper_bound=1000.0, lower_bound=-1000) #reversible reaction

# AgGPPS2_reaction = Reaction(id="MVA8", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
# mFPS144_GPP_reaction = Reaction(id="MVA9", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
# mFPS144_FPP_reaction = Reaction(id="MVA10", name="1.0 IPP + 1.0 DMAPP --> 1.0 FPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
# SiNPPS1_reaction = Reaction(id="MVA11", name="1.0 IPP + 1.0 DMAPP --> 1.0 NPP + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
# APS_npp_reaction = Reaction(id="MVA12", name="1.0 NPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
# APS_gpp_reaction = Reaction(id="MVA13", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0, lower_bound=OVEREXPRESSION_LOWER_BOUND)
# APinene_con_reaction = Reaction(id="r_apinene_con", name="1.0 Alpha-Pinene ->", upper_bound=1000.0, lower_bound=0)


# # Add metabolites and reaction stoichiometry to reaction obchange_info_offsetects
# Erg10_reaction.add_metabolites({
#     Acetyl_CoA: -2.0,
#     Acetoacetyl_CoA: 1.0,
#     CoA: 1.0
# })

# Erg13_reaction.add_metabolites({
#     Acetyl_CoA: -1.0,
#     Acetoacetyl_CoA: -1.0,
#     HMG_CoA: 1.0,
#     CoA: 1.0
# })

# tHMGR_reaction.add_metabolites({
#     HMG_CoA: -1.0,
#     NADPH: -2.0,
#     H_p: -2.0,
#     Mevalonate: 1.0,
#     NADP: 2.0,
#     CoA: 1.0
# })

# Erg12_reaction.add_metabolites({
#     Mevalonate: -1.0,
#     ATP: -1.0,
#     Mevalonate_P: 1.0,
#     ADP: 1.0,
#     H_p: 1.0
# })

# Erg8_reaction.add_metabolites({
#     Mevalonate_P: -1.0,
#     ATP: -1.0,
#     Mevalonate_PP: 1.0,
#     ADP: 1.0
# })

# Erg19_reaction.add_metabolites({
#     Mevalonate_PP: -1.0,
#     ATP: -1.0,
#     IPP: 1.0,
#     ADP: 1.0,
#     Phosphate: 1.0,
#     CO2: 1.0
# })

# Idi1_reaction.add_metabolites({
#     IPP: -1.0,
#     DMAPP: 1.0
# }, reversibly=True)

# AgGPPS2_reaction.add_metabolites({
#     IPP: -1.0,
#     DMAPP: -1.0,
#     GPP: 1.0,
#     Diphosphate: 1.0
# })

# mFPS144_GPP_reaction.add_metabolites({
#     IPP: -1.0,
#     DMAPP: -1.0,
#     GPP: 1.0,
#     Diphosphate: 1.0
# })

# mFPS144_FPP_reaction.add_metabolites({
#     GPP: -1.0,
#     IPP: -1.0,
#     FPP: 1.0,
#     Diphosphate: 1.0
# })

# APS_gpp_reaction.add_metabolites({
#     GPP: -1.0,
#     Alpha_pinene: 1.0,
#     Diphosphate: 1.0
# })

# APS_npp_reaction.add_metabolites({
#     NPP: -1.0,
#     Alpha_pinene: 1.0,
#     Diphosphate: 1.0
# })

# SiNPPS1_reaction.add_metabolites({
#     IPP: -1.0,
#     DMAPP: -1.0,
#     NPP: 1.0,
#     Diphosphate: 1.0
# })

# APinene_con_reaction.add_metabolites({
#     Alpha_pinene: -1.0
# })

# # Add reactions to model
# reactionlist = [
#     Erg10_reaction,
#     Erg13_reaction,
#     tHMGR_reaction,
#     Erg12_reaction,
#     Erg8_reaction,
#     Erg19_reaction,
#     Idi1_reaction,
#     mFPS144_GPP_reaction, #7
#     mFPS144_FPP_reaction,
#     AgGPPS2_reaction,
#     SiNPPS1_reaction,
#     APS_gpp_reaction,
#     APS_npp_reaction,
#     APinene_con_reaction
# ]

# #write_yeast_model(model)

# # --- Step 2: FBA of both models ---
# referenceState = reference.optimize()
# engineeredState = model.optimize()

# # --- Step 3: Compare fluxes of reference and model with heterologous reactions ---
# df = pd.concat([referenceState.fluxes, engineeredState.fluxes], axis=1)
# df.columns = ["Fluxes_ref", "Fluxes_engineered"]

# #calculate change in Flux and clean data
# df["Fold_change"] = df["Fluxes_ref"] / df["Fluxes_engineered"]
# df["Fold_change"] = df["Fold_change"].replace([np.inf, -np.inf], np.nan)
# df["Fold_change"] = df["Fold_change"].fillna(1)
# df["Fold_change"] = df["Fold_change"].apply(lambda x: round(x, 2))
# print(df[df["Fold_change"] != 1])

# print("max: {}".format(np.max(df["Fold_change"])))
# print("min: {}".format(np.min(df["Fold_change"])))
# print("Recations with changes in flux: {}".format(len(df[(df["Fold_change"] != 1)])))

# plt.hist(df["Fold_change"], bins=5)
# plt.show()