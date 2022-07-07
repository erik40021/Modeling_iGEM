from Data.yeastGEM_io import read_yeast_model, write_yeast_model
from Model.model_utils import summarize_properties
from cobra import Model, Reaction, Metabolite
import optlang
import xlsxwriter
from matplotlib import pyplot as plt
import numpy as np

#constants:
OVEREXPRESSION_LOWER_BOUND = 0.2
KNOCK_DOWN_HIGHER_BOUND = 500
APINENE_OBJECTIVE_COEFFICIENT = 0.2
GROWTH_OBJECTIVE_COEFFICIENT = 0.8

def simulate_variant_1():
    ''' simulates first sc_cyto variant with two overexpressions (Erg13 and tHMGR) and one knock-in (APS) '''
    model, results = build_basic_model()
    # add manipulations one by one and print reactions fluxes every time
    model.reactions.get_by_id("r_0559").lower_bound = OVEREXPRESSION_LOWER_BOUND #erg13
    results.append(solve_and_get_reaction_fluxes(model, "Erg13 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_0558").lower_bound = OVEREXPRESSION_LOWER_BOUND #tHMGR
    results.append(solve_and_get_reaction_fluxes(model, "tHMGR ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "APS ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    #write_to_excel(results)
    plot_fluxes(results, 0.1, (10,7))
    return results

def simulate_variant_2():
    ''' simulates second sc_cyto variant with three knock-ins (APS, AgGPPS2, mFPS144) and one knock-out (Erg20)'''
    model = build_basic_model()
    # add manipulations one by one and print reactions fluxes every time
    add_erg20_ko_and_alternative(model)
    solve_and_get_reaction_fluxes(model, "knocked out Erg20 and replaced with mFPS144 and AgGPPs")
    # add overexpressions of new genes
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed APS/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_dmapp_aggpps2_gpp").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed AgGPPS2/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_ipp_dmapp_mfps144_gpp").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed mFPS144/set lower bound of IPP + DMAPP -> GPP to 0.2")
    model.reactions.get_by_id("r_ipp_gpp_mfps144_fpp").higher_bound = 500
    solve_and_get_reaction_fluxes(model, "implement low affinity of mFPS144 for IPP + GPP -> FPP/set higher bound to 500")


def simulate_variant_3():
    ''' simulates third sc_cyto variant with two knock-ins (APS, SiNPPS2) and cupper-dependant promoter pCTR3 for
        Erg20 regulation '''
    model = build_basic_model()
    # add manipulations one by one and print reactions fluxes every time
    add_npp_pathway(model)
    solve_and_get_reaction_fluxes(model, "added NPP pathway")
    # add overexpressions of new genes
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = 0.2
    model.reactions.get_by_id("r_npp_aps_apinene").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed APS/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_dmapp_sinpps2_npp").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed SiNPPS2/set lower bound of corresponding reaction(s) to 0.2")
    # add downregulation of Erg20
    model.reactions.get_by_id("r_0355").higher_bound = 500 # gpp production
    model.reactions.get_by_id("r_0462").higher_bound = 500 # fpp production
    solve_and_get_reaction_fluxes(model, "implement downregulation of Erg20 by pCTR3/set higher bound of reactions to 500")


def simulate_fpp_knockout():
    model = build_basic_model()
    model.reactions.r_0462.knock_out()
    solve_and_get_reaction_fluxes(model, "knocked out fpp production reaction")

def simulate_all_changes():
    model = build_basic_model()
    model.reactions.get_by_id("r_0559").lower_bound = 0.2 #erg13
    solve_and_get_reaction_fluxes(model, "overexpressed Erg13/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_0558").lower_bound = 0.2 #tHMGR
    solve_and_get_reaction_fluxes(model, "overexpressed tHMGR/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed APS/set lower bound of corresponding reaction(s) to 0.2")
    add_erg20_ko_and_alternative(model)
    solve_and_get_reaction_fluxes(model, "knocked out Erg20 and replaced with mFPS144 and AgGPPs")
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed APS/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_dmapp_aggpps2_gpp").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed AgGPPS2/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_ipp_dmapp_mfps144_gpp").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed mFPS144/set lower bound of IPP + DMAPP -> GPP to 0.2")
    model.reactions.get_by_id("r_ipp_gpp_mfps144_fpp").higher_bound = 500
    solve_and_get_reaction_fluxes(model, "implement low affinity of mFPS144 for IPP + GPP -> FPP/set higher bound to 500")
    add_npp_pathway(model)
    solve_and_get_reaction_fluxes(model, "added NPP pathway")
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = 0.2
    model.reactions.get_by_id("r_npp_aps_apinene").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed APS/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_dmapp_sinpps2_npp").lower_bound = 0.2
    solve_and_get_reaction_fluxes(model, "overexpressed SiNPPS2/set lower bound of corresponding reaction(s) to 0.2")
    model.reactions.get_by_id("r_0355").higher_bound = 500 # gpp production
    model.reactions.get_by_id("r_0462").higher_bound = 500 # fpp production
    solve_and_get_reaction_fluxes(model, "implement downregulation of Erg20 by pCTR3/set higher bound of reactions to 500")


def build_basic_model():
    results = []
    model = read_yeast_model()
    results.append(solve_and_get_reaction_fluxes(model, "before"))
    reactions_list = list()
    # new metabolites:
    aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')
    # biological reaction: 1.0 GPP -> 1.0 aPinene
    react_gpp_alphapinene = Reaction('r_gpp_aps_apinene')
    react_gpp_alphapinene.name = 'GPP -> aPinene'
    react_gpp_alphapinene.subsystem = 'Cytoplasm'
    react_gpp_alphapinene.lower_bound = 0.  # This is the default
    react_gpp_alphapinene.upper_bound = 1000.  # This is the default
    react_gpp_alphapinene.add_metabolites({
        model.metabolites.get_by_id("s_0745"): -1.0, # GPP (C10H17O7P2)
        model.metabolites.get_by_id("s_0633"): 1.0, # diphosphate/ppi (HO7P2)
        aPinene: 1.0
    })
    react_gpp_alphapinene.gene_reaction_rule = '(APS)'
    
    # reaction that consumes aPinene
    react_alphapinen_con = Reaction('r_apinene_con')
    react_alphapinen_con.name = 'aPinene -> '
    react_alphapinen_con.subsystem = 'Cytoplasm'
    react_alphapinen_con.lower_bound = 0.  # This is the default
    react_alphapinen_con.upper_bound = 1000.  # This is the default
    react_alphapinen_con.add_metabolites({
        aPinene: -1.0
    })
    reactions_list.append(react_gpp_alphapinene)
    reactions_list.append(react_alphapinen_con)
    model.add_reactions(reactions_list)

    model.objective = {model.reactions.get_by_id("r_apinene_con"): APINENE_OBJECTIVE_COEFFICIENT, 
                    model.reactions.get_by_id("r_2111"): GROWTH_OBJECTIVE_COEFFICIENT}
    results.append(solve_and_get_reaction_fluxes(model, "GPP -> aPinene", "and objective to aPinene:"
                    +str(APINENE_OBJECTIVE_COEFFICIENT) + "/growth:" + str(GROWTH_OBJECTIVE_COEFFICIENT)))
    return model, results


def add_erg20_ko_and_alternative(model):
    reactions_list = list()
    react_dmapp_gpp = Reaction('r_dmapp_aggpps2_gpp')
    react_dmapp_gpp.name = 'DMAPP -> GPP'
    react_dmapp_gpp.subsystem = 'Cytoplasm'
    react_dmapp_gpp.lower_bound = 0.  # This is the default
    react_dmapp_gpp.upper_bound = 1000.  # This is the default
    react_dmapp_gpp.add_metabolites({
        model.metabolites.get_by_id("s_1376"): -2.0, # DMAPP (C5H9O7P2)
        model.metabolites.get_by_id("s_0745"): 1.0, # GPP (C10H17O7P2)
        model.metabolites.get_by_id("s_0633"): 1.0 # diphosphate/ppi (HO7P2)
    })
    react_dmapp_gpp.gene_reaction_rule = '(AgGPPS2)'

    react_ipp_dmapp_gpp = Reaction('r_ipp_dmapp_mfps144_gpp')
    react_ipp_dmapp_gpp.name = 'IPP + DMAPP -> GPP'
    react_ipp_dmapp_gpp.subsystem = 'Cytoplasm'
    react_ipp_dmapp_gpp.lower_bound = 0.  # This is the default
    react_ipp_dmapp_gpp.upper_bound = 1000.  # This is the default
    react_ipp_dmapp_gpp.add_metabolites({
        model.metabolites.get_by_id("s_0934"): -1.0, # IPP (C5H9O7P2)
        model.metabolites.get_by_id("s_1376"): -1.0, # DMAPP (C5H9O7P2)
        model.metabolites.get_by_id("s_0745"): 1.0, # GPP (C10H17O7P2)
        model.metabolites.get_by_id("s_0633"): 1.0 # diphosphate/ppi (HO7P2)
    })
    react_ipp_dmapp_gpp.gene_reaction_rule = '(mFPS144)'

    react_ipp_gpp_fpp = Reaction('r_ipp_gpp_mfps144_fpp')
    react_ipp_gpp_fpp.name = 'IPP + GPP -> FPP'
    react_ipp_gpp_fpp.subsystem = 'Cytoplasm'
    react_ipp_gpp_fpp.lower_bound = 0.  # This is the default
    react_ipp_gpp_fpp.upper_bound = 1000.  # This is the default
    react_ipp_gpp_fpp.add_metabolites({
        model.metabolites.get_by_id("s_0934"): -1.0, # IPP (C5H9O7P2)
        model.metabolites.get_by_id("s_0745"): -1.0, # GPP (C10H17O7P2)
        model.metabolites.get_by_id("s_0190"): 1.0, # FPP (C15H25O7P2)
        model.metabolites.get_by_id("s_0633"): 1.0 # diphosphate/ppi (HO7P2)
    }) # TODO: set upper bound to something > 1000 because of catalytic preference of mFPS144 for GPP
    react_ipp_gpp_fpp.gene_reaction_rule = '(mFPS144)'

    reactions_list.append(react_dmapp_gpp)
    reactions_list.append(react_ipp_dmapp_gpp)
    reactions_list.append(react_ipp_gpp_fpp)
    model.add_reactions(reactions_list)

    model.genes.YJL167W.knock_out() # knocks out Erg20


def add_npp_pathway(model):
    reactions_list = list()
    npp = Metabolite('NPP', formula='C10H17O7P2', name='Neryl pyrophosphate', charge=-3, compartment='c')
    react_dmapp_npp = Reaction('r_dmapp_sinpps2_npp')
    react_dmapp_npp.name = 'DMAPP -> NPP'
    react_dmapp_npp.subsystem = 'Cytoplasm'
    react_dmapp_npp.lower_bound = 0.  # This is the default
    react_dmapp_npp.upper_bound = 1000.  # This is the default
    react_dmapp_npp.add_metabolites({
        model.metabolites.get_by_id("s_0934"): -1.0, # IPP (C5H9O7P2)
        model.metabolites.get_by_id("s_1376"): -1.0, # DMAPP (C5H9O7P2)
        npp: 1.0,
        model.metabolites.get_by_id("s_0633"): 1.0 # diphosphate/ppi (HO7P2)
    })
    react_dmapp_npp.gene_reaction_rule = '(SiNPPS2)'

    react_npp_alphapinene = Reaction('r_npp_aps_apinene')
    react_npp_alphapinene.name = 'GPP -> aPinene'
    react_npp_alphapinene.subsystem = 'Cytoplasm'
    react_npp_alphapinene.lower_bound = 0.  # This is the default
    react_npp_alphapinene.upper_bound = 1000.  # This is the default
    react_npp_alphapinene.add_metabolites({
        npp: -1.0, # (C10H17O7P2)
        model.metabolites.get_by_id("aPinene"): 1.0, # alphapinene (C10H16)
        model.metabolites.get_by_id("s_0633"): 1.0 # diphosphate/ppi (HO7P2)
    })
    react_npp_alphapinene.gene_reaction_rule = '(APS)'

    reactions_list.append(react_dmapp_npp)
    reactions_list.append(react_npp_alphapinene)
    model.add_reactions(reactions_list)
    
    
def solve_and_get_reaction_fluxes(model, change_biological="", change_technical=""):
    solution = model.optimize()
    print("\n--- Change made to model: " + change_biological + "/" + change_technical + " ---")
    aps_gpp_apinene = solution.fluxes.get('r_gpp_aps_apinene')
    aps_npp_apinene = solution.fluxes.get('r_npp_aps_apinene')
    erg20_gpp = solution.fluxes.get('r_0355')
    aggpps2_gpp = solution.fluxes.get('r_dmapp_aggpps2_gpp')
    mfps144_gpp = solution.fluxes.get('r_ipp_dmapp_mfps144_gpp')
    npp_prod = solution.fluxes.get("r_dmapp_sinpps2_npp")
    fpp_prod = solution.fluxes.get('r_0462')
    biomass_prod = solution.fluxes.get('r_2111')
    l = [biomass_prod, aps_gpp_apinene, aps_npp_apinene, aggpps2_gpp, mfps144_gpp, erg20_gpp, fpp_prod, npp_prod]
    for i in range(0,len(l)):
        if l[i] == None:
            l[i] = 0.0
    ap_prod = l[1] + l[2]
    gpp_prod = l[3] + l[4] + l[5]
    if solution.fluxes.get('r_apinene_con') != ap_prod:
        print("-- Warning: aPinene consumption != production!")
    print("aPinene production flux: " + str(ap_prod) + " of which")
    print(" - "+str(aps_gpp_apinene)+" by aps with gpp"); print(" - "+str(aps_npp_apinene)+" by aps with npp")
    print("fpp production flux: " + str(fpp_prod))
    print("biomass flux: " + str(biomass_prod)) #biomass
    print("gpp production flux: "+ str(gpp_prod) + " of which")
    print(" - "+str(erg20_gpp)+" by Erg20"); print(" - "+str(aggpps2_gpp)+" by AgGPPS2")
    print(" - "+str(mfps144_gpp)+" by mFPS144")
    print("npp production flux: " + str(npp_prod))
    return l + [change_biological]
    

def plot_fluxes(res, bar_width=0.1, size=(10,7)):
    labels = ["biomass", "aps_gpp_apinene", "aps_npp_apinene", "aggpps2_gpp",
    "mfps144_gpp", "erg20_gpp", "fpp_prod", "npp_prod"]
    subtitles = [i[8] for i in res]
    x = np.arange(len(subtitles))  # the label locations
    width = bar_width # the width of the bars
    fig, ax = plt.subplots()
    offset = -3.5
    for i in range(0,7):
        ax.bar(x + (offset+i)*width, [j[i] for j in res], width, label=labels[i])
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Flux [mmol/gdcw/h]')
    ax.set_title('Fluxes after Manipulations')
    ax.set_xticks(x, subtitles)
    ax.legend()
    fig.tight_layout()
    fig.set_size_inches(size[0], size[1])
    plt.show()


def write_to_excel(res):
    results = [["biomass", "aps_gpp_apinene", "aps_npp_apinene", "aggpps2_gpp",
    "mfps144_gpp", "erg20_gpp", "fpp_prod", "npp_prod", "change_info"]] + res
    with xlsxwriter.Workbook('Output/scc.xlsx') as workbook:
        worksheet = workbook.add_worksheet()
        for row_num, data in enumerate(results):
            worksheet.write_row(row_num, 0, data)


# ----------------- manual all-in-one built (old) -----------------------

# --- step 1: load scaffold model
model = read_yeast_model() # loads existing SBML yeastGEM model 
# check content of scaffold
# summarize_properties(model)



# --- step 2: adjust model (add new and adjust knocked-down/deleted or overexpressed 
#                       reactions and corresponding metabolites)

# How to add new reaction: Create object and add respective metabolites
# How to 'knock-out' an existing reaction: set upper and lower bound to 0 (= no flux can pass through)
# How to overexpress a reaction: ?

reactions_list = list()

# new metabolites:
aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')

# biological reaction: 1.0 GPP -> 1.0 aPinene
react_gpp_alphapinene = Reaction('r_gpp_aps_apinene')
react_gpp_alphapinene.name = 'GPP -> aPinene'
react_gpp_alphapinene.subsystem = 'Cytoplasm'
react_gpp_alphapinene.lower_bound = 0.  # This is the default
react_gpp_alphapinene.upper_bound = 1000.  # This is the default
react_gpp_alphapinene.add_metabolites({
    model.metabolites.get_by_id("s_0745"): -1.0, # GPP (C10H17O7P2)
    model.metabolites.get_by_id("s_0633"): 1.0, # diphosphate/ppi (HO7P2)
    aPinene: 1.0
})
react_gpp_alphapinene.gene_reaction_rule = '(APS)'
#TODO: add associated gene(s)!

# reaction that consumes aPinene
react_alphapinen_con = Reaction('r_apinene_con')
react_alphapinen_con.name = 'aPinene -> '
react_alphapinen_con.subsystem = 'Cytoplasm'
react_alphapinen_con.lower_bound = 0.  # This is the default
react_alphapinen_con.upper_bound = 1000.  # This is the default
react_alphapinen_con.add_metabolites({
    aPinene: -1.0
})

reactions_list.append(react_gpp_alphapinene)
reactions_list.append(react_alphapinen_con)

model.add_reactions(reactions_list)

# adjust lower bounds of reactions to simulate overexpression of underlying genes

model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = 0.2
model.reactions.get_by_id("r_0559").lower_bound = 0.2 #erg13 




# --- step 3: set model objective (reaction(s) whose fluxes shall be maximized)
# model.objective = "r_apinene_con"
# mehrere objectives definieren:
model.objective = {model.reactions.get_by_id("r_apinene_con"): 0.2, model.reactions.get_by_id("r_2111"): 0.8}

print(f"\nModel objective expression: {model.objective.expression}")
print(f"Model objective direction: {model.objective.direction}\n")




# --- step 4: run solver
# declare which solver to use (best one: Gurobipy)
model.solver = 'gurobi'

solution = model.optimize()
print(solution)
print(f"\nObjective value of solution: {solution.objective_value}")


# save solution object
#solution.fluxes.to_csv("out.csv", index=False)

# summarize solution
print(model.summary())


# last step: save model
#write_yeast_model(model)   # saving, writes SBML file


