import math
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
    plot_fluxes(results, 0.1, (10,7), save_as='scc_variant1.png')

def simulate_variant_2():
    ''' simulates second sc_cyto variant with three knock-ins (APS, AgGPPS2, mFPS144) and one knock-out (Erg20)'''
    model, results = build_basic_model()
    # add manipulations one by one and print reactions fluxes every time
    add_erg20_ko_and_alternative(model)
    results.append(solve_and_get_reaction_fluxes(model, "Erg20-KO +mFPS144 +AgGPPs"))
    # add overexpressions of new genes
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "APS ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_dmapp_aggpps2_gpp").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "AgGPPS2 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_ipp_dmapp_mfps144_gpp").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "mFPS144 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_ipp_gpp_mfps144_fpp").higher_bound = KNOCK_DOWN_HIGHER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "mFPS144 low FPP aff.", "higher bound=" + str(KNOCK_DOWN_HIGHER_BOUND)))
    plot_fluxes(results, 0.1, (12, 6), save_as='scc_variant2.png')


def simulate_variant_3():
    ''' simulates third sc_cyto variant with two knock-ins (APS, SiNPPS2) and cupper-dependant promoter pCTR3 for
        Erg20 regulation '''
    model, results = build_basic_model()
    # add manipulations one by one and print reactions fluxes every time
    add_npp_pathway(model)
    results.append(solve_and_get_reaction_fluxes(model, "+SiNPPS2"))
    # add overexpressions of new genes
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    model.reactions.get_by_id("r_npp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "APS ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_dmapp_sinpps2_npp").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "SiNPPS2 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    # add downregulation of Erg20
    model.reactions.get_by_id("r_0355").higher_bound = KNOCK_DOWN_HIGHER_BOUND # gpp production
    model.reactions.get_by_id("r_0462").higher_bound = KNOCK_DOWN_HIGHER_BOUND # fpp production
    results.append(solve_and_get_reaction_fluxes(model, "Erg20-KD by pCTR3", "higher bound=" + str(KNOCK_DOWN_HIGHER_BOUND)))
    plot_fluxes(results, 0.1, (12,6), save_as='scc_variant3.png')


def simulate_fpp_knockout():
    model = build_basic_model()
    model.reactions.r_0462.knock_out()
    solve_and_get_reaction_fluxes(model, "knocked out fpp production reaction")

def simulate_all_changes():
    model, results = build_basic_model()
    # part 1:
    model.reactions.get_by_id("r_0559").lower_bound = OVEREXPRESSION_LOWER_BOUND #erg13
    results.append(solve_and_get_reaction_fluxes(model, "Erg13 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_0558").lower_bound = OVEREXPRESSION_LOWER_BOUND #tHMGR
    results.append(solve_and_get_reaction_fluxes(model, "tHMGR ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "APS ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    # part 2:
    add_erg20_ko_and_alternative(model)
    results.append(solve_and_get_reaction_fluxes(model, "Erg20-KO +mFPS144 +AgGPPs"))
    # add overexpressions of new genes
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "APS ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_dmapp_aggpps2_gpp").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "AgGPPS2 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_ipp_dmapp_mfps144_gpp").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "mFPS144 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_ipp_gpp_mfps144_fpp").higher_bound = KNOCK_DOWN_HIGHER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "mFPS144 low FPP aff.", "higher bound=" + str(KNOCK_DOWN_HIGHER_BOUND)))
    # part 3:
    add_npp_pathway(model)
    results.append(solve_and_get_reaction_fluxes(model, "+SiNPPS2"))
    # add overexpressions of new genes
    model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    model.reactions.get_by_id("r_npp_aps_apinene").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "APS ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    model.reactions.get_by_id("r_dmapp_sinpps2_npp").lower_bound = OVEREXPRESSION_LOWER_BOUND
    results.append(solve_and_get_reaction_fluxes(model, "SiNPPS2 ovexp", "lower bound=" + str(OVEREXPRESSION_LOWER_BOUND)))
    # add downregulation of Erg20
    model.reactions.get_by_id("r_0355").higher_bound = KNOCK_DOWN_HIGHER_BOUND # gpp production
    model.reactions.get_by_id("r_0462").higher_bound = KNOCK_DOWN_HIGHER_BOUND # fpp production
    results.append(solve_and_get_reaction_fluxes(model, "Erg20-KD by pCTR3", "higher bound=" + str(KNOCK_DOWN_HIGHER_BOUND)))
    plot_fluxes(results, 0.1, (15,8), save_as='scc_all_manipulations.png')

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
    l = [solution.fluxes.get('r_2111'), solution.fluxes.get('r_0462'), solution.fluxes.get('r_gpp_aps_apinene'), 
        solution.fluxes.get('r_0355'), solution.fluxes.get("r_dmapp_sinpps2_npp"), solution.fluxes.get('r_npp_aps_apinene'),
        solution.fluxes.get('r_dmapp_aggpps2_gpp'), solution.fluxes.get('r_ipp_dmapp_mfps144_gpp'), solution.fluxes.get('r_apinene_con')]
    for i in range(0,len(l)):
        if l[i] == None:
            l[i] = 0.0
    ap_prod = l[2] + l[5]
    gpp_prod = l[3] + l[6] + l[7]
    if l[8] != ap_prod:
        print("-- Warning: aPinene consumption != production!")
    print("aPinene production flux: " + str(ap_prod) + " of which")
    print(" - "+str(l[2])+" by aps with gpp"); print(" - "+str(l[5])+" by aps with npp")
    print("fpp production flux: " + str(l[1]))
    print("biomass flux: " + str(l[0])) #biomass
    print("gpp production flux: "+ str(gpp_prod) + " of which")
    print(" - "+str(l[6])+" by AgGPPS2"); print(" - "+str(l[7])+" by mFPS144")
    print(" - "+str(l[3])+" by Erg20"); 
    print("npp production flux: " + str(l[4]))
    return l[:len(l)-1] + [change_biological]
    

def plot_fluxes(res, bar_width=0.1, size=(10,7), save_as=None):
    labels = ["biomass", "fpp_prod", "aps_gpp_apinene", "erg20_gpp", "npp_prod", "aps_npp_apinene", "aggpps2_gpp","mfps144_gpp"]
    colors = ['seagreen', 'yellow', 'red', 'royalblue', 'sandybrown', 'darkred', 'cornflowerblue', 'lightsteelblue']
    subtitles = [i[8] for i in res]
    ys = []
    for i in range(0, 8):
        ys += [j[i] for j in res]
    for i in range(0, len(res)):
        ys[len(res)*2 + i] += ys[len(res)*5 + i]
        ys[len(res)*3 + i] += ys[len(res)*6 + i] + ys[len(res)*7 + i]
    max_y = max(ys)
    for i in range(0,len(subtitles)):
        if len(subtitles[i]) > 14:
            a = str.split(subtitles[i])
            subtitles[i] = ' '.join(a[:math.ceil(len(a)/2)]) + "\n" + ' '.join(a[math.ceil(len(a)/2):])
    x = np.arange(len(subtitles))  # the label locations
    width = bar_width # the width of the bars
    fig, ax = plt.subplots()
    offset = -2
    ax.bar(x + (offset)*width, [j[0] for j in res], width, label=labels[0], color=colors[0])
    ax.bar(x + (offset+1)*width, [j[1] for j in res], width, label=labels[1], color=colors[1])
    ax.bar(x + (offset+2)*width, [j[2] for j in res], width, label=labels[2], color=colors[2])
    ax.bar(x + (offset+2)*width, [j[5] for j in res], width, bottom=[j[2] for j in res], label=labels[5], color=colors[5])
    ax.bar(x + (offset+3)*width, [j[3] for j in res], width, label=labels[3], color=colors[3])
    ax.bar(x + (offset+3)*width, [j[6] for j in res], width, bottom=[j[3] for j in res], label=labels[6], color=colors[6])
    ax.bar(x + (offset+3)*width, [j[7] for j in res], width, bottom=[j[6] for j in res], label=labels[7], color=colors[7])
    ax.bar(x + (offset+4)*width, [j[4] for j in res], width, label=labels[4], color=colors[4])
     # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Flux [mmol/gdcw/h]')
    ax.set_title('Fluxes after Manipulations')
    if len(res) > 5:
        ax.set_xticks(x, subtitles, rotation=45)
    else:
        ax.set_xticks(x, subtitles)
    ax.set_ylim(top=max_y + 0.05*max_y)
    ax.legend()
    fig.tight_layout()
    fig.set_size_inches(size[0], size[1])
    if save_as is not None:
        plt.savefig("Output/" + save_as, dpi=500)
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

# reactions_list = list()

# # new metabolites:
# aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')

# # biological reaction: 1.0 GPP -> 1.0 aPinene
# react_gpp_alphapinene = Reaction('r_gpp_aps_apinene')
# react_gpp_alphapinene.name = 'GPP -> aPinene'
# react_gpp_alphapinene.subsystem = 'Cytoplasm'
# react_gpp_alphapinene.lower_bound = 0.  # This is the default
# react_gpp_alphapinene.upper_bound = 1000.  # This is the default
# react_gpp_alphapinene.add_metabolites({
#     model.metabolites.get_by_id("s_0745"): -1.0, # GPP (C10H17O7P2)
#     model.metabolites.get_by_id("s_0633"): 1.0, # diphosphate/ppi (HO7P2)
#     aPinene: 1.0
# })
# react_gpp_alphapinene.gene_reaction_rule = '(APS)'
# #TODO: add associated gene(s)!

# # reaction that consumes aPinene
# react_alphapinen_con = Reaction('r_apinene_con')
# react_alphapinen_con.name = 'aPinene -> '
# react_alphapinen_con.subsystem = 'Cytoplasm'
# react_alphapinen_con.lower_bound = 0.  # This is the default
# react_alphapinen_con.upper_bound = 1000.  # This is the default
# react_alphapinen_con.add_metabolites({
#     aPinene: -1.0
# })

#reactions_list.append(react_gpp_alphapinene)
#reactions_list.append(react_alphapinen_con)

#model.add_reactions(reactions_list)

# adjust lower bounds of reactions to simulate overexpression of underlying genes

#model.reactions.get_by_id("r_gpp_aps_apinene").lower_bound = 0.2
#model.reactions.get_by_id("r_0559").lower_bound = 0.2 #erg13 

# --- step 3: set model objective (reaction(s) whose fluxes shall be maximized)
# model.objective = "r_apinene_con"
# mehrere objectives definieren:
#model.objective = {model.reactions.get_by_id("r_apinene_con"): 0.2, model.reactions.get_by_id("r_2111"): 0.8}

# --- step 4: run solver
# declare which solver to use (best one: Gurobipy)
#model.solver = 'gurobi'
#solution = model.optimize()
#print(solution)
#print(f"\nObjective value of solution: {solution.objective_value}")
# save solution object
#solution.fluxes.to_csv("out.csv", index=False)

# last step: save model
#write_yeast_model(model)   # saving, writes SBML file


