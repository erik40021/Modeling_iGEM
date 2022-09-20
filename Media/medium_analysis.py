from cProfile import label
from cmath import nan
import math
from operator import index
import pandas as pd
import xlsxwriter
from matplotlib import pyplot as plt

NUTRIENTS_YL=[                 # all essential non-carbon 'nutrients' for YL
"EX_o2_LPAREN_e_RPAREN_",
"EX_h2o_LPAREN_e_RPAREN_",
"EX_h_LPAREN_e_RPAREN_",
"EX_k_LPAREN_e_RPAREN_",
"EX_na1_LPAREN_e_RPAREN_",
"EX_nh4_LPAREN_e_RPAREN_",
"EX_pi_LPAREN_e_RPAREN_",
"EX_so4_LPAREN_e_RPAREN_",
]

NUTRIENTS_SC=['r_1992','r_1654','r_1832','r_1861','r_2005','r_2020','r_2049','r_2060',
    'r_2100','r_4593','r_4594','r_4595','r_4596','r_4597','r_4600']




# -> for example usage, see 'run_media_analysis()' in sc_cytosol.py

def run_medium_test(model, exchange_reactions, organism, compare_factor):
    '''
    calculates list of media with their respective performances regarding the pre-set (!)
    objective value in model and stores in excel file
    @param exchange_reactions: list of all exchange reaction objects
    '''
    if organism == "sc":
        nutrients = NUTRIENTS_SC
    elif organism == "yl":
        nutrients = NUTRIENTS_YL
    if compare_factor == "mass":
        compare_factor_index = 2
    elif compare_factor == "carbon":
        compare_factor_index = 4
    carbon_source_solution = []
    
    for r in exchange_reactions:                           #test model for every source
        medium=model.medium
        if r.id in nutrients:
            continue
        for i in exchange_reactions:                       #set all sources in medium to 0
            medium[i.id]=0
        for n in nutrients:                     #add nutrients
            medium[n] = 10000
        
        equal_weight, formula, equal_carbon = calculate_amount_per_mass(r)                           #add only one carbon source
        medium[r.id] = equal_weight
        model.medium = medium
        solution_equal_mass = model.slim_optimize()
        medium[r.id] = equal_carbon
        model.medium = medium
        solution_equal_carbon_number = model.slim_optimize()
        if solution_equal_mass is None or not solution_equal_mass > 0:
            continue

        solution = (r.id, r.name, solution_equal_mass, equal_weight, solution_equal_carbon_number, equal_carbon, solution_equal_carbon_number/60, formula)
        carbon_source_solution.append(solution)          #save objective value
        
    carbon_source_solution.sort(reverse = True, key = lambda l: l[compare_factor_index])
    oxygen_solution = run_oxygen_test(model, carbon_source_solution, nutrients, compare_factor_index)
    carbon_source_solution.insert(0, ["id", "reaction name", "apinene flux [equal mass]", "mmol used", "apinene flux [equal carbon]",
        "mmol used", "carbon efficiency", "metabolite formula"])
    return carbon_source_solution, oxygen_solution
    
def run_oxygen_test(model, carbon_source_solution, nutrients, compare_factor_index, oxygen_analysis_number=20, upper_limit=100, resolution=20):
    oxygen_exchange_id = nutrients[0]
    carbon_source_solution.sort(reverse = True, key = lambda l: l[compare_factor_index])
    selection_to_be_analysed = carbon_source_solution[:oxygen_analysis_number]
    reaction_ids = [sublist[0] for sublist in selection_to_be_analysed]
    reaction_names = [sublist[1] for sublist in selection_to_be_analysed]
    amounts_mmol = [sublist[compare_factor_index + 1] for sublist in selection_to_be_analysed]
    step_size = upper_limit/resolution
    solutions = []
    legend = ["name"]
    for i in range(0, resolution + 1):
        legend.append(str(upper_limit - i*step_size))
    solutions.append(legend)
    medium = model.medium
    for m in medium:
        medium[m] = 0           # reset all previous media components to 0
    for n in nutrients:
        medium[n] = 10000
    for i in range(0, len(reaction_ids)):
        medium[reaction_ids[i]] = amounts_mmol[i]   # set carbon source to previously calculated, comparable amount
        solution = [reaction_names[i]]
        # medium[oxygen_exchange_id] = 10000          # reset oxygen for first solving to inital value (10000) 
        # model.medium = medium
        # solution.append(model.slim_optimize())      # first entry: flux for initial oxygen 
        for j in range(0, resolution + 1):
            medium[oxygen_exchange_id] = upper_limit - j*step_size
            model.medium = medium
            s = model.slim_optimize()
            if s is None or math.isnan(s):
                solution.append(0)
                continue
            solution.append(s)
        solutions.append(solution)
        medium[reaction_ids[i]] = 0     # reset tested carbon source to 0
    return solutions


def calculate_amount_per_mass(reaction):
    '''
    note: formula of metabolite is expected to have HiII-system compliant element order, e.g.: C -> H -> N -> O -> P -> S
    '''
    try:
        formula = reaction.reactants[0].formula
    except:
        # in case reaction has no reactant, look for product (=same metabolite)
        formula = reaction.products[0].formula
    if "R" in formula:      # filter out rests
        formula = formula.split("R")[0]
    elements = ["H", "N", "O", "P", "S", "X"]
    elem_counts = [0, 0, 0, 0, 0, 0]
    try:
        divstr = formula.split("C")
        if len(divstr) == 1: # filter out formulas without carbon
            return 0, formula, 0
        i = 0
        while(i < len(elem_counts)):
            divstr = divstr[1].split(elements[i])
            if len(divstr) == 1:
                if i > len(elements)-2: #no more dividers there
                    if len(divstr[0]) == 0:
                        elem_counts[i] = 1
                    else:
                        elem_counts[i] = int(divstr[0])
                    break
                k = 1
                partitioning_done = True
                for e in elements[i+1:len(elements)-1]:
                    k += 1
                    if e in divstr[0]:
                        divstr = divstr[0].split(e)
                        if len(divstr[0]) == 0:
                            elem_counts[i] = 1
                        else:
                            elem_counts[i] = int(divstr[0])
                        i += k
                        partitioning_done = False
                        break
                if partitioning_done:
                    if len(divstr[0]) == 0:
                        elem_counts[i] = 1
                    else:
                        elem_counts[i] = int(divstr[0])
                    break
            elif len(divstr[0]) == 0:
                elem_counts[i] = 1
                i += 1
            else:
                elem_counts[i] = int(divstr[0])
                i += 1
    except IndexError as e:
        print(e)
        print("occured for ", formula)
    except:
        print("Error occured for ", formula)
    # print("Result for ", formula, ":  C:", elem_counts[0], "H:", elem_counts[1], "N:", elem_counts[2],
    # "O:", elem_counts[3], "P:", elem_counts[4], "S:", elem_counts[5])
    # print("mass of molecule:", (16*elem_counts[0] + 1*elem_counts[1] + 14*elem_counts[2] 
    #     + 16*elem_counts[3] + 31*elem_counts[4] + 32*elem_counts[5]))

    amount = 1800 / (12*elem_counts[0] + 1*elem_counts[1] + 14*elem_counts[2] 
        + 16*elem_counts[3] + 31*elem_counts[4] + 32*elem_counts[5])
    carbon = 60 / elem_counts[0]    # means: we test with 60 mmol of carbon atoms per gdcw/h (= 60*6*10^20 = 3.6*10^22 carbon atoms)
    return amount, formula, carbon


def plot_oxygen_curves(filename, verbose=False):
    # TODO: adapt to plot-stylsheet conventions (colors!)
    df = pd.read_excel('Output/Media/Oxygen/' + filename + '.xlsx', index_col=0)
    
    metabolites = list(df.index)
    oxygenUptake = list(df.columns)

    #create figure
    fig = plt.figure(figsize=(14,10))
    ax = fig.add_subplot()

    for metabolite in metabolites:
        ax.plot(oxygenUptake, df.loc[metabolite].values, label=metabolite)
    
    ax.set_xlabel("Oxygen uptake [mmol/gcdw/h]")
    ax.set_ylabel("Alpha-pinene Flux [mmol/gdcw/h]")
    ax.legend()
    plt.savefig('Output/Media/Oxygen/' + filename + '.png')
    if not verbose:
        plt.show()



def media_results_to_excel(solutions, filename):
    workbook = xlsxwriter.Workbook('Output/Media/' + filename + '.xlsx')      #create .xlsx
    worksheet = workbook.add_worksheet()
    for i in range(0, len(solutions)):
        for j in range(0, len(solutions[i])):
            worksheet.write(i, j, solutions[i][j])
    workbook.close()
    print(filename + '.xlsx was written')
