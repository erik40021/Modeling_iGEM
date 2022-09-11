from operator import index
import pandas as pd
import xlsxwriter

# Food_YL=[                      #list carbon sources here   #not used atm
# "EX_Fat_LPAREN_e_RPAREN_",
# "EX_glc_LPAREN_e_RPAREN_",
# "EX_inost_LPAREN_e_RPAREN_",
# "EX_tre_LPAREN_e_RPAREN_",
# "EX_xyl_D_LPAREN_e_RPAREN_",
# "EX_fru_LPAREN_e_RPAREN_",
# "EX_glyc_LPAREN_e_RPAREN_"
# ]

NUTRIENTS_YL=[                 # all essential non-carbon 'nutrients' for YL
"EX_h2o_LPAREN_e_RPAREN_",
"EX_h_LPAREN_e_RPAREN_",
"EX_k_LPAREN_e_RPAREN_",
"EX_na1_LPAREN_e_RPAREN_",
"EX_nh4_LPAREN_e_RPAREN_",
"EX_o2_LPAREN_e_RPAREN_",
"EX_pi_LPAREN_e_RPAREN_",
"EX_so4_LPAREN_e_RPAREN_",
]

NUTRIENTS_SC=['r_1654','r_1832','r_1861','r_1992','r_2005','r_2020','r_2049','r_2060',
    'r_2100','r_4593','r_4594','r_4595','r_4596','r_4597','r_4600']




# -> for example usage, see 'run_media_analysis()' in sc_cytosol.py

def run_medium_test(model, exchange_reactions, nutrients):
    '''
    calculates list of media with their respective performances regarding the pre-set (!)
    objective value in model and stores in excel file
    @param exchange_reactions: list of all exchange reaction objects
    '''
    solutions=[]
    
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
        solutions.append(solution)          #save objective value
        
    solutions.sort(reverse = True, key = lambda l: l[2])
    solutions.insert(0, ["id", "reaction name", "apinene flux [equal mass]", "mmol used", "apinene flux [equal carbon]",
        "mmol used", "carbon efficiency", "metabolite formula"])
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

#calculate_amount_per_mass("C10H15N2O3S")

def medium_objectivevalue_xlsx(solutions, name):
    workbook = xlsxwriter.Workbook('Output/Media/' + name + '.xlsx')      #create .xlsx
    worksheet = workbook.add_worksheet()
    # for i in range(0, len(solutions[0])):
    #     worksheet.write(0,i,solutions[0][i])         #head of table
    for i in range(0, len(solutions)):
        for j in range(0, len(solutions[i])):
            worksheet.write(i, j, solutions[i][j])
    workbook.close()
    print(name + '.xlsx was written')
