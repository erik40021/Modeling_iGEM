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
        medium[r.id] = calculate_amount_per_mass(r)                           #add only one carbon source
        model.medium = medium

        solution = r.id, model.slim_optimize(), r.name, r.reactants[0].formula
        if solution[1] > 0:
            solutions.append(solution)          #save objective value
        
    solutions.sort(reverse = True, key = lambda l: l[1])
    return solutions
    

def calculate_amount_per_mass(reaction):
    '''
    calculates amount of molecule needed to reach 10 mg.
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
            return 0
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
        print("occured for ", formula)
    print("carbon", elem_counts[0], "hydrogen", elem_counts[1], "nitrogen", elem_counts[2],
    "oxygen", elem_counts[3], "phosphor", elem_counts[4], "sulphur", elem_counts[5])
    print("mass of molecule:", (16*elem_counts[0] + 1*elem_counts[1] + 14*elem_counts[2] 
        + 16*elem_counts[3] + 31*elem_counts[4] + 32*elem_counts[5]))
    amount = 10000 / (16*elem_counts[0] + 1*elem_counts[1] + 14*elem_counts[2] 
        + 16*elem_counts[3] + 31*elem_counts[4] + 32*elem_counts[5])
    return amount

#calculate_amount_per_mass("C10H15N2O3S")

def medium_objectivevalue_xlsx(solutions, name):
    workbook = xlsxwriter.Workbook('Output/Media/' + name + '.xlsx')      #create .xlsx
    worksheet = workbook.add_worksheet()
    worksheet.write(0,0,"Carbon source reaction")         #head of table
    worksheet.write(0,1,"Objective value")
    worksheet.write(0,2,"Reaction name")
    worksheet.write(0,3,"Metabolite formula")
    for s in range(0, len(solutions)):
        worksheet.write(s+1, 0, solutions[s][0])             # write id of C-Source in first column
        worksheet.write(s+1, 1, solutions[s][1])             # objective value for this c-source
        worksheet.write(s+1, 2, solutions[s][2])             # name 
        worksheet.write(s+1, 3, solutions[s][3])             # formula
    workbook.close()
    print(name + '.xlsx was written')
