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
    @param exchange_reactions: list of all exchange reactions
    '''
    solutions=[]
    
    for r in exchange_reactions:                           #test model for every source
        medium=model.medium
        for i in exchange_reactions:                       #set all sources in medium to 0
            medium[i.id]=0
        for n in nutrients:                     #add nutrients
            medium[n]=10000
        medium[r.id]=10                            #add only one carbon source
        model.medium = medium

        solution = r.id, model.slim_optimize(), r.name, r.reactants[0].formula
        if solution[1] > 0:
            solutions.append(solution)          #save objective value
        
    solutions.sort(reverse = True, key = lambda l: l[1])
    return solutions
    

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

