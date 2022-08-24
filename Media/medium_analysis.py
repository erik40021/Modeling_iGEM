from operator import index
from yl_cytosol import run_medium_test          #still needs rund_medium_test eqivalents from other Yeasts
from Model.model_utils import medium_objectivevalue_xlsx
import pandas as pd
Food=[                      #list carbon sources here
"EX_glc_LPAREN_e_RPAREN_",
"EX_inost_LPAREN_e_RPAREN_",
"EX_tre_LPAREN_e_RPAREN_",
"EX_xyl_D_LPAREN_e_RPAREN_",
"EX_fru_LPAREN_e_RPAREN_",
"EX_glyc_LPAREN_e_RPAREN_"
]
Nutrients=[                 #list other nutrients here
"EX_h2o_LPAREN_e_RPAREN_",
"EX_h_LPAREN_e_RPAREN_",
"EX_k_LPAREN_e_RPAREN_",
"EX_na1_LPAREN_e_RPAREN_",
"EX_nh4_LPAREN_e_RPAREN_",
"EX_o2_LPAREN_e_RPAREN_",
"EX_pi_LPAREN_e_RPAREN_",
"EX_so4_LPAREN_e_RPAREN_",
]

List=[]
all=pd.read_excel("Output/EX_Reactions.xlsx")          #List is used as C-Source and is generated from all EX_reactions
cnt=0
while cnt<len(all):
    List.append(all.iat[cnt,0])
    cnt+=1

Solutions,Formula=run_medium_test(List,Nutrients)       #run analysis

medium_objectivevalue_xlsx(List,Solutions)      #print to xlsx
