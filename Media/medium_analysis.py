from operator import index
from yl_cyto_old import run_medium_test          #still needs rund_medium_test eqivalents from other Yeasts
from media_utils import medium_objectivevalue_xlsx
import pandas as pd
Food=[                      #list carbon sources here   #not used atm
"EX_Fat_LPAREN_e_RPAREN_",
"EX_glc_LPAREN_e_RPAREN_",
"EX_inost_LPAREN_e_RPAREN_",
"EX_tre_LPAREN_e_RPAREN_",
"EX_xyl_D_LPAREN_e_RPAREN_",
"EX_fru_LPAREN_e_RPAREN_",
"EX_glyc_LPAREN_e_RPAREN_"
]
NUTRIENTS_YL=[                 #list other nutrients here
"EX_h2o_LPAREN_e_RPAREN_",
"EX_h_LPAREN_e_RPAREN_",
"EX_k_LPAREN_e_RPAREN_",
"EX_na1_LPAREN_e_RPAREN_",
"EX_nh4_LPAREN_e_RPAREN_",
"EX_o2_LPAREN_e_RPAREN_",
"EX_pi_LPAREN_e_RPAREN_",
"EX_so4_LPAREN_e_RPAREN_",
]
NUTRIENTS_SC=[
'r_1654',
'r_1832',
'r_1861',
'r_1992',
'r_2005',
'r_2020',
'r_2049',
'r_2060',
'r_2100',
'r_4593',
'r_4594',
'r_4595',
'r_4596',
'r_4597',
'r_4600'
]

def calculate_candidates(import_reactions, nutrients, filename):
    Solutions=run_medium_test(import_reactions,nutrients)       #run analysis

    medium_objectivevalue_xlsx(import_reactions,Solutions, filename)      #print to xlsx
