from yl_cytosol import run_full
from Model.model_utils import medium_objectivevalue_xlsx

Food=[                      #list carbon sources here
"EX_glc_LPAREN_e_RPAREN_",
"EX_inost_LPAREN_e_RPAREN_",
"trehalose_c_tp",
"EX_xyl_D_LPAREN_e_RPAREN_"
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

Solutions=run_full(Food,Nutrients)              #run analysis
medium_objectivevalue_xlsx(Food,Solutions)      #print to xlsx