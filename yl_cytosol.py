from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model

from GSMM import GSMM
from Media.medium_analysis import NUTRIENTS_YL, medium_objectivevalue_xlsx, run_medium_test

#constants:
OVEREXPRESSION_LOWER_BOUND = 0.0
KNOCK_DOWN_HIGHER_BOUND = 500
APINENE_OBJECTIVE_COEFFICIENT = 1.0
GROWTH_OBJECTIVE_COEFFICIENT = 1.0 - APINENE_OBJECTIVE_COEFFICIENT
GROWTH_LOWER_BOUND = 0.0

class Yl_cyto(GSMM):
    def __init__(self) -> None:
        super().__init__()
        self.model = None
        self.build_model()


    def run_media_analysis(self, filename="yl_cyto_equalmass"):
        self.model.objective = self.model.reactions.get_by_id("aPinene_ex")
        exchange_reactions = self.model.exchanges
        solutions = run_medium_test(self.model, exchange_reactions, NUTRIENTS_YL) # run analysis
        medium_objectivevalue_xlsx(solutions, filename) #store as excel file

    def build_model(self):
        self.model = read_sbml_model("Data/iYLI647.xml")
        reactionlist =[]
        aPinene = Metabolite('aPinene', formula='C10H16', name='Alphapinene', compartment='c')
        NPP = Metabolite('NPP', formula='C5H9O7P2', name='NPP', compartment='c' )
        GPP = self.model.metabolites.get_by_id("grdp_c")
        FPP = self.model.metabolites.get_by_id("frdp_c")
        IPP = self.model.metabolites.get_by_id("ipdp_c")
        DMAPP = self.model.metabolites.get_by_id("dmpp_c")
        HMGR = self.model.metabolites.get_by_id("hmgcoa_c")
        MEV = self.model.metabolites.get_by_id("mev_R_c")
        NADPH = self.model.metabolites.get_by_id("nadph_c")
        NADP = self.model.metabolites.get_by_id("nadp_c")
        H = self.model.metabolites.get_by_id("h_c")
        COA = self.model.metabolites.get_by_id("coa_c")
        Diphosphate = self.model.metabolites.get_by_id("ppi_c")
        Fat = self.model.metabolites.get_by_id("triglyc_SC_e")

        #EX_Fat
        EX_Fat_reaction=Reaction(id="EX_Fat_LPAREN_e_RPAREN_",name="-->Fat extracellular",subsystem="Extracellular")
        EX_Fat_reaction.add_metabolites({
            Fat: 1.0
        })
        reactionlist.append(EX_Fat_reaction)

        #AP_Syn
        APS_GPP_reaction = Reaction(id="aPinen_gpp_s", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
        APS_NPP_reaction = Reaction(id="aPinen_fpp_s", name="1.0 NPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
        Extract_aPinene =Reaction(id="aPinene_ex", name="extract Alpha-Pinene", subsystem="Cytoplasm", upper_bound=1000.0)
        APS_GPP_reaction.add_metabolites({
            GPP: -1.0,
            aPinene: 1.0,
            Diphosphate: 1.0
        })
        APS_NPP_reaction.add_metabolites({
            NPP: -1.0,
            aPinene: 1.0,
            Diphosphate: 1.0
        })
        Extract_aPinene.add_metabolites({
            aPinene: -1.0
        })
        reactionlist.append(APS_GPP_reaction)
        reactionlist.append(APS_NPP_reaction)
        reactionlist.append(Extract_aPinene)

        #MEV_Pathway
        tHMGR_reaction = Reaction(id="tHMGR", name="1.0 HMG-CoA + 2.0 NADPH + 2 h_c --> 1.0 Mevalonate + 2.0 NADP + 1.0 CoA",subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND,upper_bound=1000.0)
        tHMGR_reaction.add_metabolites({
            HMGR:-1.0,
            NADPH:-2.0,
            H: -2.0,
            MEV: 1.0,
            NADP: 2.0,
            COA: 1.0
        })
        reactionlist.append(tHMGR_reaction)

        #GPP_Pathway
        GPPS_reaction = Reaction(id="GPP_s", name= "1.0 IPP+1.0 DMAPP--> 1.0 GPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound= OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
        mFPPS_gpp_reaction = Reaction(id="FPP_GPP_s",name="1.0 IPP+1.0 DMAPP--> 1.0 GPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND,upper_bound=1000.0)
        mFPPS_fpp_reaction = Reaction(id="FPP_FPP_s",name="1.0 GPP+1.0 DMAPP--> 1.0 FPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound=OVEREXPRESSION_LOWER_BOUND,upper_bound=1000.0)
        GPPS_reaction.add_metabolites({
            IPP: -1.0,
            DMAPP: -1.0,
            Diphosphate: 1.0,
            GPP: 1.0
        })
        mFPPS_gpp_reaction.add_metabolites({
            IPP: -1.0,
            DMAPP: -1.0,
            Diphosphate: 1.0,
            GPP: 1.0
        })
        mFPPS_gpp_reaction.gene_reaction_rule='mFPPs'
        mFPPS_fpp_reaction.add_metabolites({
            GPP: -1.0,
            DMAPP: -1.0,
            Diphosphate: 1.0,
            FPP: 1.0
        })
        mFPPS_fpp_reaction.gene_reaction_rule='mFPPs'

        reactionlist.append(GPPS_reaction)
        reactionlist.append(mFPPS_gpp_reaction)
        reactionlist.append(mFPPS_fpp_reaction)

        #NPP_Pathway
        SINPPS_reaction = Reaction(id="SINPP_s", name= "1.0 IPP+1.0 DMAPP--> 1.0 NPP + 1.0 Diphosphate",subsystem="Cytoplasm",lower_bound= OVEREXPRESSION_LOWER_BOUND, upper_bound=1000.0)
        SINPPS_reaction.add_metabolites({
            IPP: -1.0,
            DMAPP: -1.0,
            Diphosphate: 1.0,
            NPP: 1.0
        })
        reactionlist.append(SINPPS_reaction)

        #erg20_knockdown
        self.model.reactions.get_by_id("DMATT").higher_bound=KNOCK_DOWN_HIGHER_BOUND
        self.model.reactions.get_by_id("GRTT").higher_bound=KNOCK_DOWN_HIGHER_BOUND

        #erg13_overexpression
        self.model.reactions.get_by_id("HMGCOAS").lower_bound=OVEREXPRESSION_LOWER_BOUND

        for reaction in reactionlist:
            self.model.add_reaction(reaction)

    def add_formulas(self):
        import pandas as pd
        df = pd.read_excel('Output/Media/name_Modeling.xlsx', sheet_name=0) # can also index sheet by name or fetch all sheets
        formulas = df['formula'].tolist()
        ids = df['id'].tolist()
        for i in range(0, len(ids)):
            try:
                self.model.reactions.get_by_id(ids[i]).reactants[0].formula = formulas[i]
            except:
                self.model.reactions.get_by_id(ids[i]).products[0].formula = formulas[i]
        write_sbml_model(self.model,"Data/iYLI647+.xml") #remove + in name if model is correct

