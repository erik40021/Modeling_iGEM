from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from cobra.io import read_sbml_model
from GSMM import GSMM
from Media.medium_analysis import media_results_to_excel, run_medium_test


class Yl_perox(GSMM):
    def __init__(self) -> None:
        super().__init__()
        self.model = None
        self.build_model()


    def run_media_analysis(self, compare_factor):
        self.model.objective = self.model.reactions.get_by_id("r_apinene_con")
        exchange_reactions = self.model.exchanges
        s1, s2 = run_medium_test(self.model, exchange_reactions, "yl", compare_factor) # run analysis
        return s1, s2

    def build_model(self):
        self.model = read_sbml_model("Data/iYLI647.xml")
        # initilize metabolite objects for heterologous reactions in peroxisomes
        Acetyl_CoA = self.model.metabolites.get_by_id("accoa_x")
        Acetoacetyl_CoA = Metabolite(id="aacoa_x", formula="C25H36N7O18P3S", name="acetoacetyl-CoA[p]", charge=-4, compartment="x") 
        CoA = self.model.metabolites.get_by_id("coa_x")
        HMG_CoA = Metabolite(id="hmgcoa_x", formula="C27H39N7O20P3S", name="3-hydroxy-3-methylglutaryl-CoA[p]", charge=-5, compartment="x")
        NADPH = self.model.metabolites.get_by_id("nadph_x")
        NADP = self.model.metabolites.get_by_id("nadp_x")
        H_p = self.model.metabolites.get_by_id("h_x")
        Mevalonate = Metabolite(id="mev_x", formula="C6H11O4", name="(R)-mevalonate[p]", charge=-1, compartment="x")
        Mevalonate_P = Metabolite(id="5pmev_x", formula="C6H10O7P", name="(R)-5-phosphomevalonic acid[p]", charge=-3, compartment="x")
        Mevalonate_PP = Metabolite(id="5dpmev_x", formula="C6H10O10P2", name="(R)-5-diphosphomevalonic acid[p]", charge=-4, compartment="x")
        ATP = self.model.metabolites.get_by_id("atp_x")
        ADP = self.model.metabolites.get_by_id("adp_x")
        IPP = Metabolite(id="ipdp_x", formula="C5H9O7P2", name="isopentyl diphosphate[p]", charge=-3, compartment="x")
        DMAPP = Metabolite(id="dmpp_x", formula="C5H9O7P2", name="dimethylallyl diphosphate[p]", charge=-3, compartment="x")
        GPP = Metabolite(id="grdp_x", formula="C10H17O7P2", name="geranyl diphosphate[p]", charge=-3, compartment="x")
        Alpha_pinene = Metabolite(id="0000_x", formula="C10H16", name="(+)-alpha-pinene[p]", charge=0, compartment="x")#metabolite not excisting in any other compartment
        Phosphate = self.model.metabolites.get_by_id("pi_x")
        FPP = Metabolite(id="s_4280", formula="C15H25O7P2", name="farnesyl diphosphate[p]", charge=-3, compartment="x")
        NPP = Metabolite(id="s_4281",formula='C10H17O7P2', name='neryl diphosphate', charge=-3, compartment='x')
        CO2 = self.model.metabolites.get_by_id("co2_x")
        Diphosphate = self.model.metabolites.get_by_id("ppi_x")

        # initilize reaction objects for heterologous reactions
        Erg10_reaction = Reaction(id="MVA1", name="2.0 Acetyl-CoA --> 1.0 Acetoacetyl-CoA + 1.0 CoA", upper_bound=1000.0)
        Erg13_reaction = Reaction(id="MVA2", name="1.0 Acetyl-CoA + 1.0 Acetoacetyl-CoA --> 1.0 HMG-CoA + 1.0 CoA", upper_bound=1000.0)
        tHMGR_reaction = Reaction(id="MVA3", name="1.0 HMG-CoA + 2.0 NADPH + 2 h_p --> 1.0 Mevalonate + 2.0 NADP + 1.0 CoA", upper_bound=1000.0)

        Erg12_reaction = Reaction(id="MVA4", name="1.0 Mevalonate + 1.0 ATP --> 1.0 Mevalonat-P + 1.0 ADP + 1.0 h_p", upper_bound=1000.0)
        Erg8_reaction = Reaction(id="MVA5", name="1.0 Mevalonate-P + 1.0 ATP --> 1.0 Mevalonat-PP + 1.0 ADP", upper_bound=1000.0)
        Erg19_reaction = Reaction(id="MVA6", name="1.0 Mevalonate-PP + 1.0 ATP --> 1.0 IPP + 1.0 ADP + 1.0 PO4 + 1.0 CO2", lower_bound=None, upper_bound=1000.0) #irreversible reactions, so lower_bound = None
        Idi1_reaction = Reaction(id="MVA7", name="1.0 IPP <--> 1.0 DMAPP", upper_bound=1000.0)

        AgGPPS2_reaction = Reaction(id="MVA8", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0)
        mFPS144_GPP_reaction = Reaction(id="MVA9", name="1.0 IPP + 1.0 DMAPP --> 1.0 GPP + 1.0 Diphosphate", upper_bound=1000.0)
        mFPS144_FPP_reaction = Reaction(id="MVA10", name="1.0 IPP + 1.0 DMAPP --> 1.0 FPP + 1.0 Diphosphate", upper_bound=1000.0)
        SiNPPS1_reaction = Reaction(id="MVA11", name="1.0 IPP + 1.0 DMAPP --> 1.0 NPP + 1.0 Diphosphate", upper_bound=1000.0)
        APS_npp_reaction = Reaction(id="MVA12", name="1.0 NPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0)
        APS_gpp_reaction = Reaction(id="MVA13", name="1.0 GPP --> 1.0 Alpha-Pinene + 1.0 Diphosphate", upper_bound=1000.0)
        APinene_con_reaction = Reaction(id="r_apinene_con", name="1.0 Alpha-Pinene ->", upper_bound=1000.0)
        Diphosphate_con_reaction = Reaction(id="r_diphosphate_con", name="1.0 Diphosphate ->", upper_bound=1000.0, lower_bound=0) #artificial diphosphate consumption, to avoid accumulation

        # Add metabolites and reaction stoichiometry to reaction objects
        Erg10_reaction.add_metabolites({
            Acetyl_CoA: -2.0,
            Acetoacetyl_CoA: 1.0,
            CoA: 1.0
        })

        Erg13_reaction.add_metabolites({
            Acetyl_CoA: -1.0,
            Acetoacetyl_CoA: -1.0,
            HMG_CoA: 1.0,
            CoA: 1.0
        })

        tHMGR_reaction.add_metabolites({
            HMG_CoA: -1.0,
            NADPH: -2.0,
            H_p: -2.0,
            Mevalonate: 1.0,
            NADP: 2.0,
            CoA: 1.0
        })

        Erg12_reaction.add_metabolites({
            Mevalonate: -1.0,
            ATP: -1.0,
            Mevalonate_P: 1.0,
            ADP: 1.0,
            H_p: 1.0
        })

        Erg8_reaction.add_metabolites({
            Mevalonate_P: -1.0,
            ATP: -1.0,
            Mevalonate_PP: 1.0,
            ADP: 1.0
        })

        Erg19_reaction.add_metabolites({
            Mevalonate_PP: -1.0,
            ATP: -1.0,
            IPP: 1.0,
            ADP: 1.0,
            Phosphate: 1.0,
            CO2: 1.0
        })

        Idi1_reaction.add_metabolites({
            IPP: -1.0,
            DMAPP: 1.0
        }, reversibly=True)

        Diphosphate_con_reaction.add_metabolites({
            Diphosphate: -1.0
        })

        mFPS144_GPP_reaction.add_metabolites({
            IPP: -1.0,
            DMAPP: -1.0,
            GPP: 1.0,
            Diphosphate: 1.0
        })

        mFPS144_FPP_reaction.add_metabolites({
            GPP: -1.0,
            IPP: -1.0,
            FPP: 1.0,
            Diphosphate: 1.0
        })

        APS_gpp_reaction.add_metabolites({
            GPP: -1.0,
            Alpha_pinene: 1.0,
            Diphosphate: 1.0
        })

        APS_npp_reaction.add_metabolites({
            NPP: -1.0,
            Alpha_pinene: 1.0,
            Diphosphate: 1.0
        })

        SiNPPS1_reaction.add_metabolites({
            IPP: -1.0,
            DMAPP: -1.0,
            NPP: 1.0,
            Diphosphate: 1.0
        })

        APinene_con_reaction.add_metabolites({
            Alpha_pinene: -1.0
        })

        # Add reactions to model
        reactionlist = [
            Erg10_reaction,
            Erg13_reaction,
            tHMGR_reaction,
            Erg12_reaction,
            Erg8_reaction,
            Erg19_reaction,
            Idi1_reaction,
            mFPS144_GPP_reaction, #7
            mFPS144_FPP_reaction,
            AgGPPS2_reaction,
            SiNPPS1_reaction,
            APS_gpp_reaction,
            APS_npp_reaction,
            APinene_con_reaction,
            Diphosphate_con_reaction
        ]
        for reaction in reactionlist:
            self.model.add_reaction(reaction)
