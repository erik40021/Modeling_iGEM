from cobra import Model, Reaction, Metabolite
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model

#constants:
OVEREXPRESSION_LOWER_BOUND = 0.2
KNOCK_DOWN_HIGHER_BOUND = 500
APINENE_OBJECTIVE_COEFFICIENT = 0.0
GROWTH_OBJECTIVE_COEFFICIENT = 1.0

model = read_sbml_model("Data/yl_perox_manipulated.xml")



def simple_run():
    model.objective={model.reactions.get_by_id("Biomass_Climit"):GROWTH_OBJECTIVE_COEFFICIENT, model.reactions.get_by_id("r_apinene_con"):APINENE_OBJECTIVE_COEFFICIENT}
    solution = model.optimize()
    print(solution.objective_value)
#simple_run()
