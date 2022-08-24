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

model = read_sbml_model("Data/yl_cyto_manipulated.xml")

def run_medium_test(CSource,Nutrients):
    Solutions=[]
    #objective function
    model.objective={model.reactions.get_by_id("Biomass_Climit"):GROWTH_OBJECTIVE_COEFFICIENT, model.reactions.get_by_id("aPinene_ex"):APINENE_OBJECTIVE_COEFFICIENT}
    
    for f in CSource:                          #test model for every food
        medium=model.medium
        for i in CSource:                      #set food in medium to 0
            medium[i]=0
        for n in Nutrients:                 #add nutrients
            medium[n]=10000
        medium[f]=10                   #add only one food
        model.medium = medium

        solution = model.optimize()
        Solutions.append(solution.objective_value)          #save objective value
        
    print(Solutions)
    return Solutions

def simple_run():
    model.objective={model.reactions.get_by_id("Biomass_Climit"):GROWTH_OBJECTIVE_COEFFICIENT, model.reactions.get_by_id("aPinene_ex"):APINENE_OBJECTIVE_COEFFICIENT}
    solution = model.optimize()
    print(solution.objective_value)

