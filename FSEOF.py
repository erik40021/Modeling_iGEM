from cobra.io import read_sbml_model
import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
from scipy.optimize import curve_fit

class FSEOF():

    def __init__(self, path, biomassID, objectiveID):
        
        
        self.model = read_sbml_model(path)
        self.objectiveID = objectiveID
        self.biomassID = biomassID
        
        
        self.initial_fluxes = None
        self.maxProductFlux = None
        self.initialProductFlux = None
        self.optimalGrowth = None
        self.solution = None

    def calculate_intitial_fluxes(self):
        self.initial_fluxes = self.model.optimize()
        self.initialProductFlux = self.initial_fluxes.get_primal_by_id(self.objectiveID)
        self.optimalGrowth = self.initial_fluxes.objective_value
    
    def calculate_maxmimal_product_flux(self):
        self.model.objective = self.model.reactions.get_by_id(self.objectiveID)
        self.solution = self.model.optimize() 
        self.model.objective = self.model.reactions.get_by_id(self.biomassID)
        self.maxProductFlux = self.solution.objective_value
    
    def find_targets(self, steps, maxFluxCutoff=0.95):
        
        fluxes = pd.DataFrame()
        fluxCapacity = pd.DataFrame()
        enforcedFluxes = list()

        #constrain growth rate to 95% of optimal growth rate
        growthConstraint = self.model.problem.Constraint(
            self.model.reactions.get_by_id(self.biomassID).flux_expression, lb = maxFluxCutoff * self.optimalGrowth
        )
        self.model.add_cons_vars(growthConstraint)
        

        for i in range(steps - 1):
            
            #calulate current enforced flux and add it as constraint to the model    
            lb = self.initialProductFlux + (i/steps) * (self.maxProductFlux - self.initialProductFlux)
            constraint = self.model.problem.Constraint(
                self.model.reactions.get_by_id(self.objectiveID).flux_expression, lb=lb
            )
            self.model.add_cons_vars(constraint)

            #store enforced flux in enforcedFluxes list for x-values of linear regression
            enforcedFluxes.append(lb)

            #calculate FVA, Vavg, and flux capacity for all reaction with the current enforced flux
            FVA = flux_variability_analysis(self.model)
            fluxes["Average_{}".format(i)] = (FVA["maximum"] + FVA["minimum"]) / 2
            fluxCapacity["Flux_capacity_{}".format(i)] = FVA["maximum"] - FVA["minimum"]
        
        print("Enforced Fluxes: \n", enforcedFluxes)
        print("Fluxes:\n", fluxes)
        print("fluxCapacity:\n", fluxCapacity)

        #define helper functions for linear regression
        def f(x, m, b):
            return m*x+b
        
        def getSlopes(y):
            popt, _ = curve_fit(f, enforcedFluxes, y)
            return popt[0]
        
        #calculate slopes of linear regressions for fluxes and flux capactiy
        fluxes["q_slope"] = fluxes.apply(lambda y: getSlopes(y), axis=1)
        
        fluxCapacity["l_sol"] = fluxCapacity.apply(lambda y: getSlopes(y), axis=1)
        
        targets = pd.concat([fluxes["q_slope"], fluxCapacity["l_sol"]], axis=1)

        return targets

# TO-DO:
# -why step -1 ?
# -remove uneccessary prints
# - limit output of FVA? Why is is the model read so often?

        
        









