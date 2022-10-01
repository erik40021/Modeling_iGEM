from multiprocessing import freeze_support
import pwd
from cobra.io import read_sbml_model
import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
from scipy.optimize import curve_fit
from FSEOF import FSEOF
import os


def main():
    
    x = FSEOF("Data/sc_cyto_gpp_manipulated.xml", "r_2111", "r_apinene_con")
    x.find_targets(11)
    x.sort_targets()

    print(x.targets)
    x.targets.to_excel("targets_sc_cyto_gpp_FVA.xlsx")

if __name__ == "__main__":
    freeze_support()
    main()
