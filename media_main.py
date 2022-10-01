

from Media.medium_analysis import media_results_to_excel, plot_oxygen_curves
from sc_cytosol import Sc_cyto
from sc_peroxisome import Sc_perox
from yl_cytosol import Yl_cyto
from yl_peroxisome import Yl_perox


def run_all_analyses():
    #create model objects
    sc_cyto_gpp = Sc_cyto("GPP")
    sc_cyto_npp = Sc_cyto("NPP")
    sc_perox = Sc_perox()
    yl_cyto = Yl_cyto()
    yl_perox = Yl_perox()

    # run analyses with equal mass
    s1, s2 = sc_cyto_gpp.run_media_analysis("mass")
    # media_results_to_excel(s1, "Csource/sc_cyto_gpp_csource")
    media_results_to_excel(s2, "Oxygen/sc_cyto_gpp_oxygen_equalmass")
    s1, s2 = sc_cyto_npp.run_media_analysis("mass")
    # media_results_to_excel(s1, "Csource/sc_cyto_npp_csource")
    media_results_to_excel(s2, "Oxygen/sc_cyto_npp_oxygen_equalmass")
    s1, s2 = sc_perox.run_media_analysis("mass")
    # media_results_to_excel(s1, "Csource/sc_perox_csource")
    media_results_to_excel(s2, "Oxygen/sc_perox_oxygen_equalmass")
    s1, s2 = yl_cyto.run_media_analysis("mass")
    # media_results_to_excel(s1, "Csource/yl_cyto_csource")
    media_results_to_excel(s2, "Oxygen/yl_cyto_oxygen_equalmass")
    s1, s2 = yl_perox.run_media_analysis("mass")
    # media_results_to_excel(s1, "Csource/yl_perox_csource")
    media_results_to_excel(s2, "Oxygen/yl_perox_oxygen_equalmass") 

    # and (for oxygen) also with equal carbon
    s1, s2 = sc_cyto_gpp.run_media_analysis("carbon")
    media_results_to_excel(s2, "Oxygen/sc_cyto_gpp_oxygen_equalcarbon")
    s1, s2 = sc_cyto_npp.run_media_analysis("carbon")
    media_results_to_excel(s2, "Oxygen/sc_cyto_npp_oxygen_equalcarbon")
    s1, s2 = sc_perox.run_media_analysis("carbon")
    media_results_to_excel(s2, "Oxygen/sc_perox_oxygen_equalcarbon")
    s1, s2 = yl_cyto.run_media_analysis("carbon")
    media_results_to_excel(s2, "Oxygen/yl_cyto_oxygen_equalcarbon")
    s1, s2 = yl_perox.run_media_analysis("carbon")
    media_results_to_excel(s2, "Oxygen/yl_perox_oxygen_equalcarbon")


def plot_oxygen_results(verbose):
    plot_oxygen_curves("sc_cyto_gpp_oxygen_equalmass", verbose, title=r"$\it{S. cerivisiae}$ cytosol")
    plot_oxygen_curves("sc_cyto_npp_oxygen_equalmass", verbose, title=r"$\it{S. cerivisiae}$ cytosol")   
    plot_oxygen_curves("sc_perox_oxygen_equalmass", verbose, title=r"$\it{S. cerivisiae}$ peroxisome")
    plot_oxygen_curves("yl_cyto_oxygen_equalmass", verbose, title=r"$\it{Y. lipolytica}$ cytosol")
    plot_oxygen_curves("yl_perox_oxygen_equalmass", verbose, title=r"$\it{Y. lipolytica}$ peroxisome")

    plot_oxygen_curves("sc_cyto_gpp_oxygen_equalcarbon", verbose)
    plot_oxygen_curves("sc_cyto_npp_oxygen_equalcarbon", verbose)   
    plot_oxygen_curves("sc_perox_oxygen_equalcarbon", verbose)
    plot_oxygen_curves("yl_cyto_oxygen_equalcarbon", verbose)
    plot_oxygen_curves("yl_perox_oxygen_equalcarbon", verbose)

#run_all_analyses()
plot_oxygen_results(True)
