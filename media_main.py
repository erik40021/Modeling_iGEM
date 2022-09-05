

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
    #run analyses
    sc_cyto_gpp.run_media_analysis("sc_cyto_gpp_equalmass")
    sc_cyto_npp.run_media_analysis("sc_cyto_npp_equalmass")
    sc_perox.run_media_analysis("sc_perox_equalmass")
    yl_cyto.run_media_analysis("yl_cyto_equalmass")
    yl_perox.run_media_analysis("yl_perox_equalmass")
    
run_all_analyses()