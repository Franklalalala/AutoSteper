# This script is applied to scenario when one is not interested in mapping end-state isomers,
# or these pathways have identical end-state products.
import os

import matplotlib.pyplot as plt
from autosteper.plotter import FullereneDataParser_Plotter
from autosteper.tools import strip_extraFullerene


src_pathways = r'src_pathways'
new_pathway_workbase = r'new_pathway_workbase'
###################################################################################
os.makedirs(new_pathway_workbase, exist_ok=True)
new_pathway_workbase = os.path.abspath(new_pathway_workbase)

a_plotter = FullereneDataParser_Plotter()
os.chdir(src_pathways)
for a_pathway in os.listdir(r'./'):
    a_plotter.plot_pathway_unit(src_pathway_root=a_pathway,
                                # this is the lowest-energy pathway generated by AutoSteper for #4169 C66Cl10
                                new_pathway_workbase=new_pathway_workbase,
                                group_symbol='Cl',
                                proj_ring_seq=[61, 62, 63, 64, 65, 66])


