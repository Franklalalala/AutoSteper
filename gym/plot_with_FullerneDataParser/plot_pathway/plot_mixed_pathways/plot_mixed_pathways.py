# These pathways have mixed end-state products
import os

import matplotlib.pyplot as plt
from autosteper.plotter import FullereneDataParser_Plotter
from autosteper.tools import strip_extraFullerene


src_pathways = r'src_pathways'
new_pathway_workbase = r'new_pathway_workbase'
max_adduct_path = r'10_addons_1.xyz'
###################################################################################
os.makedirs(new_pathway_workbase, exist_ok=True)
new_pathway_workbase = os.path.abspath(new_pathway_workbase)
max_adduct_path = os.path.abspath(max_adduct_path)

a_plotter = FullereneDataParser_Plotter()
os.chdir(src_pathways)
for a_pathway in os.listdir(r'./'):
    a_plotter.plot_pathway_unit(src_pathway_root=a_pathway,
                                # 1, 2 are the lowest-energy pathway generated by AutoSteper for #4169 C66Cl10
                                # 34, 37 are the pathways proposed in article: DOI: 10.1038/nchem.549
                                new_pathway_workbase=new_pathway_workbase,
                                group_symbol='Cl',
                                proj_ring_seq=[61, 62, 63, 64, 65, 66],
                                is_match_max_adduct = True,
                                max_adduct_path=max_adduct_path,
                                diff_len=1
                                )
