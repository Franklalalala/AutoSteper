import matplotlib.pyplot as plt
from autosteper.plotter import FullereneDataParser_Plotter
from autosteper.tools import get_low_e_ranks, strip_extraFullerene
import pandas as pd


a_plotter = FullereneDataParser_Plotter()
# Here we take example on dihept-C66H4
info = pd.read_pickle(r'path/to/passed_info.pickle')
cutoff_para = {
    'mode': 'rank',
    'rank': 5
}
for a_rank in get_low_e_ranks(e_arr=info['energy'], para=cutoff_para):
    a_xyz_path = info['xyz_path'][a_rank]
    cage, addon_set = strip_extraFullerene(coord_file_path=a_xyz_path, group='H')
    a_plotter.load_cage(atoms=cage)
    a_plotter.load_addons(addon_set=addon_set)
    plt.savefig(f'rank_{a_rank+1}_2D.png', dpi=400)
    plt.close()

a=1

