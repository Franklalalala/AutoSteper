import matplotlib.pyplot as plt
from autosteper.plotter import FullereneDataParser_Plotter
from autosteper.tools import strip_extraFullerene

a_plotter = FullereneDataParser_Plotter()
pristine_cage, addon_set = strip_extraFullerene(coord_file_path=r'dihept_C66H4.xyz')
# ####################### Pure cage ##################################
a_plotter.load_cage(atoms=pristine_cage)
plt.savefig('no_addon.png', dpi=400)
plt.close()

####################### Simple Plot ##################################
a_plotter.load_cage(atoms=pristine_cage)
a_plotter.load_addons(addon_set=addon_set, group_symbol='H')
plt.savefig('simple_plot.png', dpi=400)
plt.close()

# ####################### see atom labels ##################################
a_plotter.load_cage(atoms=pristine_cage, show_C_label=True, C_label_transparency=0.5, C_label_color='orange')
a_plotter.load_addons(addon_set=addon_set, group_symbol='H', show_addon_nums=True)
plt.savefig('plot_with_label.png', dpi=400)
plt.close()
# ####################### try all hexagons, pick your favorite ##################################
a_plotter.try_every_hexagon(dump_pic_folder='./try_all_hexagons',
                            addon_set=addon_set,
                            atoms=pristine_cage, group_symbol='H')

####################### Re-plot to project from a specific ring ##################################
# assume [10, 30, 31, 33, 21, 50] is an ideal ring to project from
a_plotter.load_cage(atoms=pristine_cage, show_C_label=True, C_label_transparency=0.5, C_label_color='orange',
                    proj_ring_seq=[10, 30, 31, 33, 21, 50])
a_plotter.load_addons(addon_set=addon_set, group_symbol='H', show_addon_nums=True)
plt.savefig('re_plot_with_label.png', dpi=400)
plt.close()

####################### Clean  ##################################
a_plotter.load_cage(atoms=pristine_cage, proj_ring_seq=[10, 30, 31, 33, 21, 50])
a_plotter.load_addons(addon_set=addon_set, group_symbol='H')
plt.savefig('Clean.png', dpi=400)
plt.close()

####################### Plot pathways for publication  ##################################
pristine_cage, addon_set = strip_extraFullerene(coord_file_path=r'C66Cl10_4169_exp.xyz')
# This sequence could be observed in ./../plot_pathway/plot_one_pathway/re_label/path_rank_1
original_seq = [1, 2, 3, 12, 7, 23, 38, 35, 18, 15]
new_seq = list(range(1, 11, 1))
replace_addon_map = dict(zip(original_seq, new_seq))

a_plotter.load_cage(atoms=pristine_cage, proj_ring_seq=[61, 62, 63, 64, 65, 66], sphere_ratio=4)
a_plotter.load_addons(addon_set=addon_set, group_symbol='Cl', show_addon_nums=True, fontsize=25,
                      addon_label_size=750, addon_nums_rotation=-17, replace_addon_map=replace_addon_map)
plt.savefig('pathway_on_one_pic.png', dpi=400)
plt.close()



