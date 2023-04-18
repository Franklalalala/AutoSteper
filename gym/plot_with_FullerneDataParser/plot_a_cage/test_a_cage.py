import matplotlib.pyplot as plt
from ase.build.molecule import molecule
from autosteper.plotter import FullereneDataParser_Plotter


a_plotter = FullereneDataParser_Plotter()
####################### Plot on customized figure ##################################
C60 = molecule(name='C60')
fig = plt.figure(figsize=[10, 10])
ax = fig.add_subplot(111)
a_plotter.load_cage(atoms=C60, ax=ax)
plt.savefig('customized_figure.png', dpi=400)
plt.close()

####################### Simple Plot ##################################
# If no figure specified, plotter will set a new figure with figure size [10, 10]
C60 = molecule(name='C60')
a_plotter.load_cage(atoms=C60)
plt.savefig('simple.png', dpi=400)
plt.close()

####################### see atom labels ##################################
C60 = molecule(name='C60')
a_plotter.load_cage(atoms=C60, show_C_label=True, C_label_transparency=0.5, C_label_color='orange')
plt.savefig('see_labels.png', dpi=400)
plt.close()

####################### Re-plot to project from a specific ring ##################################
# assume [48, 49, 50, 51, 52, 45] is an ideal ring to project from
C60 = molecule(name='C60')
a_plotter.load_cage(atoms=C60, show_C_label=True, C_label_transparency=0.5, C_label_color='orange',
                    proj_ring_seq=[48, 49, 50, 51, 52, 45])
plt.savefig('re_plot.png', dpi=400)
plt.close()
# See labels on the out-most ring has changed to the specified ring,
# that means this ring is the start of projection.
