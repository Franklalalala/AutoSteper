from autosteper.plotter import FullereneDataParser_Plotter


a_plotter = FullereneDataParser_Plotter()

a_plotter.plot_swr_unit(src_swr_root=r'11_to_12', group_symbol='Cl',
                        proj_ring_seq=[71, 72, 73, 74, 75, 76])
