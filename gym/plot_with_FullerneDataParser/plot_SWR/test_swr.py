from autosteper.plotter import FullereneDataParser_Plotter


a_plotter = FullereneDataParser_Plotter()
# Make sure the projection ring stays in the *target* atoms. DOes NOT CONTAIN THE ROTATED BOND IN SWR.
# To reproduce pictures in paper, try [3,7,12,14,16,18].
a_plotter.plot_swr(src_swr_root=r'q_11_to_tgt_12', group_symbol='Cl',
                        proj_ring_seq=[1, 5, 9, 13, 15, 17])
