import os
from autosteper.tools import read_xtb_log
from autosteper.cage import name2seq, Cage
from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from ase.atoms import Atoms
import numpy as np
from ase.io.gaussian import read_gaussian_out


class Checker():
    def __init__(self, chk_skin: float, group: str, cage: Cage):
        self.skin = chk_skin
        self.group = group
        self.cage = cage


    def check_ADJ(self, dok_matrix):
        for idx, an_atom in enumerate(self.last_image):
            dok_adj = dok_matrix[idx]
            adj_array = dok_adj.tocoo().col
            len_adj_arr = len(adj_array)
            if an_atom.symbol == self.cage.symbol:
                if len_adj_arr == 4:
                    flag = 0
                    for ii in adj_array:
                        if self.last_image[ii].symbol in [self.group, 'H', 'O']:
                            flag = 1
                            break
                    if flag:
                        continue
                    else:
                        return 3
                elif len_adj_arr < 3:
                    return 4
                elif len_adj_arr == 3:
                    for ii in adj_array:
                        if self.last_image[ii].symbol in [self.group, 'H', 'O']:
                            return 4
                elif len_adj_arr > 4:
                    return 6
            elif len_adj_arr >= 2:
                return 5
            elif len_adj_arr == 0:
                return 1
        return 0

    def check_group_intactness(self):
        dummy_image = Atoms()
        neighborlist = build_neighbor_list(self.last_image, bothways=True, self_interaction=False, skin=self.skin)
        dok_matrix = neighborlist.get_connectivity_matrix()
        if self.group == 'OH':
            for idx, an_atom in enumerate(self.last_image):
                if an_atom.symbol == 'H':
                    dok_adj = dok_matrix[idx]
                    adj_array = dok_adj.tocoo().col
                    if not len(adj_array) == 1:
                        return 7, None
                else:
                    dummy_image.append(an_atom)

        if self.group == 'CH3':
            for idx, an_atom in enumerate(self.last_image):
                if an_atom.symbol == 'H':
                    dok_adj = dok_matrix[idx]
                    adj_array = dok_adj.tocoo().col
                    if not len(adj_array) == 1:
                        return 7, None
                elif idx > self.cage.size - 1:
                    an_atom.symbol = 'H'
                    dummy_image.append(an_atom)
                else:
                    dummy_image.append(an_atom)

        if self.group == 'CF3':
            for idx, an_atom in enumerate(self.last_image):
                if an_atom.symbol == 'F':
                    dok_adj = dok_matrix[idx]
                    adj_array = dok_adj.tocoo().col
                    if not len(adj_array) == 1:
                        return 7, None
                elif idx > self.cage.size - 1:
                    an_atom.symbol = 'H'
                    dummy_image.append(an_atom)
                else:
                    dummy_image.append(an_atom)

        return None, dummy_image



    def check_last_image(self):
        if self.group in ['CH3', 'CF3', 'OH']:
            group_status, self.last_image = self.check_group_intactness()
            if group_status:
                return None, group_status
        if self.group in ['OH', 'Cl', 'F', 'H', 'Br']:
            cutoffs = natural_cutoffs(self.last_image)
        else:
            cutoffs = [natural_cutoffs(self.cage.atoms)[0]] * len(self.last_image)
        neighborlist = NeighborList(cutoffs, skin=self.skin, bothways=True, self_interaction=False)
        neighborlist.update(self.last_image)
        dok_matrix = neighborlist.get_connectivity_matrix()
        ADJ_status = self.check_ADJ(dok_matrix)
        if ADJ_status:
            return None, ADJ_status
        else:
            bonds = np.array([*dok_matrix.keys()])
            map_dict = dict()
            for i in range(bonds.shape[0]):
                bond = bonds[i]
                symbol_0 = self.last_image[bond[0]].symbol
                symbol_1 = self.last_image[bond[1]].symbol
                if not symbol_0 == symbol_1:
                    if bond[0] < bond[1]:
                        if symbol_0 == self.cage.symbol:
                            cage_element_idx = bond[0]
                            addon_idx = bond[1]
                        else:
                            cage_element_idx = bond[1]
                            addon_idx = bond[0]
                        map_dict.update({addon_idx: cage_element_idx})
            addon_set = set(map_dict.values())
            return addon_set, None

        # Caution!! opt_root should be abs path
    def check(self, opt_mood: str, opt_root: str, is_init: bool, init_cycle: int=None):
        # status_code == None means the topology stays unchanged during optimization.
        # status_code == 1 means at least one addon atom breaks the bond with the cage and becomes a radical during optimization.
        # status_code == 2 means at least one addon atom breaks the bond with the original addon site and changed to another during optimization.
        # status_code == 3 means at least one 3 membered carbon ring formed during optimization, which is extremely unstable. (carbon for instance, this can be changed for other elements.)
        # status_code == 4 means at least one carbon atom only has 2 neighbor carbons or less, which means the cage is broken.
        # status_code == 5 means at least one addon atom binds with 2 or more atoms.
        # status_code == 6 means at least one carbon atom has 5 or more neighbors, which means a small cluster is formed.
        # status_code == 7 means the inner intactness of at least one addon group has broken.
        yes_list = []
        failed_list = []
        wrong_list = []
        init_yes_list = []
        opt_root = os.path.abspath(opt_root)
        for a_folder in os.listdir(opt_root):
            if opt_mood == 'xtb':
                log_path = os.path.join(opt_root, a_folder, 'xtbopt.log')
                if not os.path.exists(log_path):
                    wrong_list.append(os.path.join(opt_root, a_folder, f'{a_folder}.xyz') + '\n')
                    continue
                if is_init:
                    nimages, self.last_image = read_xtb_log(log_path=log_path)
                else:
                    self.last_image = read(log_path, format='xyz')

            elif opt_mood == 'ase':
                xyz_path = os.path.join(opt_root, a_folder, 'opt.xyz')
                if not os.path.exists(xyz_path):
                    wrong_list.append(str(xyz_path) + '\n')
                    continue
                self.last_image = read(xyz_path)
                log_path = os.path.join(opt_root, a_folder, 'opt.log')
                if is_init:
                    with open(log_path, 'r') as file:
                        nimages = len(file.readlines()) - 2
            elif opt_mood == 'gaussian':
                log_path = os.path.join(opt_root, a_folder, 'gau.log')
                if not os.path.exists(log_path):
                    wrong_list.append(os.path.join(opt_root, a_folder, f'{a_folder}.xyz') + '\n')
                    continue
                traj = read_gaussian_out(fd=log_path) # The read_gaussian_out function is slightly changed to get all traj.
                self.last_image = traj[-1]
                if is_init:
                    nimages = len(traj)
                new_log_path = os.path.join(opt_root, a_folder, 'gau_simple.log')
                last_xyz_path = os.path.join(opt_root, a_folder, 'gau.xyz')
                if not os.path.exists(new_log_path):
                    write(filename=new_log_path, images=traj, format='xyz')
                    e = self.last_image.get_potential_energy()
                    write(filename=last_xyz_path, images=self.last_image, format='xyz', comment=f'Energy: {str(e)}')

            opted_addonset, status_code = self.check_last_image()
            _, init_addonset, bin_arr = name2seq(a_folder, cage_size=self.cage.size)
            if status_code:
                failed_list.append(f'{status_code}' + '          ' + str(log_path) + '\n')
                self.cage.failed_bin_arr = np.vstack((self.cage.failed_bin_arr, bin_arr))
            else:
                if init_addonset == opted_addonset:
                    if is_init == True and nimages < init_cycle:
                        init_yes_list.append(str(log_path) + '\n')
                    else:
                        yes_list.append(str(log_path) + '\n')

                elif len(init_addonset) == len(opted_addonset):
                    status_code = 2
                    failed_list.append(f'{status_code}' + '          ' + str(log_path) + '\n')
                    self.cage.failed_bin_arr = np.vstack((self.cage.failed_bin_arr, bin_arr))
                else:
                    status_code = 1
                    failed_list.append(f'{status_code}' + '          ' + str(log_path) + '\n')
                    self.cage.failed_bin_arr = np.vstack((self.cage.failed_bin_arr, bin_arr))

        if is_init:
            with open('init_yes_paths', 'a') as f:
                f.writelines(init_yes_list)
        with open('yes_paths', 'w') as f:
            f.writelines(yes_list)

        with open('failed_paths', 'a') as f:
            f.writelines(failed_list)

        with open('wrong_paths', 'w') as f:
            f.writelines(wrong_list)

