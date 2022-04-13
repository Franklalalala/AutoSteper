import os
from autosteper.tools import read_xtb_log
from autosteper.cage import name2seq
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs, build_neighbor_list
from ase.atoms import Atoms
import numpy as np


class Checker():
    def __init__(self, chk_skin: float, group: str, cage_size: int):
        self.skin = chk_skin
        self.group = group
        self.cage_size = cage_size


    def check_ADJ(self, dok_matrix):
        for idx, an_atom in enumerate(self.last_image):
            dok_adj = dok_matrix[idx]
            adj_array = dok_adj.tocoo().col
            len_adj_arr = len(adj_array)
            if an_atom.symbol == 'C':
                if len_adj_arr == 4:
                    flag = 0
                    for ii in adj_array:
                        if self.last_image[ii].symbol == self.group or  self.last_image[ii].symbol == 'H':
                            flag = 1
                            break
                    if flag:
                        continue
                    else:

                        aa=1
                        return 3
                elif len_adj_arr < 3:
                    return 4
                elif len_adj_arr == 3:
                    for ii in adj_array:
                        if self.last_image[ii].symbol == self.group or self.last_image[ii].symbol == 'H':
                            return 4
            elif len_adj_arr == 2 and self.last_image[adj_array[0]].symbol == 'C' and self.last_image[adj_array[1]].symbol == 'C':
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
                        return 6, None
                else:
                    dummy_image.append(an_atom)

        if self.group == 'CH3':
            for idx, an_atom in enumerate(self.last_image):
                if an_atom.symbol == 'H':
                    dok_adj = dok_matrix[idx]
                    adj_array = dok_adj.tocoo().col
                    if not len(adj_array) == 1:
                        return 6, None
                elif idx > self.cage_size - 1:
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
                        return 6, None
                elif idx > self.cage_size - 1:
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
        if self.group in ['OH', 'Cl', 'F', 'H']:
            cutoffs = natural_cutoffs(self.last_image)
        else:
            cutoffs = [0.76] * len(self.last_image)
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
                        if symbol_0 == 'C':
                            carbon_idx = bond[0]
                            addon_idx = bond[1]
                        else:
                            carbon_idx = bond[1]
                            addon_idx = bond[0]
                        map_dict.update({addon_idx: carbon_idx})
            addon_set = set(map_dict.values())
            return addon_set, None

        # Caution!! opt_root should be abs path
    def check(self, opt_mood: str, opt_root: str, is_init: bool, init_cycle: int=None):
        # status_code == None means the topology stays unchanged during optimization.
        # status_code == 1 means at least one addon atom breaks the bond with the cage and becomes a radical during optimization.
        # status_code == 2 means at least one addon atom breaks the bond with the original addon site and changed to another during optimization.
        # status_code == 3 means at least one 3 membered carbon ring formed during optimization, which is extremely unstable.
        # status_code == 4 means at least one carbon atom has 2 neighbors, which means the cage is broken.
        # status_code == 5 means at least one addon atom binds with 2 carbon atoms.
        # status_code == 6 means the inner intactness of at least one addon group has broken.
        yes_list = []
        failed_list = []
        wrong_list = []
        init_yes_list = []
        print(opt_root)
        opt_root = os.path.abspath(opt_root)
        print(opt_root)
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
                if is_init:
                    log_path = os.path.join(opt_root, a_folder, 'opt.log')
                    with open(log_path, 'r') as file:
                        nimages = len(file.readlines()) - 2

            opted_addonset, status_code = self.check_last_image()
            _, init_addonset = name2seq(a_folder, cage_size=self.cage_size)
            if status_code:
                failed_list.append(f'{status_code}' + '          ' + str(log_path) + '\n')
            else:
                if init_addonset == opted_addonset:
                    if is_init == True and nimages < init_cycle:
                        init_yes_list.append(str(log_path) + '\n')
                    else:
                        yes_list.append(str(log_path) + '\n')

                elif len(init_addonset) == len(opted_addonset):
                    status_code = 2
                    failed_list.append(f'{status_code}' + '          ' + str(log_path) + '\n')
                else:
                    status_code = 1
                    failed_list.append(f'{status_code}' + '          ' + str(log_path) + '\n')

        if is_init:
            with open('init_yes_paths', 'a') as f:
                f.writelines(init_yes_list)
        with open('yes_paths', 'w') as f:
            f.writelines(yes_list)

        with open('failed_paths', 'a') as f:
            f.writelines(failed_list)

        with open('wrong_paths', 'w') as f:
            f.writelines(wrong_list)

