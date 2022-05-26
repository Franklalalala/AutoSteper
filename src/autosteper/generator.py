import os
import numpy as np
from numpy import pi
from ase.atoms import Atom, Atoms
from ase.io import read, write
from ase.neighborlist import build_neighbor_list
from autosteper.cage import seq2name, Cage
from autosteper.tools import rotate_around_axis


pre_def_geom_para = {
    'len_cage_Cl': 1.81,
    'len_cage_Br': 2.00,
    'len_cage_F': 1.37,
    'len_cage_H': 1.10,
    'len_cage_CH3': 1.54,
    'len_cage_CF3': 1.55,
    'len_cage_OH': 1.41,
    'len_CH3': 1.09,
    'len_CF3': 1.33,
    'len_OH': 0.97,
    'ang_CH3': 110,
    'ang_CF3': 112,
    'ang_OH': 108
}


class Generator():
    def __init__(self, gen_para: dict):
        self.gen_core_path = gen_para['gen_core_path']
        self.group = gen_para['group']
        self.skin = gen_para['skin']
        self.geom_mode = gen_para['geom_mode']
        if self.geom_mode == 'pre_defined':
            self.geom_para = pre_def_geom_para
        else:
            self.geom_para = gen_para['geom_para']


    def gen_seq(self, cage: Cage, mode: str, gen_out_path: str, prev_seq: str=None, random_num: int=None):
        commandline = f'{self.gen_core_path} --graph6str \"{cage.graph6str}\" --addnum \"{cage.add_num}\" -o \"{gen_out_path}\" -m \"{mode}\" '
        if mode == 'step':
            if os.name == 'posix':
                prev_seq = '_'.join(prev_seq.split())
            commandline = commandline + prev_seq
        elif mode == 'random':
            commandline = commandline + f' -r {random_num}'
        os.system(commandline)


    def build(self, is_first: bool, gen_out_path: str, dump_folder: str, cage: Cage, prev_xyz_path: str, calc=None,
              parent_name: str=None, parent_info: dict=None, prev_addon_set: set=None, pre_scan_map: dict=None):

        def _build_unit(old_cage: Atoms, old_coords: np.array):
            # after this loop, old cage will become new cage, but the name stays.
            for i in addon_set:
                ADJ_coords = []
                dok_adj = dok_matrix[i]
                adj_array = dok_adj.tocoo().col
                if len(adj_array) > 3:
                    dist = []
                    for ii in adj_array:
                        dist.append(np.linalg.norm(old_coords[i] - old_coords[ii]))
                    map_dict = dict(zip(dist, adj_array))
                    sorted_dist = sorted(map_dict.keys())
                    for ii in sorted_dist[0:3]:
                        ADJ_coords.append(old_coords[map_dict[ii]])
                else:
                    for ii in adj_array:
                        ADJ_coords.append(old_coords[ii])
                norm_vec_unnormalized = np.cross(ADJ_coords[0] - ADJ_coords[1], ADJ_coords[0] - ADJ_coords[2])
                normal_vec = norm_vec_unnormalized / np.linalg.norm(norm_vec_unnormalized)

                if (old_coords[i] - cage.centre) @ normal_vec > 0:
                    addon_coord_0 = old_coords[i] + normal_vec * self.geom_para[f'len_cage_{self.group}']
                else:
                    addon_coord_0 = old_coords[i] - normal_vec * self.geom_para[f'len_cage_{self.group}']

                if self.group in ['Cl', 'F', 'H', 'Br']:
                    new_atom_0 = Atom(symbol=self.group, position=addon_coord_0)
                    old_cage.append(new_atom_0)
                else:
                    new_atom_0 = Atom(symbol=self.group[0], position=addon_coord_0)
                    old_cage.append(new_atom_0)
                    addon_symbol_1 = self.group[1]
                    new_atom_1 = Atom(symbol=addon_symbol_1)
                    old_cage.append(new_atom_1)
                    old_cage.set_angle(angle=self.geom_para[f'ang_{self.group}'], a1=i, a2=-2, a3=-1)
                    old_cage.set_distance(distance=self.geom_para[f'len_{self.group}'], a0=-2, a1=-1, fix=0)
                    if self.group in ['CF3', 'CH3']:
                        vec_new_atom_0_1 = old_cage[-1].position - old_cage[-2].position
                        addon_coord_2 = old_cage[-2].position + rotate_around_axis(norm_axis=normal_vec,
                                                                                old_vec=vec_new_atom_0_1,
                                                                                theta=2 / 3 * pi)
                        addon_coord_3 = old_cage[-2].position + rotate_around_axis(norm_axis=normal_vec,
                                                                                old_vec=vec_new_atom_0_1,
                                                                                theta=4 / 3 * pi)
                        new_atom_2 = Atom(symbol=self.group[1], position=addon_coord_2)
                        new_atom_3 = Atom(symbol=self.group[1], position=addon_coord_3)
                        old_cage.append(new_atom_2)
                        old_cage.append(new_atom_3)
            if not calc == None:
                old_cage.calc = calc
                # different api may have different ways to get e
                e = old_cage.get_total_energy()
                pre_scan_map.update({e: name})
                write(filename=os.path.join(dump_folder, f'{name}.xyz'), format='xyz', images=old_cage, comment=str(e))
            else:
                write(filename=os.path.join(dump_folder, f'{name}.xyz'), format='xyz', images=old_cage)


        old_cage = read(prev_xyz_path, format='xyz')
        old_coords = old_cage.get_positions()
        neighborList = build_neighbor_list(old_cage, bothways=True, self_interaction=False)
        dok_matrix = neighborList.get_connectivity_matrix()

        if is_first:
            with open(gen_out_path, 'r') as file:
                for a_seq in file.readlines():
                    name, addon_set, _ = seq2name(seq=a_seq, cage=cage)
                    _build_unit(old_cage=old_cage.copy(), old_coords=old_coords.copy())
            return pre_scan_map
        # Do not have a black list or blak list is empty
        elif cage.has_blk_list == False or len(cage.blk_list.blk_list_arr.shape) == 1:
            with open(gen_out_path, 'r') as file:
                for a_seq in file.readlines():
                    name, new_addon_set, a_bin_arr = seq2name(seq=a_seq, cage=cage)
                    if name in parent_info.keys():
                        parent_info[name][0].append(parent_name)
                    else:
                        addon_set = new_addon_set - prev_addon_set
                        _build_unit(old_cage=old_cage.copy(), old_coords=old_coords.copy())
                        parent_info.update({name: [[parent_name]]})
            return parent_info, pre_scan_map
        # Have a black list
        else:
            addon_set_list = []
            q_bin_arr_list = []
            name_list = []
            with open(gen_out_path, 'r') as file:
                for a_seq in file.readlines():
                    name, new_addon_set, a_bin_arr = seq2name(seq=a_seq, cage=cage)
                    if name in parent_info.keys():
                        parent_info[name][0].append(parent_name)
                    else:
                        addon_set_list.append(new_addon_set - prev_addon_set)
                        q_bin_arr_list.append(a_bin_arr)
                        name_list.append(name)
                        parent_info.update({name: [[parent_name]]})
            # Check black list in a batch
            if len(addon_set_list) > 0:
                uncutted_idxes = cage.blk_list.chk_blk(q_bin_arr_list=q_bin_arr_list)
                for an_idx in uncutted_idxes:
                    addon_set = addon_set_list[an_idx]
                    name = name_list[an_idx]
                    _build_unit(old_cage=old_cage.copy(), old_coords=old_coords.copy())
            return parent_info, pre_scan_map

