import os
from typing import Union

import numpy as np
from ase.atoms import Atom
from ase.io import read, write
from ase.neighborlist import build_neighbor_list
from autosteper.cage import seq2name, Cage
from autosteper.tools import rotate_around_axis
from numpy import pi

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
    def __init__(self, gen_para: dict, pst_cage: Cage):
        if 'gen_core_path' in gen_para.keys():
            self.gen_core_path = gen_para['gen_core_path']
        self.group = gen_para['group']
        self.geom_mode = gen_para['geom_mode']
        self.pst_cage = pst_cage
        self.ps_map = dict()
        if self.geom_mode == 'pre_defined':
            self.geom_para = pre_def_geom_para
            self.skin = 0.15
        else:
            self.geom_para = gen_para['geom_para']
            self.skin = gen_para['geom_para']['skin']

    def gen_seq(self, mode: str, gen_out_path: str, prev_seq: str = None, random_num: int = None):
        commandline = f'{self.gen_core_path} --graph6str \"{self.pst_cage.graph6str}\" --addnum \"{self.pst_cage.add_num}\" -o \"{gen_out_path}\" -m \"{mode}\" '
        if mode == 'step':
            if os.name == 'posix':
                prev_seq = '_'.join(prev_seq.split())
            commandline = commandline + prev_seq
        elif mode == 'random':
            commandline = commandline + f' -r {random_num}'
        os.system(commandline)

    def pre_build_unit(self, prev_xyz_path: str):
        self.old_cage = read(prev_xyz_path, format='xyz')
        self.old_coords = self.old_cage.get_positions()
        neighborList = build_neighbor_list(self.old_cage, bothways=True, self_interaction=False)
        self.dok_matrix = neighborList.get_connectivity_matrix()

    # new_addon_sites start from 0
    def build_unit(self, new_addon_sites: Union[list, set]):
        new_cage = self.old_cage.copy()
        for i in new_addon_sites:
            ADJ_coords = []
            dok_adj = self.dok_matrix[i]
            adj_array = dok_adj.tocoo().col
            if len(adj_array) > 3:
                dist = []
                for ii in adj_array:
                    dist.append(np.linalg.norm(self.old_coords[i] - self.old_coords[ii]))
                map_dict = dict(zip(dist, adj_array))
                sorted_dist = sorted(map_dict.keys())
                for ii in sorted_dist[0:3]:
                    ADJ_coords.append(self.old_coords[map_dict[ii]])
            else:
                for ii in adj_array:
                    ADJ_coords.append(self.old_coords[ii])
            norm_vec_unnormalized = np.cross(ADJ_coords[0] - ADJ_coords[1], ADJ_coords[0] - ADJ_coords[2])
            normal_vec = norm_vec_unnormalized / np.linalg.norm(norm_vec_unnormalized)

            if (self.old_coords[i] - self.pst_cage.centre) @ normal_vec > 0:
                addon_coord_0 = self.old_coords[i] + normal_vec * self.geom_para[f'len_cage_{self.group}']
            else:
                addon_coord_0 = self.old_coords[i] - normal_vec * self.geom_para[f'len_cage_{self.group}']

            if self.group in ['Cl', 'F', 'H', 'Br']:
                new_atom_0 = Atom(symbol=self.group, position=addon_coord_0)
                new_cage.append(new_atom_0)
            else:
                new_atom_0 = Atom(symbol=self.group[0], position=addon_coord_0)
                new_cage.append(new_atom_0)
                addon_symbol_1 = self.group[1]
                new_atom_1 = Atom(symbol=addon_symbol_1)
                new_cage.append(new_atom_1)
                new_cage.set_angle(angle=self.geom_para[f'ang_{self.group}'], a1=i, a2=-2, a3=-1)
                new_cage.set_distance(distance=self.geom_para[f'len_{self.group}'], a0=-2, a1=-1, fix=0)
                if self.group in ['CF3', 'CH3']:
                    vec_new_atom_0_1 = new_cage[-1].position - new_cage[-2].position
                    addon_coord_2 = new_cage[-2].position + rotate_around_axis(norm_axis=normal_vec,
                                                                               old_vec=vec_new_atom_0_1,
                                                                               theta=2 / 3 * pi)
                    addon_coord_3 = new_cage[-2].position + rotate_around_axis(norm_axis=normal_vec,
                                                                               old_vec=vec_new_atom_0_1,
                                                                               theta=4 / 3 * pi)
                    new_atom_2 = Atom(symbol=self.group[1], position=addon_coord_2)
                    new_atom_3 = Atom(symbol=self.group[1], position=addon_coord_3)
                    new_cage.append(new_atom_2)
                    new_cage.append(new_atom_3)
        return new_cage

    def build(self, is_first: bool, gen_out_path: str, dyn_cage: Cage,
              dump_folder: str, prev_xyz_path: str,
              parent_name: str = None, parent_info: dict = None, prev_addon_set: set = None):
        print('start build:\n')

        def _lazy_dump():
            if dyn_cage.is_pre_scan:
                new_cage.calc = dyn_cage.calc
                # different api may have different ways to get e
                e = new_cage.get_total_energy()[0][0]
                self.ps_map.update({name: e})
                write(filename=os.path.join(dump_folder, f'{name}.xyz'), format='xyz', images=new_cage, comment=str(e))
            else:
                write(filename=os.path.join(dump_folder, f'{name}.xyz'), format='xyz', images=new_cage)

        self.pre_build_unit(prev_xyz_path=prev_xyz_path)
        if is_first:
            with open(gen_out_path, 'r') as file:
                for a_seq in file.readlines():
                    name, addon_set, _ = seq2name(seq=a_seq, cage=self.pst_cage)
                    new_cage = self.build_unit(new_addon_sites=addon_set)
                    _lazy_dump()

        # Do not have a black list or blak list is empty
        elif dyn_cage.has_blk_list == False or len(dyn_cage.blk_list.blk_list_arr.shape) == 1:
            with open(gen_out_path, 'r') as file:
                for a_seq in file.readlines():
                    name, new_addon_set, a_bin_arr = seq2name(seq=a_seq, cage=self.pst_cage)
                    if name in parent_info.keys():
                        parent_info[name][0].append(parent_name)
                    else:
                        addon_set = new_addon_set - prev_addon_set
                        new_cage = self.build_unit(new_addon_sites=addon_set)
                        _lazy_dump()
                        parent_info.update({name: [[parent_name]]})
            return parent_info
        # Have a black list
        else:
            print('blk start:\n')
            addon_set_list = []
            q_bin_arr_list = []
            name_list = []
            with open(gen_out_path, 'r') as file:
                for a_seq in file.readlines():
                    name, new_addon_set, a_bin_arr = seq2name(seq=a_seq, cage=self.pst_cage)
                    if name in parent_info.keys():
                        parent_info[name][0].append(parent_name)
                    else:
                        addon_set_list.append(new_addon_set - prev_addon_set)
                        q_bin_arr_list.append(a_bin_arr)
                        name_list.append(name)
                        parent_info.update({name: [[parent_name]]})
            # Check black list in a batch
            if len(addon_set_list) > 0:
                uncutted_idxes = dyn_cage.blk_list.chk_blk(q_bin_arr_list=q_bin_arr_list)
                print(len(addon_set_list) - len(uncutted_idxes))
                for an_idx in uncutted_idxes:
                    addon_set = addon_set_list[an_idx]
                    name = name_list[an_idx]
                    new_cage = self.build_unit(new_addon_sites=addon_set)
                    _lazy_dump()
            return parent_info
