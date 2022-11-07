import os
import shutil

import networkx as nx
import numpy as np
import pandas as pd
from ase.atoms import Atoms, Atom
from ase.io import write, read
from ase.io.gaussian import read_gaussian_out
from ase.neighborlist import build_neighbor_list
from ase.units import Hartree
from networkx import isomorphism
from numpy import sin, cos


def rotate_around_axis(norm_axis: np.array, old_vec: np.array, theta: float):
    a, b, c = norm_axis
    rotate_mat = np.array([[(b ** 2 + c ** 2) * cos(theta) + a ** 2, a * b * (1 - cos(theta)) - c * sin(theta),
                            a * c * (1 - cos(theta)) + b * sin(theta)],
                           [a * b * (1 - cos(theta)) + c * sin(theta), b ** 2 + (1 - b ** 2) * cos(theta),
                            b * c * (1 - cos(theta)) - a * sin(theta)],
                           [a * c * (1 - cos(theta)) - b * sin(theta), b * c * (1 - cos(theta)) + a * sin(theta),
                            c ** 2 + (1 - c ** 2) * cos(theta)]])
    new_vec = old_vec @ rotate_mat
    return new_vec


def read_xtb_log(log_path):
    with open(log_path, "r") as file:
        lines = file.readlines()
        natoms = int(lines[0])
        nimages = len(lines) // (natoms + 2)
        n = (nimages - 1) * (natoms + 2) + 2
        symbols = []
        positions = []
        for line in lines[n:n + natoms]:
            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
        last_image = Atoms(symbols=symbols, positions=positions)
    return nimages, last_image


def get_last_from_log(log_path: str):
    """
    read xyz format logs(xtb opt logs)
    return last image and its energy
    """
    with open(log_path, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    nimages = len(lines) / (natoms + 2)
    final_cursor = int((nimages - 1) * (natoms + 2))
    final_e_line = lines[final_cursor + 1]
    e = float(final_e_line.split()[1]) * Hartree
    last_image = Atoms()
    for a_line in lines[final_cursor + 2:]:
        sym, x, y, z = a_line.split()
        an_atom = Atom(symbol=sym, position=[float(x), float(y), float(z)])
        last_image.append(an_atom)
    return e, last_image


def sort_atomic(atoms: Atoms):
    cage = Atoms()
    addends = Atoms()
    for an_atom in atoms:
        if an_atom.symbol == 'C':
            cage.append(an_atom)
        else:
            addends.append(an_atom)
    pristine_cage = cage.copy()
    cage.extend(addends)
    return pristine_cage, cage


# CAUTION!!! Pristine cage needs to be in the top position of coordinates.
def deal_group(atoms: Atoms, group: str, cage_size: int):
    new_atoms = Atoms()
    pristine_cage = atoms[:cage_size]
    if group == 'OH':
        for an_atom in atoms:
            if an_atom.symbol != 'H':
                new_atoms.append(an_atom)
    else:
        new_atoms = pristine_cage.copy()
        for an_atom in atoms[cage_size:]:
            if an_atom.symbol == 'C':
                an_atom.symbol = 'H'
                new_atoms.append(an_atom)
    return pristine_cage, new_atoms


def strip_extraFullerene(atoms: Atoms, group: str = None, cage_size: int = None):
    if group:
        pristine_cage, atoms = deal_group(atoms=atoms, group=group, cage_size=cage_size)
    else:
        pristine_cage, atoms = sort_atomic(atoms=atoms)
    if group in ['CH3', 'CF3']:
        cutoffs = [natural_cutoffs(pristine_cage)[0]] * len(atoms)
        neighborlist = NeighborList(cutoffs, skin=0.15, bothways=True, self_interaction=False)
    else:
        neighborlist = build_neighbor_list(atoms=atoms, skin=0.15, bothways=True, self_interaction=False)
    dok_matrix = neighborlist.get_connectivity_matrix()
    bonds = np.array([*dok_matrix.keys()])
    map_dict = dict()
    for i in range(bonds.shape[0]):
        bond = bonds[i]
        symbol_0 = atoms[bond[0]].symbol
        symbol_1 = atoms[bond[1]].symbol
        if not symbol_0 == symbol_1:
            if bond[0] < bond[1]:
                if symbol_0 == 'C':
                    cage_element_idx = bond[0]
                    addon_idx = bond[1]
                else:
                    cage_element_idx = bond[1]
                    addon_idx = bond[0]
                map_dict.update({addon_idx: cage_element_idx})
    addon_set = set(map_dict.values())
    return pristine_cage, addon_set


# e_arr needs to be sorted
def get_low_e_ranks(e_arr: np.ndarray, para: dict, is_reverse: bool = False):
    if para == None:
        para = {'mode': 'None'}
    assert para['mode'] in ['None', 'rank', 'value',
                            'value_and_rank', 'rank_and_value',
                            'value_or_rank', 'rank_or_value'], f'Please check your run cutoff mode keyword.'
    assert isinstance(e_arr, np.ndarray) == True

    if para['mode'] == 'None':
        rank_list = range(len(e_arr))
    elif is_reverse:
        if para['mode'] == 'rank':
            rank_list = range(len(e_arr))[::-1][:para['rank']]
        else:
            e_arr = e_arr - e_arr[0]
            if e_arr[-1] < para['value']:
                idx = len(e_arr)
            else:
                for idx, a_e in enumerate(e_arr[::-1]):
                    if a_e <= e_arr[-1] - para['value']:
                        break
            if para['mode'] == 'value_and_rank' or para['mode'] == 'rank_and_value':
                idx = min(idx, para['rank'])
            if para['mode'] == 'value_or_rank' or para['mode'] == 'rank_or_value':
                idx = max(idx, para['rank'])
            rank_list = range(len(e_arr))[::-1][:idx]
    else:
        if para['mode'] == 'rank':
            rank_list = range(len(e_arr))[:para['rank']]
        else:
            e_arr = e_arr - min(e_arr)
            if e_arr[-1] < para['value']:
                idx = len(e_arr)
            else:
                for idx, a_e in enumerate(e_arr):
                    if a_e > para['value']:
                        break

            if para['mode'] == 'value_and_rank' or para['mode'] == 'rank_and_value':
                idx = min(idx, para['rank'])
            if para['mode'] == 'value_or_rank' or para['mode'] == 'rank_or_value':
                idx = max(idx, para['rank'])
            rank_list = range(len(e_arr))[:idx]

    if 'nimg_th' in para.keys():
        new_rank_list = []
        for idx, a_nimg in enumerate(para['nimages']):
            if a_nimg >= para['nimg_th']:
                new_rank_list.append(rank_list[idx])
        rank_list = new_rank_list
    for a_rank in rank_list:
        yield a_rank


def get_low_e_xyz(dump_folder: str, add_num: int, old_workbase: str, cutoff_para: dict):
    dump_folder = os.path.abspath(dump_folder)
    cwd_ = os.getcwd()
    addon_path = str(add_num) + 'addons'
    os.chdir(old_workbase)
    os.chdir(addon_path)
    passed_info = pd.read_pickle('passed_info.pickle')
    e_arr = np.array(passed_info['energy'])
    xyz_paths = passed_info['xyz_path']

    os.chdir(dump_folder)
    for a_rank in get_low_e_ranks(e_arr=e_arr, para=cutoff_para):
        a_xyz_path = xyz_paths[a_rank]
        shutil.copy(src=a_xyz_path, dst=f'{add_num}_addons_rank_{a_rank + 1}.xyz')
    os.chdir(cwd_)


def get_G(atoms: Atoms):
    nb_list = build_neighbor_list(atoms, self_interaction=False)
    mat = nb_list.get_connectivity_matrix()
    adj = mat.toarray()
    G = nx.from_numpy_array(adj)
    assert isinstance(G, object)
    return G


def lazy_file_mid_name(file_mid_name: str = None):
    if file_mid_name != None:
        file_mid_name = '_' + file_mid_name + '_'
    else:
        file_mid_name = '_'
    return file_mid_name


def simple_parse_logs(dump_root: str, src_root: str, mode: str,
                      group: str = None, cage_size: int = None,
                      file_mid_name: str = None):
    """
    Turn disordered output logs to standard infos, xyz and logs.

    read_gaussian_out is edited to get all the trajectory atoms
    dump_root: path to dump
    src_root: path of the original root
    mode: 'gauss' for gaussian format log, 'xyz' for xyz format log.
    group: in case the group have more than one atom
    cage_size: combined usage for strip_extraFullerene
    """
    cwd_ = os.getcwd()
    os.makedirs(dump_root, exist_ok=True)
    os.chdir(dump_root)
    os.makedirs('xyz', exist_ok=True)
    os.makedirs('log', exist_ok=True)
    os.makedirs('info', exist_ok=True)
    dump_xyz = os.path.abspath('xyz')
    dump_log = os.path.abspath('log')
    dump_info = os.path.abspath('info')
    e_map = {}
    e_name_map = {}
    e_atom_map = {}
    add_num_list = []
    os.chdir(src_root)
    for a_file in os.listdir('./'):
        if mode == 'gauss':
            try:
                atoms = read_gaussian_out(fd=a_file)[-1]
                e = atoms.get_total_energy()
            except:
                # for the errored optimization file
                atoms = read_gaussian_out(fd=a_file)[-2]
                e = atoms.get_total_energy()
        elif mode == 'xyz':
            try:
                e, atoms = get_last_from_log(log_path=a_file)
            except Exception as e:
                raise RuntimeError(f'Error: {str(e)}\nPlease check the xyz format log file.')

        _, addon_set = strip_extraFullerene(atoms=atoms, group=group, cage_size=cage_size)
        add_num = len(addon_set)

        if add_num not in e_map.keys():
            add_num_list.append(add_num)
            dummy_list = []
            dummy_list.append(e)
            e_map.update({add_num: dummy_list})
        else:
            old_dummy_list = e_map[add_num]
            old_dummy_list.append(e)
            e_map.update({add_num: old_dummy_list})
        e_name_map.update({e: a_file})
        e_atom_map.update({e: atoms})
    max_len = 0
    for a_key in add_num_list:
        e_list = e_map[a_key]
        if max_len < len(e_list):
            max_len = len(e_list)

    file_mid_name = lazy_file_mid_name(file_mid_name)

    for a_key in add_num_list:
        e_list = e_map[a_key]
        e_list = sorted(e_list)
        for idx, a_e in enumerate(e_list):
            write(filename=os.path.join(dump_xyz, f'{str(a_key)}_addons{file_mid_name}{str(idx + 1)}.xyz'),
                  images=e_atom_map[a_e], format='xyz')
            shutil.copy(src=e_name_map[a_e],
                        dst=os.path.join(dump_log, f'{str(a_key)}_addons{file_mid_name}{str(idx + 1)}.log'))
        for i in range(max_len - len(e_list)):
            e_list.append(None)
        e_map.update({a_key: e_list})

    info = pd.DataFrame(columns=sorted(add_num_list))
    for an_add in info.columns:
        info[an_add] = e_map[an_add]
    info.to_pickle(path=os.path.join(dump_info, 'info.pickle'))
    info.to_excel(os.path.join(dump_info, 'info.xlsx'))
    os.chdir(cwd_)


def get_connection(xyz_root: str, connection_dump: str, step: int, file_mid_name: str = None):
    """
    Get connection infos for simple_parse_logs function outputs.
    """
    xyz_root = os.path.abspath(xyz_root)
    connection_dump = os.path.abspath(connection_dump)
    cwd_ = os.getcwd()
    os.makedirs(connection_dump, exist_ok=True)
    file_mid_name = lazy_file_mid_name(file_mid_name)

    os.chdir(xyz_root)
    max_rank = 1
    max_add_num = 1
    for a_file in os.listdir('./'):
        a_name = a_file.split('.')[0]
        an_add_num = int(a_name.split('_')[0])
        a_rank = int(a_name.split('_')[-1])
        if a_rank > max_rank:
            max_rank = a_rank
        if an_add_num > max_add_num:
            max_add_num = an_add_num

    prev_G_list = []
    for a_rank in range(max_rank):
        file_name = f'{max_add_num}_addons{file_mid_name}{a_rank + 1}.xyz'
        if not os.path.exists(file_name):
            continue
        a_prev_atoms = read(file_name)
        prev_G_list.append(get_G(a_prev_atoms))

    for an_add in np.arange(max_add_num - step, -step, -step):
        new_G_list = []
        for a_rank in range(max_rank):
            cnt_list = []
            file_name = f'{an_add}_addons{file_mid_name}{a_rank + 1}.xyz'
            if not os.path.exists(file_name):
                continue
            new_atoms = read(file_name)
            new_G = get_G(new_atoms)
            new_G_list.append(new_G)
            for prev_G in prev_G_list:
                GM = isomorphism.GraphMatcher(prev_G, new_G)
                if GM.subgraph_is_isomorphic():
                    cnt_list.append(1)
                else:
                    cnt_list.append(0)
            np.save(file=os.path.join(connection_dump, f'{an_add}_addons{file_mid_name}{a_rank + 1}.npy'),
                    arr=np.array(cnt_list))
        prev_G_list = new_G_list
    os.chdir(cwd_)


def simple_dump_pathway(pathway_info_path: str, dump_root: str, src_xyz_root: str, top_K: int = None):
    cwd_ = os.getcwd()
    pathway_info = pd.read_pickle(pathway_info_path)
    dump_root = os.path.abspath(dump_root)
    os.makedirs(dump_root, exist_ok=True)
    src_xyz_root = os.path.abspath(src_xyz_root)
    for idx, a_pathway in enumerate(pathway_info['name'][:top_K]):
        a_sub_dump = os.path.join(dump_root, f'path_rank_{idx + 1}')
        os.makedirs(a_sub_dump)
        os.chdir(src_xyz_root)
        for a_name in a_pathway:
            shutil.copy(src=a_name, dst=os.path.join(a_sub_dump, a_name))
        os.chdir(a_sub_dump)
        for a_name in a_pathway:
            with open('traj.log', 'a') as f_w, open(a_name, 'r') as f_r:
                for a_line in f_r.readlines():
                    f_w.write(a_line)
    os.chdir(cwd_)


def get_pathway_info(e_info_path: str, xyz_root: str, cnt_root: str, dump_info_path: str, file_mid_name: str = None):
    cwd_ = os.getcwd()

    def _get_path_way_unit(idx: int, name_list: list, rel_e_list: list):
        if idx == -1:
            all_name_list.append(name_list)
            all_rel_e_list.append(rel_e_list)
            all_e_area_list.append(sum(rel_e_list))
            return
        else:
            an_add = add_num_list[idx]
            a_min_e = min_e_list[idx]
            for a_rank in range(max_rank):
                file_name = f'{an_add}_addons{file_mid_name}{a_rank + 1}.xyz'
                if not os.path.exists(file_name):
                    continue
                cnt = np.load(os.path.join(cnt_root, f'{an_add}_addons{file_mid_name}{a_rank + 1}.npy'))
                prev_rank = int(name_list[0].split('_')[-1][:-4]) - 1
                if not len(cnt) > prev_rank:
                    continue
                if cnt[prev_rank] == 1:
                    name_list = [file_name, *name_list]
                    rel = e_info[an_add][a_rank] - a_min_e
                    rel_e_list = [rel, *rel_e_list]
                    _get_path_way_unit(idx=idx - 1, name_list=name_list.copy(), rel_e_list=rel_e_list.copy())
                    del rel_e_list[0]
                    del name_list[0]

    file_mid_name = lazy_file_mid_name(file_mid_name)
    cnt_root = os.path.abspath(cnt_root)
    dump_info_path = os.path.abspath(dump_info_path)
    os.chdir(xyz_root)
    e_info = pd.read_pickle(e_info_path)
    add_num_list = sorted(list(e_info.columns))
    min_e_list = []
    max_rank = 0
    for a_max_add in add_num_list:
        min_e_list.append(e_info[a_max_add][0])
        max_rank = max(max_rank, len(e_info[a_max_add]))
    all_name_list = []
    all_rel_e_list = []
    all_e_area_list = []
    for idx, a_max_add_e in enumerate(e_info[a_max_add]):
        a_max_rel_e = a_max_add_e - min_e_list[-1]
        a_max_name = f'{a_max_add}_addons{file_mid_name}{idx + 1}.xyz'
        _get_path_way_unit(idx=len(add_num_list) - 2, name_list=[a_max_name], rel_e_list=[a_max_rel_e])
    all_e_area_list = np.array(all_e_area_list) - min(all_e_area_list)
    info = pd.DataFrame({'name': all_name_list, 'rel_e': all_rel_e_list, 'e_area': all_e_area_list})
    info = info.sort_values(by='e_area')
    info.index = sorted(info.index)
    info.to_pickle(path=dump_info_path)
    os.chdir(cwd_)


def map_SWR(q_atoms: Atoms, tgt_atoms: Atoms):
    """
    Currently only support single SWR scenario.
    q_atoms: atoms that to be queried
    tgt_atoms: the target atoms
    swr_container: List that contains the SWR positions,
    first element corresponds to the queried graph, second for the target graph.
    """
    q_G = get_G(q_atoms)
    tgt_G = get_G(tgt_atoms)
    swr_container = []
    for an_edge in q_G.edges:
        i, j = an_edge
        dummy_q_G = q_G.copy()
        dummy_q_G.remove_node(i)
        dummy_q_G.remove_node(j)
        GM = isomorphism.GraphMatcher(tgt_G, dummy_q_G)
        if GM.subgraph_is_isomorphic():
            swr_container.append([i, j])
            for an_edge in tgt_G.edges:
                ii, jj = an_edge
                dummy_tgt_G = tgt_G.copy()
                dummy_tgt_G.remove_node(ii)
                dummy_tgt_G.remove_node(jj)
                if nx.is_isomorphic(dummy_tgt_G, dummy_q_G):
                    swr_container.append([ii, jj])
                    return swr_container
