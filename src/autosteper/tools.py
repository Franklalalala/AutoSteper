import os
import shutil
from typing import Union

import networkx as nx
import numpy as np
import pandas as pd
from ase.atoms import Atoms, Atom
from ase.io import write, read
from ase.io.gaussian import read_gaussian_out
from ase.neighborlist import build_neighbor_list
from ase.units import Hartree
from networkx import isomorphism, Graph
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


def sort_complex_group(atoms: Atoms):
    sorted_atoms = Atoms()
    cage_idx_list = []
    old_addons_list = []
    neighborList = build_neighbor_list(atoms, bothways=True, self_interaction=False)
    dok_matrix = neighborList.get_connectivity_matrix()
    for idx, an_atom in enumerate(atoms):
        if not an_atom.symbol == 'C':
            continue
        flag = 0
        dok_adj = dok_matrix[idx]
        adj_array = dok_adj.tocoo().col
        for a_ngr in adj_array:
            if atoms[a_ngr].symbol != 'C':
                flag = 1
            else:
                carbon_ngr_idx = a_ngr
        if flag == 0:
            cage_idx_list.append(idx)
        else:
            old_addons_list.append(carbon_ngr_idx)
    new_addons_set = set()
    cage_size_count = 0
    for idx, an_atom in enumerate(atoms):
        if idx in old_addons_list:
            new_addons_set.add(cage_size_count)
        if idx in cage_idx_list:
            sorted_atoms.append(an_atom)
            cage_size_count = cage_size_count + 1
    for idx, an_atom in enumerate(atoms):
        if not idx in cage_idx_list:
            sorted_atoms.append(an_atom)
    return sorted_atoms, cage_size_count, new_addons_set


def sort_simple_group(atoms: Atoms):
    sorted_atoms = Atoms()
    new_addons_set = set()
    neighborList = build_neighbor_list(atoms, bothways=True, self_interaction=False)
    dok_matrix = neighborList.get_connectivity_matrix()
    cage_size_count = 0
    for idx, an_atom in enumerate(atoms):
        if not an_atom.symbol == 'C':
            continue
        sorted_atoms.append(an_atom)
        dok_adj = dok_matrix[idx]
        adj_array = dok_adj.tocoo().col
        for a_ngr in adj_array:
            if atoms[a_ngr].symbol != 'C':
                new_addons_set.add(cage_size_count)
        cage_size_count = cage_size_count + 1
    for an_atom in atoms:
        if an_atom.symbol != 'C':
            sorted_atoms.append(an_atom)
    return sorted_atoms, cage_size_count, new_addons_set


def strip_extraFullerene(atoms: Atoms = None, coord_file_path: str = None, group: str = None):
    if coord_file_path:
        atoms = read(coord_file_path)
    elif not atoms:
        raise RuntimeError('Please input geometry-containing file, or input an ASE Atoms object.')

    if group in ['CH3', 'CF3']:
        sorted_atoms, cage_size, addon_set = sort_complex_group(atoms=atoms)
        return sorted_atoms[:cage_size], addon_set
    else:
        sorted_atoms, cage_size, addon_set = sort_simple_group(atoms=atoms)
        return sorted_atoms[:cage_size], addon_set


# e_arr needs to be sorted
def get_low_e_ranks(e_arr: np.ndarray, para: dict, is_reverse: bool = False):
    if para == None:
        para = {'mode': 'None'}
    assert para['mode'] in ['None', 'rank', 'value',
                            'value_and_rank', 'rank_and_value',
                            'value_or_rank', 'rank_or_value'], f'Please check your run cutoff mode keyword.'
    if not isinstance(e_arr, np.ndarray):
        e_arr = np.array(e_arr)

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


def rm_addons_from_G(G: Graph, rm_addon_sites: Union[list, set]):
    for a_rm_addon_site in rm_addon_sites:
        for a_ngb in G[a_rm_addon_site]:
            if len(G[a_ngb]) < 3:
                G.remove_node(a_ngb)
                break
    return G


def rm_addons_from_atoms(src_atoms: Atoms, rm_site_list: Union[list, set], group: str=None):
    rm_atoms_idx_list = []
    neighborList = build_neighbor_list(src_atoms, bothways=True, self_interaction=False)
    dok_matrix = neighborList.get_connectivity_matrix()
    non_C_count = 0
    if group in ['CH3', 'CF3']:
        for an_atom in src_atoms:
            if an_atom.symbol != 'C':
                non_C_count = non_C_count + 1
        cage_size = len(src_atoms) - non_C_count / 3 * 4
    for a_site in rm_site_list:
        dok_adj = dok_matrix[a_site]
        adj_array = dok_adj.tocoo().col
        for a_ngr in adj_array:
            if group in ['CH3', 'CF3']:
                if a_ngr >= cage_size:
                    break
            else:
                if src_atoms[a_ngr].symbol != 'C':
                    break
        rm_atoms_idx_list.append(a_ngr)
        if group in ['OH', 'CH3', 'CF3']:
            dok_adj = dok_matrix[a_ngr]
            adj_array = dok_adj.tocoo().col
            for a_ngr_ngr in adj_array:
                if src_atoms[a_ngr_ngr].symbol != 'C':
                    rm_atoms_idx_list.append(a_ngr_ngr)
    removed_atoms = Atoms()
    for idx, an_atom in enumerate(src_atoms):
        if not idx in rm_atoms_idx_list:
            removed_atoms.append(an_atom)
    return removed_atoms


def get_a_seq(all_sites: Union[list, set], tgt_seq_len: int, tgt_seq: set = None):
    if type(all_sites) == list:
        all_sites = set(all_sites)
    if tgt_seq is None:
        tgt_seq = set()
    if len(tgt_seq) == tgt_seq_len:
        yield tgt_seq
    else:
        for i in all_sites - tgt_seq:
            tgt_seq.add(i)
            for a_seq in get_a_seq(all_sites=all_sites, tgt_seq_len=tgt_seq_len, tgt_seq=tgt_seq):
                if len(a_seq) == tgt_seq_len:
                    yield a_seq
            tgt_seq.remove(i)

def get_max_adduct_filename(root: str):
    max_add = 0
    max_filename = str()
    for a_file in os.listdir(root):
        if a_file.endswith('xyz'):
            if 'addons' in a_file.split('_'):
                an_add = int(a_file.split('_')[0])
                if max_add < an_add:
                    max_add = an_add
                    max_filename = a_file
    return max_filename


def re_label_a_pathway(src_root: str, re_label_root: str, group_symbol: str, has_header: bool = True):
    def _get_re_label_filename(addon_set, has_header):
        if has_header:
            re_label = f'{len(addon_set)}_addons_' + \
                       '_'.join([str(x + 1) for x in sorted(list(addon_set))]) + '.xyz'
        else:
            re_label = '_'.join([str(x + 1) for x in sorted(list(addon_set))]) + '.xyz'
        return re_label
    cwd_re_label = os.getcwd()
    re_label_root = os.path.abspath(re_label_root)
    add_list = []
    add_atoms_map = {}
    os.chdir(src_root)
    for a_file in os.listdir('./'):
        if a_file.endswith('xyz'):
            an_add = int(a_file.split('_')[0])
            add_list.append(an_add)
            add_atoms_map.update({an_add: read(a_file)})
    add_list.sort()
    max_adduct = add_atoms_map[add_list[-1]]
    pristine_cage, max_addon_set = strip_extraFullerene(atoms=max_adduct, group=group_symbol)
    os.chdir(re_label_root)
    write(filename=r'pristine_cage.xyz', images=pristine_cage, format='xyz',
          comment=r'This file is for visualization only.')
    re_label = _get_re_label_filename(addon_set=max_addon_set, has_header=has_header)
    write(filename=re_label, images=max_adduct, format='xyz')
    prev_G = get_G(max_adduct)
    prev_add_set = max_addon_set
    prev_atoms = max_adduct.copy()
    for nxt_add in add_list[-2::-1]:
        nxt_atoms = add_atoms_map[nxt_add]
        nxt_G = get_G(nxt_atoms)
        _, nxt_add_set = strip_extraFullerene(atoms=nxt_atoms, group=group_symbol)
        diff_add_set = prev_add_set - nxt_add_set
        diff_len = len(prev_add_set) - len(nxt_add_set)
        # in case the addon set matches directly
        if len(diff_add_set) == diff_len:
            prev_add_set = nxt_add_set
            prev_G = nxt_G
            prev_atoms = rm_addons_from_atoms(src_atoms=prev_atoms, rm_site_list=diff_add_set, group=group_symbol)
            re_label = _get_re_label_filename(addon_set=nxt_add_set, has_header=has_header)
            write(filename=re_label, images=prev_atoms.copy(), format='xyz',
                  comment=r'This file is for visualization only.')
            add_atoms_map.update({len(prev_add_set): prev_atoms.copy()})
        else:  # in case the addon set does not match directly
            for a_seq in get_a_seq(all_sites=prev_add_set, tgt_seq_len=diff_len):
                removed_G = rm_addons_from_G(G=prev_G.copy(), rm_addon_sites=a_seq)
                GM = isomorphism.GraphMatcher(removed_G, nxt_G)
                if GM.is_isomorphic():
                    prev_add_set = prev_add_set - a_seq
                    prev_G = removed_G
                    prev_atoms = rm_addons_from_atoms(src_atoms=prev_atoms, rm_site_list=a_seq, group=group_symbol)
                    re_label = _get_re_label_filename(addon_set=prev_add_set, has_header=has_header)
                    write(filename=re_label, images=prev_atoms.copy(), format='xyz',
                          comment=r'This file is for visualization only.')
                    add_atoms_map.update({len(prev_add_set): prev_atoms.copy()})
                    break
    write(filename='traj.log', images=pristine_cage, format='xyz',
          comment=r'This file is for visualization only.')
    for an_add in add_list:
        write(filename='traj.log', images=add_atoms_map[an_add], format='xyz', append=True)
    os.chdir(cwd_re_label)


def match_max_adduct(query_atoms_path: str, tgt_atoms_path: str, diff_len: int, group_symbol: str):
    query_atoms = read(query_atoms_path)
    tgt_atoms = read(tgt_atoms_path)
    query_G = get_G(query_atoms)
    tgt_G = get_G(tgt_atoms)
    # files need to be sorted first, carbon atoms of cage should stay on top of coordinates.
    _, tgt_addon_set = strip_extraFullerene(atoms=tgt_atoms, group=group_symbol)
    pristine_cage, query_addon_set = strip_extraFullerene(atoms=query_atoms, group=group_symbol)
    # ablation of target atoms to match the query atom
    for a_seq in get_a_seq(all_sites=tgt_addon_set, tgt_seq_len=diff_len):
        removed_G = rm_addons_from_G(G=tgt_G.copy(), rm_addon_sites=a_seq)
        GM = isomorphism.GraphMatcher(query_G, removed_G)
        if GM.subgraph_is_isomorphic():
            half_G = removed_G
            half_matched_atoms = rm_addons_from_atoms(src_atoms=tgt_atoms, rm_site_list=a_seq)
            break
    dummy_all_adduct_set = set(range(len(query_atoms)))
    dummy_all_cage_set = set(range(len(pristine_cage)))
    append_nodes_list = list(dummy_all_adduct_set - set(half_G.nodes))
    # re-build half-matched graph to match the query atom
    for a_seq in get_a_seq(all_sites=dummy_all_cage_set - tgt_addon_set, tgt_seq_len=diff_len):
        for idx, an_addon in enumerate(a_seq):
            half_G.add_node(append_nodes_list[idx])
            half_G.add_edge(an_addon, append_nodes_list[idx], weight=1)
        GM = isomorphism.GraphMatcher(query_G, half_G)
        if GM.is_isomorphic():
            match_map = next(GM.match())
            re_matched_atoms = Atoms()
            for a_tgt_idx in range(len(dummy_all_adduct_set)):
                for idx, a_query_atom in enumerate(query_atoms):
                    mapped_idx = match_map[idx]
                    if mapped_idx == a_tgt_idx:
                        re_matched_atoms.append(a_query_atom)
            break
        else:
            for an_addon in a_seq:
                half_G.remove_node(an_addon)
    return half_matched_atoms, re_matched_atoms


def lazy_file_mid_name(file_mid_name: str = None):
    if file_mid_name != None:
        file_mid_name = '_' + file_mid_name + '_'
    else:
        file_mid_name = '_'
    return file_mid_name


def simple_parse_logs(dump_root: str, src_root: str, mode: str,
                      group: str = None, file_mid_name: str = None):
    """
    Turn disordered output logs to standard infos, xyz and logs.

    read_gaussian_out is edited to get all the trajectory atoms
    dump_root: path to dump
    src_root: path of the original root
    mode: 'gauss' for gaussian format log, 'xyz' for xyz format log.
    group: in case the group have more than one atom
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

        _, addon_set = strip_extraFullerene(atoms=atoms, group=group)
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
        os.makedirs(a_sub_dump, exist_ok=True)
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
    info[:1000].to_excel(dump_info_path[:-6]+'xlsx')
    os.chdir(cwd_)


# To do: a better way to get swr maps
def map_SWR(q_atoms: Atoms, tgt_atoms: Atoms):
    """
    Currently only support single SWR scenario.
    q_atoms: atoms that to be queried
    tgt_atoms: the target atoms
    swr_container: List that contains the SWR positions,
    first element corresponds to the queried graph, second for the target graph.
    """
    Print('You have entered a SWR match program, this may take a few minutes, take a break for coffee please.')
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


def match_swr_cages(query_cage: Atoms, target_cage: Atoms):
    # swr = map_SWR(q_atoms=query_cage, tgt_atoms=target_cage)
    swr = [[29, 31], [54, 55]] # for test
    query_G = get_G(query_cage)
    dummy_query_G = query_G.copy()
    dummy_query_G.remove_nodes_from(swr[0])

    target_G = get_G(target_cage)
    dummy_target_G = target_G.copy()
    dummy_target_G.remove_nodes_from(swr[1])

    GM = isomorphism.GraphMatcher(dummy_query_G, dummy_target_G)
    match_map = next(GM.match())
    new_query_cage = Atoms()
    flag = 0
    for a_tgt_idx in range(len(query_cage)):
        if not a_tgt_idx in match_map.values():
            new_query_cage.append(query_cage[swr[0][flag]])
            flag = 1
        for a_q_idx, a_query_atom in enumerate(query_cage):
            if not a_q_idx in match_map.keys():
                continue
            mapped_idx = match_map[a_q_idx]
            if mapped_idx == a_tgt_idx:
                new_query_cage.append(a_query_atom)
    return new_query_cage, match_map, swr[1]
