# To refactor!

import os
import shutil

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import tqdm
from ase import Atoms
from ase.io.gaussian import read_gaussian_out
from ase.neighborlist import build_neighbor_list
from ase.units import Hartree, kJ, mol, kcal
from autosteper.cage import name2seq, seq2name, Cage
from autosteper.optimizers import *
from autosteper.tools import strip_extraFullerene
from networkx.algorithms import isomorphism
from tqdm import tqdm
import math


def get_G(atoms: Atoms):
    nb_list = build_neighbor_list(atoms, self_interaction=False)
    mat = nb_list.get_connectivity_matrix()
    adj = mat.toarray()
    G = nx.from_numpy_array(adj)
    return G


def get_connection(workbase: str, dump:str, max_add: int, file_mid_formats: str, step: int, max_rank: int):
    os.chdir(workbase)
    dump = os.path.join(dump, 'connection')
    os.makedirs(dump, exist_ok=True)
    prev_atoms = read(f'{max_add}_{file_mid_formats}_0.xyz')
    prev_G_list = [get_G(prev_atoms)]
    for an_add in np.arange(max_add-step, -step, -step):
        new_G_list = []
        for a_rank in range(max_rank):
            cnt_list = []
            file_name = f'{an_add}_addons_dft_{a_rank}.xyz'
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
            np.save(file=os.path.join(dump, f'{an_add}_addons_dft_{a_rank}.npy'),
                    arr=np.array(cnt_list))
        prev_G_list = new_G_list


def simple_gauss_rename(dump_root: str, src_root: str, cage_size: int, group_size: int=None):
    """
    parse the gaussian output file, get energy info, and dump atoms in order.
    every file need to be started with 'addon num'+'_'
    read_gaussian_out is edited to get all the trajectory atoms
    dump_root: path to dump
    src_root: path of the original root
    cage_size: to get addon num
    group_size: in case the group have more than one atom
    """
    os.chdir(dump_root)
    os.makedirs('xyz')
    os.makedirs('log')
    os.makedirs('info')
    dump_xyz = os.path.abspath('xyz')
    dump_log = os.path.abspath('log')
    dump_info = os.path.abspath('info')
    e_map = {}
    e_name_map = {}
    add_num_list = []
    os.chdir(src_root)
    for a_file in os.listdir('./'):
        try:
            atoms = read_gaussian_out(fd=a_file)[-1]
            e = atoms.get_total_energy()
        except:
            # for the errored optimization file
            atoms = read_gaussian_out(fd=a_file)[-2]
            e = atoms.get_total_energy()
        if group_size:
            add_num = (len(atoms) - cage_size)/ group_size
        else:
            add_num = len(atoms) - cage_size

        if add_num not in e_map.keys():
            add_num_list.append(add_num)
            dummy_list = [None] * 10
            dummy_list[0] = e
            e_map.update({add_num: dummy_list})
        else:
            old_dummy_list = e_map[add_num]
            for idx, ii in enumerate(old_dummy_list):
                if ii == None:
                    break
            old_dummy_list[idx] = e
            e_map.update({add_num: old_dummy_list})
        e_name_map.update({e: a_file})
    new_e_map = dict()
    for a_key in sorted(add_num_list):
        e_list = e_map[a_key]
        new_e_list = []
        for a_e in e_list:
            if a_e:
                new_e_list.append(a_e)
        new_e_list.sort()
        # # For rename
        for idx, a_e in enumerate(new_e_list):
            atoms = read_gaussian_out(fd=e_name_map[a_e])[-1]
            write(filename=os.path.join(dump_xyz, f'{str(a_key)}_addons_dft_{str(idx)}.xyz'), images=atoms, format='xyz')
            shutil.copy(src=e_name_map[a_e],
                        dst=os.path.join(dump_log, f'{str(a_key)}_addons_dft_{str(idx)}.log'))
        for i in range(len(e_list) - len(new_e_list)):
            new_e_list.append(None)
        new_e_map.update({a_key: np.array(new_e_list)})
    info = pd.DataFrame(new_e_map)
    info.to_pickle(path=os.path.join(dump_info, 'info.pickle'))
    info.to_excel(os.path.join(dump_info, 'info.xlsx'))


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
                tgt_dummy_q_G = tgt_G.copy()
                tgt_dummy_q_G.remove_node(ii)
                tgt_dummy_q_G.remove_node(jj)
                if nx.is_isomorphic(tgt_dummy_q_G, dummy_q_G):
                    swr_container.append([ii, jj])
                    return swr_container


def find_SWR(q_root: str, tgt_root: str, dump_path: str, step: int, is_low_e:bool=None, is_unique: bool=None):
    '''
    Find SWR from atoms in q_root to atoms in tgt_root
    addon number is required in the file name, it better to be the renamed workbase.
    if is_unique is true, for every atoms in q_root, only one SWR target is outputed,
    typically for the lowest energy isomer, here we take the rank info in the name as criteria.
    if is_low_e is true, for every atoms in q_root, only one SWR target is outputed,
    and it should have a lower energy than the 'ought to be' parents.
    This feature need connection and energy info prepared.
    '''
    for a_file in os.listdir(q_root):
        if a_file.startswith('0_'):
            q_atoms = read(os.path.join(q_root, a_file))
            break
    for a_file in os.listdir(tgt_root):
        if a_file.startswith('0_'):
            tgt_atoms = read(os.path.join(tgt_root, a_file))
            break

    swr = map_SWR(q_atoms, tgt_atoms)
    q_G = get_G(q_atoms)
    tgt_G = get_G(tgt_atoms)

    q_swr = swr[0]
    q_swr_adj_nodes = []
    for a_node in q_swr:
        for an_adj in q_G[a_node]:
            if not an_adj in q_swr or an_adj in q_swr_adj_nodes:
                q_swr_adj_nodes.append(an_adj)
    q_swr_adj_nodes_set = set(q_swr_adj_nodes)

    tgt_swr = swr[1]

    dummy_q = q_G.copy()
    dummy_q.remove_nodes_from(q_swr)
    dummy_tgt = tgt_G.copy()
    dummy_tgt.remove_nodes_from(tgt_swr)
    GM = isomorphism.GraphMatcher(dummy_q, dummy_tgt)
    q_tgt_map = next(GM.match())

    swr_idx = 0
    for a_file in os.listdir(q_root):
        a_q_atoms = read(os.path.join(q_root, a_file))
        _, addon_set = strip_extraFullerene(a_q_atoms)
        ept_nodes = q_swr_adj_nodes_set - addon_set
        ept_nodes_for_chk = set()
        for a_node in ept_nodes:
            ept_nodes_for_chk.add(q_tgt_map[a_node])
        done_nodes = q_swr_adj_nodes_set - ept_nodes
        done_nodes_for_chk = set()
        for a_node in done_nodes:
            done_nodes_for_chk.add(q_tgt_map[a_node])

        if len(ept_nodes) > 0:
            a_q_G_for_chk = get_G(a_q_atoms)
            a_q_G_for_chk.remove_nodes_from(q_swr)

            if not (is_unique or is_low_e):
                for a_tgt_file in os.listdir(tgt_root):
                    tgt_add_num = int(a_tgt_file.split('_')[0])
                    if tgt_add_num == len(addon_set) + step:
                        a_tgt_atoms = read(os.path.join(tgt_root, a_tgt_file))
                        _, tgt_addon_set = strip_extraFullerene(a_tgt_atoms)
                        if len(done_nodes_for_chk - tgt_addon_set) == 0 and len(
                                ept_nodes_for_chk - tgt_addon_set) < len(
                                ept_nodes_for_chk):
                            a_tgt_G_for_chk = get_G(a_tgt_atoms)
                            a_GM = isomorphism.GraphMatcher(a_tgt_G_for_chk, a_q_G_for_chk)
                            if a_GM.subgraph_is_isomorphic():
                                os.chdir(dump_path)
                                os.makedirs(f'swr_{swr_idx}', exist_ok=True)
                                os.chdir(f'swr_{swr_idx}')
                                swr_idx += 1
                                write(filename='q_atoms.xyz', images=a_q_atoms, format='xyz')
                                write(filename='tgt_atoms.xyz', images=a_tgt_atoms, format='xyz')
            else:
                tgt_atoms_list = []
                rank_list = []
                for a_tgt_file in os.listdir(tgt_root):
                    tgt_add_num = int(a_tgt_file.split('_')[0])
                    if tgt_add_num == len(addon_set) + step:
                        a_tgt_atoms = read(os.path.join(tgt_root, a_tgt_file))
                        _, tgt_addon_set = strip_extraFullerene(a_tgt_atoms)
                        if len(done_nodes_for_chk - tgt_addon_set) == 0 and len(
                                ept_nodes_for_chk - tgt_addon_set) < len(
                                ept_nodes_for_chk):
                            a_tgt_G_for_chk = get_G(a_tgt_atoms)
                            a_GM = isomorphism.GraphMatcher(a_tgt_G_for_chk, a_q_G_for_chk)
                            if a_GM.subgraph_is_isomorphic():
                                tgt_atoms_list.append(a_tgt_atoms)
                                rank_list.append(int(a_tgt_file.split('.')[0].split('_')[-1]))
                if len(rank_list)>0:
                    if is_low_e:
                        tgt_info = pd.read_pickle(os.path.join(root, a_tgt, 'info', 'info.pickle'))
                        q_info = pd.read_pickle(os.path.join(root, a_q, 'info', 'info.pickle'))
                        tgt_e = tgt_info[len(addon_set) + step][sorted(rank_list)[0]]
                        cnt = np.load(os.path.join(root, a_q, 'connection', a_file[:-4] + '.npy'))
                        new_q_info = []
                        for idx, a_q_e in enumerate(q_info[len(addon_set) + step]):
                            if not a_q_e == None:
                                if cnt[idx] == 1:
                                    new_q_info.append(a_q_e)
                        if len(new_q_info) != 0:
                            if tgt_e > min(new_q_info):
                                continue
                    os.chdir(dump_path)
                    os.makedirs(f'swr_{swr_idx}', exist_ok=True)
                    os.chdir(f'swr_{swr_idx}')
                    swr_idx += 1
                    rank_atoms_map = dict(zip(rank_list, tgt_atoms_list))
                    write(filename='q_atoms.xyz', images=a_q_atoms, format='xyz')
                    write(filename='tgt_atoms.xyz', images=rank_atoms_map[sorted(rank_list)[0]], format='xyz')


def test_Boltzmann(workbase: str, dump_path: str, step: int, first_addon: int,
                   q_atoms: Atoms,
                   cut_rank: int, cut_value: float, cht_path):
    print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
          'The atom sequence in two cages should be identical.\n'
          'The addons should be atomic addons.\n'
          'Complex senerio could be tackled with strip_extraFullerene function in tools module.')

    pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    cage_size = len(pristine_cage)
    max_add_num = len(q_add_set)

    add_list = list(range(first_addon, max_add_num+step, step))

    sub_pcts = []
    for an_addon in add_list:
        addon_path = str(an_addon) + '_addons'
        os.chdir(cht_path)
        rel_e_all_cht = np.load(file=addon_path+'_all_rel_e.npy')
        rel_sub = np.load(file=addon_path+'_sub_rel_e.npy')
        rel_e_all_cht = np.sort(rel_e_all_cht)
        rel_sub = sorted(rel_sub)
        rel_sub_set = set(rel_sub)
        idx_list_cht = []
        for idx, a_e in enumerate(rel_e_all_cht):
            if a_e in rel_sub_set:
                idx_list_cht.append(idx)

        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]
        rel_e_all = np.array(deep_info['energy'][:cut_rank]) - min_e
        rel_e_all = rel_e_all*Hartree
        sub_idx_list = []
        print(q_add_set)


        for idx, a_name in enumerate(deep_info['name'][:cut_rank]):
            rel_e = rel_e_all[idx]

            if rel_e > cut_value:
                break
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
            print(pred_add_set)
            if len(pred_add_set - q_add_set) == 0:
                sub_idx_list.append(idx)


        print(addon_path)
        print(idx_list_cht)
        if len(idx_list_cht)==0:
            sub_pcts.append(0)
        else:
            # print(rel_e_all)
            # a_sub_pct = get_Boltzmann(idx_list=sub_idx_list, e_arr=rel_e_all)
            a_sub_pct = get_Boltzmann(idx_list=idx_list_cht, e_arr=rel_e_all_cht)
            sub_pcts.append(a_sub_pct)
    os.chdir(dump_path)
    # sub_pct_info = pd.DataFrame(dict(zip(add_list, sub_pcts)))
    # sub_pct_info.to_excel('pcts.xlsx')
    sub_pcts_arr = np.array(sub_pcts)
    np.save(file='sub_pcts.npy', arr=sub_pcts_arr)
    fig = plt.figure(dpi=400, figsize=(len(add_list), 7))
    plt.plot(add_list, sub_pcts, '+-')
    print(sub_pcts_arr)
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure Boltzmann Distribution')
    # plt.legend(loc='best')
    # plt.title(f'{dump_pic_path[:-4]}')
    plt.savefig('sub_pcts.png')
    return sub_pcts_arr



def test_iso(q_atoms: Atoms, deep_info_path: str, top_k: int=None):
    print('Please make sure the queried atoms has been optimized with at least semi-empirical level of theory.')
    target_G = get_G(q_atoms)
    deep_info = pd.read_pickle(deep_info_path)
    flag = 0
    for idx, a_xyz_path in enumerate(tqdm(deep_info['xyz_path'][:top_k])):
        an_atoms = read(a_xyz_path)
        pred_G = get_G(an_atoms)
        if nx.is_isomorphic(pred_G, target_G):
            rel_e = (deep_info['energy'][idx] - deep_info['energy'][0])*Hartree
            print('The queried atoms are found in AutoStepper scan results.\n'
                  f'The target is in rank {str(idx+1)}.\n'
                  f'Has a relative energy of {str(rel_e)}eV considering the lowest isomer in scan results.\n')
            flag = 1
            break
    if flag == 0:
        if top_k:
            print(f'The queried atoms are NOT found in top {top_k} AutoStepper scan results.\n'
                  'Please check its substructure with test_sub function or enlarge test number.')
        else:
            print('The queried atoms are NOT found in AutoStepper scan results.\n'
                  'Please check its substructure with test_sub function.')
        return None
    else:
        return idx+1


def test_sub(workbase: str, dump_file_path: str, step: int, first_addon: int, cage_size: int=None,q_atoms: Atoms=None, q_seq: str=None, q_cage: Cage=None, q_add_set: set=None):
    if q_atoms:
        print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
              'The atom sequence in two cages should be identical.\n'
              'The addons should be atomic addons.\n'
              'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
        pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    elif q_seq:
        _, q_add_set, _0 = seq2name(q_seq, q_cage)
    elif not q_add_set:
        print('Please input query addon set.\n'
              'Currently only surpport ase.atoms input, sequence input and addon set input.')
    max_add_num = len(q_add_set)
    if cage_size == None:
        if q_cage:
            cage_size = q_cage.size
        elif q_atoms:
            cage_size = len(pristine_cage)
        else:
            print('Please input cage size.')

    for an_addon in range(first_addon, max_add_num+step, step):
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)
        with open(dump_file_path, 'a') as f:
            deep_info = pd.read_pickle('deep_yes_info.pickle')
            min_e = deep_info['energy'][0]
            f.write(f'Start screen substructures in {addon_path}:\n')
            flag = 0
            for idx, a_name in enumerate(deep_info['name']):
                _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
                if len(pred_add_set - q_add_set) == 0:
                    f.write('A substructrue is found:\n')
                    f.write(f'The absent addon sites are:      {str(q_add_set - pred_add_set)}.\n')
                    f.write(f'In rank:                         {str(idx+1)}\n')
                    rel_e = (deep_info['energy'][idx] - min_e)*Hartree
                    f.write(f'Relative energy is:              {str(rel_e)}eV\n')
                    f.write(f'Name is:                         {a_name}\n')
                    f.write('\n')
                    f.write('\n')
                    flag = 1
                    # break
            if flag == 0:
                f.write(f'{addon_path} do not have a substructure.\n')

def simple_test_sub_plot(workbase: str, dump_pic_path: str, step: int, first_addon: int, cage_size: int=None,q_atoms: Atoms=None, q_seq: str=None, q_cage: Cage=None, q_add_set: set=None, cut_rank: int=None, cut_value: float=None):
    if q_atoms:
        print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
              'The atom sequence in two cages should be identical.\n'
              'The addons should be atomic addons.\n'
              'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
        _, q_add_set = strip_extraFullerene(q_atoms)
    elif q_seq:
        _, q_add_set, _0 = seq2name(q_seq, q_cage)
    elif not q_add_set:
        print('Please input query addon set.\n'
              'Currently only surpport ase.atoms input, sequence input and addon set input.')
    max_add_num = len(q_add_set)
    if cage_size == None:
        if q_cage:
            cage_size = q_cage.size
        else:
            print('Please input cage size.')

    add_list = list(range(first_addon, max_add_num+step, step))
    fig = plt.figure(dpi=400, figsize=(len(add_list), 7))
    is_first = True
    seed_info_dict = dict()
    non_seed_info_dict = dict()
    for an_addon in add_list:
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]

        rel_seed = []
        rel_non_seed = []
        old_rel_seed = []
        old_rel_non_seed = []

        for idx, a_name in enumerate(deep_info['name']):
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
            if len(pred_add_set - q_add_set) == 0:
                rel_e = (deep_info['energy'][idx] - min_e) * Hartree
                old_rel_e = (deep_info['energy'][idx]) * Hartree
                if idx < cut_rank and rel_e < cut_value:
                    old_rel_seed.append(old_rel_e)
                    rel_seed.append(rel_e)
                else:
                    old_rel_non_seed.append(old_rel_e)
                    rel_non_seed.append(rel_e)

        if len(rel_seed)+len(rel_non_seed)>0:
            if is_first:
                plt.scatter(x=[an_addon] * len(rel_seed), y=rel_seed, marker='+', c='red', label='seed')
                plt.scatter(x=[an_addon] * len(rel_non_seed), y=rel_non_seed, marker='x', c='blue', label='non-seed')
                is_first = False
            else:
                plt.scatter(x=[an_addon] * len(rel_seed), y=rel_seed, marker='+', c='red')
                plt.scatter(x=[an_addon] * len(rel_non_seed), y=rel_non_seed, marker='x', c='blue')
        seed_info_dict.update({an_addon: [old_rel_seed]})
        non_seed_info_dict.update({an_addon: [old_rel_non_seed]})
    seed_info_df = pd.DataFrame(seed_info_dict)
    non_seed_info_df = pd.DataFrame(non_seed_info_dict)
    seed_info_df.to_pickle(dump_pic_path[:-4]+'_seed.pickle')
    non_seed_info_df.to_pickle(dump_pic_path[:-4]+'_nonseed.pickle')

    plt.legend(loc='best')
    plt.ylim(None,1)
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    # plt.title(f'{dump_pic_path[:-4]}')
    plt.savefig(dump_pic_path)


def post_addon_refine(workbase: str, step: int, first_addon: int,
               ref_workbase: str, refine_optimizer: Optimizer, max_add_num: int,
               q_atoms: Atoms = None, cut_rank: int = None, cut_value: float = None):
    print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
          'The atom sequence in two cages should be identical.\n'
          'The addons should be atomic addons.\n'
          'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
    pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    q_G = get_G(q_atoms)

    cage_size = len(pristine_cage)

    add_list = list(range(first_addon, max_add_num+step, step))

    for an_addon in add_list:
        if an_addon % 2 == 1:
            dump_path = os.path.join(ref_workbase, 'odd_raw')
        else:
            dump_path = os.path.join(ref_workbase, 'even_raw')

        os.makedirs(dump_path, exist_ok=True)
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]
        os.chdir(dump_path)
        for idx, a_name in enumerate(deep_info['name']):
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
            rel_e = (deep_info['energy'][idx] - min_e) * Hartree
            old_rel_e = (deep_info['energy'][idx]) * Hartree
            if idx < cut_rank and rel_e < cut_value:
                a_xyz_path = deep_info['xyz_path'][idx]
                if len(q_add_set - pred_add_set) == 0:
                    shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
                else:
                    pred_atoms = read(a_xyz_path)
                    pred_G = get_G(pred_atoms)
                    GM = isomorphism.GraphMatcher(pred_G, q_G)
                    if GM.subgraph_is_isomorphic():
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
                    else:
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_non_sub_rank_{idx + 1}.xyz')
            elif len(q_add_set - pred_add_set) == 0:
                a_xyz_path = deep_info['xyz_path'][idx]
                shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
    os.chdir(ref_workbase)
    os.makedirs('opt_workbase', exist_ok=True)
    os.makedirs('sub_structures', exist_ok=True)
    opt_workbase = os.path.abspath('opt_workbase')
    if os.path.exists('even_raw'):
        refine_optimizer.run_a_batch(path_source='even_raw', path_destination='opt_workbase',
                                     cmd_list=refine_optimizer.cmd_list)
    if os.path.exists('odd_raw'):
        refine_optimizer.run_a_batch(path_source='odd_raw', path_destination='opt_workbase',
                                     cmd_list=[*refine_optimizer.cmd_list, '--uhf 1'])
    fig = plt.figure(dpi=400, figsize=(len(add_list)+1, 7))
    is_first = True
    for an_addon in add_list:
        e_list = []
        name_list = []
        sub_list_e = []
        non_sub_list_e = []
        os.chdir(opt_workbase)
        for a_folder in os.listdir('./'):
            if a_folder.startswith(f'{an_addon}_addons'):
                with open(os.path.join(a_folder, 'xtbopt.xyz'), 'r') as f:
                    f.readline()
                    energy_line = f.readline()
                    energy = float(energy_line.split()[1])*Hartree
                    e_list.append(energy)
                name_list.append(a_folder)
        if len(e_list) == 0:
            continue
        os.chdir(ref_workbase)
        os.chdir('sub_structures')
        os.makedirs(f'{an_addon}', exist_ok=True)
        e_arr = np.array(e_list) - min(e_list)
        e_name_map = dict(zip(e_arr, name_list))
        e_arr = sorted(e_arr)
        for idx, a_e in enumerate(e_arr):
            a_folder = e_name_map[a_e]
            if 'non' in a_folder.split('_'):
                non_sub_list_e.append(a_e)
                shutil.copy(src=os.path.join(opt_workbase, a_folder, 'xtbopt.xyz'),
                            dst=os.path.join(f'{an_addon}', f'non_sub_rank_{idx}.xyz'))
            else:
                sub_list_e.append(a_e)
                shutil.copy(src=os.path.join(opt_workbase, a_folder, 'xtbopt.xyz'),
                            dst=os.path.join(f'{an_addon}', f'is_sub_rank_{idx}.xyz'))
        np.save(file=f'{an_addon}_addons_all_rel_e.npy', arr=e_arr)
        np.save(file=f'{an_addon}_addons_sub_rel_e.npy', arr=np.array(sub_list_e))
        np.save(file=f'{an_addon}_addons_non_sub_rel_e.npy', arr=np.array(non_sub_list_e))
        if is_first:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue', label='non-sub')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red', label='sub')
            is_first = False
        else:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red')

        os.chdir(ref_workbase)

    plt.ylim(None,1)
    plt.legend(loc='best')
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    plt.savefig('sub_rel_e.png')


def post_addon(workbase: str, step: int, first_addon: int, max_add_num: int, old_workbase: str,
               q_atoms: Atoms = None, cut_rank: int = None, cut_value: float = None):
    print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
          'The atom sequence in two cages should be identical.\n'
          'The addons should be atomic addons.\n'
          'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
    pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    q_G = get_G(q_atoms)

    cage_size = len(pristine_cage)

    add_list = list(range(first_addon, max_add_num+step, step))

    os.chdir(workbase)
    os.makedirs('sub_structures', exist_ok=True)
    dump_path = os.path.abspath('sub_structures')

    fig = plt.figure(dpi=400, figsize=(len(add_list)+1, 7))
    is_first = True

    for an_addon in add_list:
        sub_list_e = []
        non_sub_list_e = []

        addon_path = str(an_addon) + 'addons'
        os.chdir(old_workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]
        os.chdir(dump_path)
        for idx, a_name in enumerate(deep_info['name']):
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
            rel_e = (deep_info['energy'][idx] - min_e) * Hartree
            old_rel_e = (deep_info['energy'][idx]) * Hartree
            if idx < cut_rank and rel_e < cut_value:
                a_xyz_path = deep_info['xyz_path'][idx]
                if len(q_add_set - pred_add_set) == 0:
                    shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
                    sub_list_e.append(rel_e)
                else:
                    pred_atoms = read(a_xyz_path)
                    pred_G = get_G(pred_atoms)
                    GM = isomorphism.GraphMatcher(pred_G, q_G)
                    if GM.subgraph_is_isomorphic():
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
                        sub_list_e.append(rel_e)
                    else:
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_non_sub_rank_{idx + 1}.xyz')
                        non_sub_list_e.append(rel_e)
            else:
                break
            # Caution!!!
            # elif len(q_add_set - pred_add_set) == 0:
            #     a_xyz_path = deep_info['xyz_path'][idx]
            #     shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
            #     sub_list_e.append(rel_e)
        os.chdir(workbase)
        np.save(file=f'{an_addon}_addons_sub_rel_e.npy', arr=np.array(sub_list_e))
        np.save(file=f'{an_addon}_addons_non_sub_rel_e.npy', arr=np.array(non_sub_list_e))
        if is_first:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue', label='non-sub')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red', label='sub')
            is_first = False
        else:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red')
    plt.ylim(None,1)
    plt.legend(loc='best')
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    plt.savefig('sub_rel_e.png')


def get_connection(workbase: str, step: int, first_addon: int, max_add_num: int, old_workbase: str,
               q_atoms: Atoms = None, cut_rank: int = None, cut_value: float = None):
    print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
          'The atom sequence in two cages should be identical.\n'
          'The addons should be atomic addons.\n'
          'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
    pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    q_G = get_G(q_atoms)

    cage_size = len(pristine_cage)

    add_list = list(range(first_addon, max_add_num+step, step))

    os.chdir(workbase)
    os.makedirs('sub_structures', exist_ok=True)
    dump_path_sub = os.path.abspath('sub_structures')
    os.makedirs('non_sub_structures', exist_ok=True)
    dump_path_non = os.path.abspath('non_sub_structures')
    os.makedirs('info', exist_ok=True)
    dump_path_info = os.path.abspath('info')
    fig = plt.figure(dpi=400, figsize=(len(add_list)+1, 7))
    is_first = True

    for an_addon in add_list:
        sub_list_e = []
        non_sub_list_e = []

        addon_path = str(an_addon) + 'addons'
        os.chdir(old_workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]

        for idx, a_name in enumerate(deep_info['name']):
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
            rel_e = (deep_info['energy'][idx] - min_e) * Hartree
            old_rel_e = (deep_info['energy'][idx]) * Hartree
            if idx < cut_rank and rel_e < cut_value:
                a_xyz_path = deep_info['xyz_path'][idx]
                if an_addon >= len(q_add_set):
                    if len(q_add_set - pred_add_set) == 0:
                        os.chdir(dump_path_sub)
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_rank_{idx + 1}.xyz')
                        sub_list_e.append(rel_e)
                    else:
                        pred_atoms = read(a_xyz_path)
                        pred_G = get_G(pred_atoms)
                        GM = isomorphism.GraphMatcher(pred_G, q_G)
                        if GM.subgraph_is_isomorphic():
                            os.chdir(dump_path_sub)
                            shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_rank_{idx + 1}.xyz')
                            sub_list_e.append(rel_e)
                        else:
                            os.chdir(dump_path_non)
                            shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_rank_{idx + 1}.xyz')
                            non_sub_list_e.append(rel_e)
                else:
                    if len(pred_add_set - q_add_set) == 0:
                        os.chdir(dump_path_sub)
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_rank_{idx + 1}.xyz')
                        sub_list_e.append(rel_e)
                    else:
                        pred_atoms = read(a_xyz_path)
                        pred_G = get_G(pred_atoms)
                        GM = isomorphism.GraphMatcher(q_G, pred_G)
                        if GM.subgraph_is_isomorphic():
                            os.chdir(dump_path_sub)
                            shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_rank_{idx + 1}.xyz')
                            sub_list_e.append(rel_e)
                        else:
                            os.chdir(dump_path_non)
                            shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_rank_{idx + 1}.xyz')
                            non_sub_list_e.append(rel_e)
            else:
                break

        os.chdir(dump_path_info)
        np.save(file=f'{an_addon}_addons_sub_rel_e.npy', arr=np.array(sub_list_e))
        np.save(file=f'{an_addon}_addons_non_sub_rel_e.npy', arr=np.array(non_sub_list_e))
        if is_first:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue', label='non-sub')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red', label='sub')
            is_first = False
        else:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red')
    plt.ylim(None,1)
    plt.legend(loc='best')
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    os.chdir(workbase)
    plt.savefig('sub_rel_e.png')


def get_low_e_xyz(dump_path: str, add_num: int, old_workbase: str, cut_rank: int, cut_value: float):
    addon_path = str(add_num) + 'addons'
    os.chdir(old_workbase)
    os.chdir(addon_path)
    deep_info = pd.read_pickle('deep_yes_info.pickle')
    min_e = deep_info['energy'][0]

    os.chdir(dump_path)
    for idx, a_name in enumerate(deep_info['name']):
        rel_e = (deep_info['energy'][idx] - min_e) * Hartree
        if idx < cut_rank and rel_e < cut_value:
            a_xyz_path = deep_info['xyz_path'][idx]
            shutil.copy(src=a_xyz_path, dst=f'{add_num}_addons_rank_{idx + 1}.xyz')


def test_sub_with_iso(workbase: str, dump_pic_path: str, step: int, first_addon: int, q_atoms: Atoms=None, cut_rank: int=None, cut_value: float=None, is_log: bool=None, dump_path: str=None):
    print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
          'The atom sequence in two cages should be identical.\n'
          'The addons should be atomic addons.\n'
          'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
    pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    q_G = get_G(q_atoms)
    max_add_num = len(q_add_set)
    cage_size = len(pristine_cage)

    add_list = list(range(first_addon, max_add_num+step, step))
    fig = plt.figure(dpi=400, figsize=(len(add_list), 7))
    is_first = True
    seed_info_dict = dict()
    non_seed_info_dict = dict()
    for an_addon in add_list:
        dump_addon_path = os.path.join(dump_path, str(an_addon))
        os.makedirs(dump_addon_path, exist_ok=True)
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]

        rel_seed = []
        rel_non_seed = []
        old_rel_seed = []
        old_rel_non_seed = []

        if is_log:
            for idx, a_name in tqdm(enumerate(deep_info['name'])):
                _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
                rel_e = (deep_info['energy'][idx] - min_e) * Hartree
                old_rel_e = (deep_info['energy'][idx]) * Hartree
                if idx < cut_rank and rel_e < cut_value:
                    a_xyz_path = deep_info['xyz_path'][idx]
                    if len(pred_add_set - q_add_set) == 0:
                        old_rel_seed.append(old_rel_e)
                        rel_seed.append(rel_e)
                        shutil.copy(src=a_xyz_path, dst=os.path.join(dump_addon_path, f'seed_rank_{idx+1}.xyz'))
                    else:
                        pred_atoms = read(a_xyz_path)
                        pred_G = get_G(pred_atoms)
                        GM = isomorphism.GraphMatcher(q_G, pred_G)
                        if GM.subgraph_is_isomorphic():
                            old_rel_seed.append(old_rel_e)
                            rel_seed.append(rel_e)
                            shutil.copy(src=a_xyz_path, dst=os.path.join(dump_addon_path, f'seed_rank_{idx + 1}.xyz'))
                elif len(pred_add_set - q_add_set) == 0:
                    a_xyz_path = deep_info['xyz_path'][idx]
                    old_rel_non_seed.append(old_rel_e)
                    rel_non_seed.append(rel_e)
                    shutil.copy(src=a_xyz_path, dst=os.path.join(dump_addon_path, f'non_seed_rank_{idx + 1}.xyz'))
        else:
            for idx, a_name in enumerate(deep_info['name']):
                _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
                rel_e = (deep_info['energy'][idx] - min_e) * Hartree
                old_rel_e = (deep_info['energy'][idx]) * Hartree
                if idx < cut_rank and rel_e < cut_value:
                    if len(pred_add_set - q_add_set) == 0:
                        old_rel_seed.append(old_rel_e)
                        rel_seed.append(rel_e)
                    else:
                        a_xyz_path = deep_info['xyz_path'][idx]
                        pred_atoms = read(a_xyz_path)
                        pred_G = get_G(pred_atoms)
                        GM = isomorphism.GraphMatcher(q_G, pred_G)
                        if GM.subgraph_is_isomorphic():
                            old_rel_seed.append(old_rel_e)
                            rel_seed.append(rel_e)
                elif len(pred_add_set - q_add_set) == 0:
                    old_rel_non_seed.append(old_rel_e)
                    rel_non_seed.append(rel_e)

        if len(rel_seed)+len(rel_non_seed)>0:
            if is_first:
                plt.scatter(x=[an_addon] * len(rel_seed), y=rel_seed, marker='+', c='red', label='seed')
                plt.scatter(x=[an_addon] * len(rel_non_seed), y=rel_non_seed, marker='x', c='blue', label='non-seed')
                is_first = False
            else:
                plt.scatter(x=[an_addon] * len(rel_seed), y=rel_seed, marker='+', c='red')
                plt.scatter(x=[an_addon] * len(rel_non_seed), y=rel_non_seed, marker='x', c='blue')
        seed_info_dict.update({an_addon: [old_rel_seed]})
        non_seed_info_dict.update({an_addon: [old_rel_non_seed]})
    seed_info_df = pd.DataFrame(seed_info_dict)
    non_seed_info_df = pd.DataFrame(non_seed_info_dict)
    seed_info_df.to_pickle(dump_pic_path[:-4]+'_seed.pickle')
    non_seed_info_df.to_pickle(dump_pic_path[:-4]+'_nonseed.pickle')

    plt.legend(loc='best')
    plt.ylim(None,1)
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    # plt.title(f'{dump_pic_path[:-4]}')
    plt.savefig(dump_pic_path)


def refine_sub(workbase: str, step: int, first_addon: int,
               ref_workbase: str, refine_optimizer: Optimizer,
               q_atoms: Atoms = None, cut_rank: int = None, cut_value: float = None):
    print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
          'The atom sequence in two cages should be identical.\n'
          'The addons should be atomic addons.\n'
          'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
    pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    q_G = get_G(q_atoms)
    max_add_num = len(q_add_set)
    cage_size = len(pristine_cage)

    add_list = list(range(first_addon, max_add_num+step, step))

    for an_addon in add_list:
        if an_addon % 2 == 1:
            dump_path = os.path.join(ref_workbase, 'odd_raw')
        else:
            dump_path = os.path.join(ref_workbase, 'even_raw')

        os.makedirs(dump_path, exist_ok=True)
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]
        os.chdir(dump_path)
        for idx, a_name in enumerate(deep_info['name']):
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
            rel_e = (deep_info['energy'][idx] - min_e) * Hartree
            old_rel_e = (deep_info['energy'][idx]) * Hartree
            if idx < cut_rank and rel_e < cut_value:
                a_xyz_path = deep_info['xyz_path'][idx]
                if len(pred_add_set - q_add_set) == 0:
                    shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
                else:
                    pred_atoms = read(a_xyz_path)
                    pred_G = get_G(pred_atoms)
                    GM = isomorphism.GraphMatcher(q_G, pred_G)
                    if GM.subgraph_is_isomorphic():
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
                    else:
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_non_sub_rank_{idx + 1}.xyz')
            elif len(pred_add_set - q_add_set) == 0:
                a_xyz_path = deep_info['xyz_path'][idx]
                shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{idx + 1}.xyz')
    os.chdir(ref_workbase)
    os.makedirs('opt_workbase', exist_ok=True)
    os.makedirs('sub_structures', exist_ok=True)
    opt_workbase = os.path.abspath('opt_workbase')
    if os.path.exists('even_raw'):
        refine_optimizer.run_a_batch(path_source='even_raw', path_destination='opt_workbase',
                                     cmd_list=refine_optimizer.cmd_list)
    if os.path.exists('odd_raw'):
        refine_optimizer.run_a_batch(path_source='odd_raw', path_destination='opt_workbase',
                                     cmd_list=[*refine_optimizer.cmd_list, '--uhf 1'])
    fig = plt.figure(dpi=400, figsize=(len(add_list), 7))
    is_first = True
    for an_addon in add_list:
        e_list = []
        name_list = []
        sub_list_e = []
        non_sub_list_e = []
        os.chdir(opt_workbase)
        for a_folder in os.listdir('./'):
            if a_folder.startswith(f'{an_addon}_addons'):
                with open(os.path.join(a_folder, 'xtbopt.xyz'), 'r') as f:
                    f.readline()
                    energy_line = f.readline()
                    energy = float(energy_line.split()[1])*Hartree
                    e_list.append(energy)
                name_list.append(a_folder)
        if len(e_list) == 0:
            continue
        os.chdir(ref_workbase)
        os.chdir('sub_structures')
        os.makedirs(f'{an_addon}', exist_ok=True)
        e_arr = np.array(e_list) - min(e_list)
        e_name_map = dict(zip(e_arr, name_list))
        e_arr = sorted(e_arr)
        for idx, a_e in enumerate(e_arr):
            a_folder = e_name_map[a_e]
            if 'non' in a_folder.split('_'):
                non_sub_list_e.append(a_e)
                shutil.copy(src=os.path.join(opt_workbase, a_folder, 'xtbopt.xyz'),
                            dst=os.path.join(f'{an_addon}', f'non_sub_rank_{idx}.xyz'))
            else:
                sub_list_e.append(a_e)
                shutil.copy(src=os.path.join(opt_workbase, a_folder, 'xtbopt.xyz'),
                            dst=os.path.join(f'{an_addon}', f'is_sub_rank_{idx}.xyz'))
        np.save(file=f'{an_addon}_addons_all_rel_e.npy', arr=e_arr)
        np.save(file=f'{an_addon}_addons_sub_rel_e.npy', arr=np.array(sub_list_e))
        np.save(file=f'{an_addon}_addons_non_sub_rel_e.npy', arr=np.array(non_sub_list_e))
        if is_first:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue', label='non-sub')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red', label='sub')
            is_first = False
        else:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red')

        os.chdir(ref_workbase)

    plt.ylim(None,1)
    plt.legend(loc='best')
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    plt.savefig('sub_rel_e.png')

def refine_sub_v2(workbase: str, step: int, first_addon: int,max_add_num:int,dump_all_xyz:str,new_name:str,
               ref_workbase: str, refine_optimizer: Optimizer,
               q_atoms: Atoms = None, cut_rank: int = None, cut_value: float = None):
    print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
          'The atom sequence in two cages should be identical.\n'
          'The addons should be atomic addons.\n'
          'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
    pristine_cage, q_add_set = strip_extraFullerene(q_atoms)
    q_G = get_G(q_atoms)
    # max_add_num = len(q_add_set)
    cage_size = len(pristine_cage)

    add_list = list(range(first_addon, max_add_num+step, step))

    for an_addon in add_list:
        if an_addon % 2 == 1:
            dump_path = os.path.join(ref_workbase, 'odd_raw')
        else:
            dump_path = os.path.join(ref_workbase, 'even_raw')

        os.makedirs(dump_path, exist_ok=True)
        addon_path = str(an_addon) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)

        deep_info = pd.read_pickle('deep_yes_info.pickle')
        min_e = deep_info['energy'][0]
        os.chdir(dump_path)
        jjj = 0
        for idx, a_name in enumerate(deep_info['name']):
            if jjj == 10:
                break
            _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
            rel_e = (deep_info['energy'][idx] - min_e) * Hartree
            old_rel_e = (deep_info['energy'][idx]) * Hartree
            if idx < cut_rank and rel_e < cut_value:
                a_xyz_path = deep_info['xyz_path'][idx]
                if len(pred_add_set - q_add_set) == 0:
                    shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{jjj}.xyz')
                    jjj = jjj + 1
                else:
                    pred_atoms = read(a_xyz_path)
                    pred_G = get_G(pred_atoms)
                    GM = isomorphism.GraphMatcher(q_G, pred_G)
                    if GM.subgraph_is_isomorphic():
                        shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{jjj}.xyz')
                        jjj = jjj + 1
                    # else:
                    #     shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_non_sub_rank_{jjj}.xyz')
            elif len(pred_add_set - q_add_set) == 0:
                a_xyz_path = deep_info['xyz_path'][idx]
                shutil.copy(src=a_xyz_path, dst=f'{an_addon}_addons_sub_rank_{jjj}.xyz')
                jjj = jjj + 1
    os.chdir(ref_workbase)
    os.makedirs('opt_workbase', exist_ok=True)
    os.makedirs('sub_structures', exist_ok=True)
    opt_workbase = os.path.abspath('opt_workbase')
    if os.path.exists('even_raw'):
        refine_optimizer.run_a_batch(path_source='even_raw', path_destination='opt_workbase',
                                     cmd_list=refine_optimizer.cmd_list)
    if os.path.exists('odd_raw'):
        refine_optimizer.run_a_batch(path_source='odd_raw', path_destination='opt_workbase',
                                     cmd_list=[*refine_optimizer.cmd_list, '--uhf 1'])
    fig = plt.figure(dpi=400, figsize=(len(add_list), 7))
    is_first = True
    for an_addon in add_list:
        e_list = []
        name_list = []
        sub_list_e = []
        non_sub_list_e = []
        os.chdir(opt_workbase)
        for a_folder in os.listdir('./'):
            if a_folder.startswith(f'{an_addon}_addons'):
                with open(os.path.join(a_folder, 'xtbopt.xyz'), 'r') as f:
                    f.readline()
                    energy_line = f.readline()
                    energy = float(energy_line.split()[1])*Hartree
                    e_list.append(energy)
                name_list.append(a_folder)
        if len(e_list) == 0:
            continue
        os.chdir(ref_workbase)
        os.chdir('sub_structures')
        os.makedirs(f'{an_addon}', exist_ok=True)
        e_arr = np.array(e_list) - min(e_list)
        e_name_map = dict(zip(e_arr, name_list))
        e_arr = sorted(e_arr)
        for idx, a_e in enumerate(e_arr):
            a_folder = e_name_map[a_e]
            if 'non' in a_folder.split('_'):
                non_sub_list_e.append(a_e)
                shutil.copy(src=os.path.join(opt_workbase, a_folder, 'xtbopt.xyz'),
                            dst=os.path.join(f'{an_addon}', f'non_sub_rank_{idx}.xyz'))
            else:
                sub_list_e.append(a_e)
                shutil.copy(src=os.path.join(opt_workbase, a_folder, 'xtbopt.xyz'),
                            dst=os.path.join(f'{an_addon}', f'is_sub_rank_{idx}.xyz'))
                aa_list = a_folder.split('_')
                shutil.copy(src=os.path.join(opt_workbase, a_folder, 'xtbopt.xyz'),
                            dst=os.path.join(f'{dump_all_xyz}', f'{aa_list[0]}_addons_{new_name}_rank_{aa_list[-1]}.xyz'))

        np.save(file=f'{an_addon}_addons_all_rel_e.npy', arr=e_arr)
        np.save(file=f'{an_addon}_addons_sub_rel_e.npy', arr=np.array(sub_list_e))
        np.save(file=f'{an_addon}_addons_non_sub_rel_e.npy', arr=np.array(non_sub_list_e))
        if is_first:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue', label='non-sub')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red', label='sub')
            is_first = False
        else:
            plt.scatter(x=[an_addon] * len(non_sub_list_e), y=non_sub_list_e, marker='+', c='blue')
            plt.scatter(x=[an_addon] * len(sub_list_e), y=sub_list_e, marker='x', c='red')

        os.chdir(ref_workbase)

    plt.ylim(None,1)
    plt.legend(loc='best')
    plt.xticks(ticks=add_list, labels=add_list)
    plt.xlabel('Addon number')
    plt.ylabel('Substructure relative energy (ev).')
    plt.savefig('sub_rel_e.png')


hartree2kjmol = Hartree / (kJ / mol)
class Path_Parser():
    def __init__(self, path_para: dict, step: int, start: int, workbase: str, q_cage: Cage, optimizer: Optimizer, refine_para: dict=None):
        self.q_cage = q_cage
        self.optimizer = optimizer
        self.workbase = os.path.join(workbase, f'path_info_{self.q_cage.name}')
        os.makedirs(self.workbase, exist_ok=True)
        self.step = step
        self.start = start
        self.q_path_num = path_para['q_path_num']
        self.q_add_num = path_para['q_add_num']
        self.q_low_e_num = path_para['q_low_e_num']
        self.log_low_e_num = path_para['log_low_e_num']
        if 'ctl_path_para' in path_para.keys():
            self.ctr_path = True
            self.ctl_parent_num = path_para['ctl_path_para']['ctl_parent_num']
            self.max_path_num = path_para['ctl_path_para']['max_path_num']
        else:
            self.ctr_path = False

        if 'refine_para' in path_para.keys():
            self.is_refine = True
            self.refine_top_k = path_para['refine_para']['refine_top_k']
            self.refine_mode = path_para['refine_para']['opt_mode']
            if self.refine_mode == 'xtb':
                self.refine_optimizer = XTB_Optimizer(path_para['refine_para']['opt_para'], checker=None, cage=None)
            elif self.refine_mode == 'ase':
                self.refine_optimizer = ASE_Optimizer(path_para['refine_para']['opt_para'], checker=None, cage=None)
            elif self.refine_mode == 'gaussian':
                self.refine_optimizer = Gaussian_Optimizer(path_para['refine_para']['opt_para'], checker=None, cage=None)
        else:
            self.is_refine = False


    def _sort_parent(self, flat_info: pd.DataFrame, info_path: str, parent_flat_info: pd.DataFrame):
        for a_cage in flat_info.keys():
            parent_names = flat_info[a_cage][0]
            name_e_map = {}
            for a_parent in parent_names:
                name_e_map.update({parent_flat_info[a_parent][1]: a_parent})
            new_parent_names = []
            for a_e in sorted(name_e_map.keys())[:self.ctl_parent_num]:
                new_parent_names.append(name_e_map[a_e])
            flat_info[a_cage][0] = new_parent_names
        flat_info.to_pickle(info_path)
        return flat_info



    def get_path_info(self):
        os.chdir(self.workbase)
        def _get_e_area_triangle(e_list: list):
            e_array = np.array(e_list)
            e_array = e_array - self.base_e
            e_area = e_array[1] / 2
            for i in range(2, len(e_array)):
                e_area += (e_array[i] + e_array[i - 1]) / 2
            return e_area * hartree2kjmol

        def _get_e_area(e_list: list):
            e_area = sum(e_list)
            return e_area

        def _get_pathway_unit(idx: int, pathway: list, cage_name: str, e_list: list, name_list: list):
            if idx == 0:
                _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                final_pathway = [init_pathway, addon_set, *pathway]
                q_isomer_e = flat_yes_list[idx][cage_name][0]
                final_e_list = [self.base_e, q_isomer_e, *e_list]
                final_name_list = [self.q_cage.name, cage_name, *name_list]

                e_areas.append(_get_e_area(final_e_list))
                pathways.append(final_pathway)
                e_lists.append(final_e_list)
                name_lists.append(final_name_list)
                return
            else:
                new_name_list = [cage_name, *name_list]
                _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                parent_name_list = flat_yes_list[idx][cage_name][0]
                q_isomer_e = flat_yes_list[idx][cage_name][1]
                new_e_list = [q_isomer_e, *e_list]
                for a_parent in parent_name_list:
                    _, parent_addon_set, _0 = name2seq(cage_size=cage_size, name=a_parent)
                    new_addon_set = addon_set - parent_addon_set
                    new_pathway = [new_addon_set, *pathway]
                    _get_pathway_unit(idx - 1, pathway=new_pathway, cage_name=a_parent, e_list=new_e_list,
                                      name_list=new_name_list)

        def _get_base_e():
            if self.optimizer.mode == 'ase':
                atoms = self.q_cage.atoms
                atoms.calc = self.optimizer.calc
                self.base_e = atoms.get_potential_energy()[0][0]
            elif self.optimizer.mode == 'xtb':
                opted_xyz = os.path.join('temp_opt', self.q_cage.name, 'xtbopt.xyz')
                if not os.path.exists(opted_xyz):
                    os.makedirs(name='temp_raw', exist_ok=True)
                    os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
                    os.makedirs(name='temp_opt', exist_ok=True)
                    self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
                                               cmd_list=self.optimizer.cmd_list)
                with open(opted_xyz, 'r') as f:
                    f.readline()
                    e_line = f.readline()
                    self.base_e = float(e_line.split()[1])
            elif self.optimizer.mode == 'gaussian':
                log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
                if not os.path.exists(log_path):
                    os.makedirs(name='temp_raw', exist_ok=True)
                    os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
                    os.makedirs(name='temp_opt', exist_ok=True)
                    self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
                                               cmd_list=self.optimizer.cmd_list, cycles=self.optimizer.init_cycle)
                log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
                self.last_image = read_gaussian_out(fd=log_path)[-1]
                self.base_e = self.last_image.get_potential_energy()

        def _refine_e(old_info: pd.DataFrame):
            if self.refine_optimizer.mode == 'ase':
                atoms = self.q_cage.atoms
                atoms.calc = self.optimizer.calc
                self.base_e = atoms.get_potential_energy()[0][0]
                return None
            elif self.refine_optimizer.mode == 'xtb':
                os.makedirs(name='refine_workbase', exist_ok=True)
                os.chdir('refine_workbase')
                os.makedirs(name='odd_raw', exist_ok=True)
                os.makedirs(name='even_raw', exist_ok=True)
                os.makedirs(name='opt_workbase', exist_ok=True)

                for add_num_idx, add_num in enumerate(add_num_list[:-1]):
                    addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                    for name_list in old_info['name']:
                        name = name_list[add_num_idx + 1]
                        if add_num % 2 == 1:
                            new_path = os.path.join('odd_raw', f'{name}.xyz')
                        else:
                            new_path = os.path.join('even_raw', f'{name}.xyz')
                        if os.path.exists(new_path):
                            continue
                        if self.optimizer.mode == 'xtb':
                            xyz_filename = 'xtbopt.xyz'
                        elif self.optimizer.mode == 'ase':
                            xyz_filename = 'opt.xyz'
                        elif self.optimizer.mode == 'gaussian':
                            xyz_filename = 'gau.xyz'
                        opt_path = os.path.join(addon_path, f'opt_final', name, xyz_filename)
                        if os.path.exists(opt_path):
                            shutil.copy(src=opt_path, dst=new_path)
                        else:
                            opt_path = os.path.join(addon_path, f'opt_{self.optimizer.init_cycle}', name, xyz_filename)
                            shutil.copy(src=opt_path, dst=new_path)
                self.refine_optimizer.run_a_batch(path_source='even_raw', path_destination='opt_workbase',
                                                  cmd_list=self.refine_optimizer.cmd_list)
                self.refine_optimizer.run_a_batch(path_source='odd_raw', path_destination='opt_workbase',
                                                  cmd_list=[*self.refine_optimizer.cmd_list, '--uhf 1'])

                if self.refine_optimizer.mode == 'xtb':
                    xyz_filename = 'xtbopt.xyz'
                elif self.refine_optimizer.mode == 'ase':
                    xyz_filename = 'opt.xyz'
                elif self.refine_optimizer.mode == 'gaussian':
                    xyz_filename = 'gau.xyz'

                os.chdir('opt_workbase')
                name_e_map = dict()
                for a_name in os.listdir('./'):
                    with open(os.path.join(a_name, xyz_filename), 'r') as f:
                        f.readline()
                        energy_line = f.readline()
                        energy = float(energy_line.split()[1])
                    name_e_map.update({a_name: energy})

                for add_num_idx, add_num in enumerate(add_num_list[:-1]):
                    addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                    for pathway_idx, name_list in enumerate(old_info['name']):
                        name = name_list[add_num_idx + 1]
                        e = name_e_map[name]
                        old_info['e_list'][pathway_idx][add_num_idx + 1] = e
                for idx, e_list in enumerate(old_info['e_list']):
                    e_area = _get_e_area(e_list=e_list)
                    old_info['e_area'][idx] = e_area

                os.chdir(self.workbase)
                refined_info = old_info.sort_values(by='e_area')
                refined_info.index = sorted(refined_info.index)
                return refined_info
            # elif self.optimizer.mode == 'gaussian':
            #     log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            #     if not os.path.exists(log_path):
            #         os.makedirs(name='temp_raw', exist_ok=True)
            #         os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
            #         os.makedirs(name='temp_opt', exist_ok=True)
            #         self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
            #                                    cmd_list=self.optimizer.cmd_list, cycles=self.optimizer.init_cycle)
            #     log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            #     self.last_image = read_gaussian_out(fd=log_path)[-1]
            #     self.base_e = self.last_image.get_potential_energy()

        self.base_e = 0

        # Dump low e isomers to a log
        deep_info_path = os.path.join(self.q_cage.workbase, f'{self.q_add_num}addons', 'deep_yes_info.pickle')
        deep_info = pd.read_pickle(deep_info_path)
        Max_q_rank = len(deep_info['energy'])
        self.log_low_e_num = min(self.log_low_e_num, Max_q_rank)

        log_path = os.path.join(f'top_{self.log_low_e_num}_isomers.log')
        for rank in range(self.log_low_e_num):
            src_path = deep_info['xyz_path'][rank]
            if rank == 0:
                shutil.copy(src=src_path, dst=log_path)
            else:
                with open(log_path, 'a') as w_file, open(src_path, 'r') as r_file:
                    for a_line in r_file.readlines():
                        w_file.write(a_line)

        add_num_list = list(range(self.start, self.q_add_num + self.step, self.step))
        flat_yes_list = []

        if self.ctr_path:
            assert len(add_num_list) > 2, print(
                'Current quary addon number is less than 3 steps away from the begining.\n'
                'Do not need control pathes.')
            for ii, i in enumerate(add_num_list):
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'flat_yes_info.pickle')
                a_flat_yes_info = pd.read_pickle(info_path)
                if ii < 2:
                    flat_yes_list.append(a_flat_yes_info)
                else:
                    parent_flat_info = pd.read_pickle(os.path.join(self.q_cage.workbase, f'{i-self.step}addons', 'flat_yes_info.pickle'))
                    sorted_flat_yes_info = self._sort_parent(flat_info=a_flat_yes_info, info_path=info_path, parent_flat_info=parent_flat_info)
                    flat_yes_list.append(sorted_flat_yes_info)
        else:
            for i in add_num_list:
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'flat_yes_info.pickle')
                a_flat_yes_info = pd.read_pickle(info_path)
                flat_yes_list.append(a_flat_yes_info)

        # Query the path of some low e isomers
        self.q_low_e_num = min(self.q_low_e_num, Max_q_rank)
        q_rank = self.q_low_e_num

        q_isomer_e = deep_info['energy'][q_rank]
        q_isomer_name = deep_info['name'][q_rank]
        q_isomer_atoms = read(deep_info['xyz_path'][q_rank])
        q_rank_workbase = os.path.join(self.workbase, f'isomer_rank_{q_rank}')
        os.makedirs(exist_ok=True, name=q_rank_workbase)
        os.chdir(q_rank_workbase)
        write(filename=q_isomer_name + '.xyz', images=q_isomer_atoms, format='xyz', comment=str(q_isomer_e))
        init_pathway = {0}

        idx = len(add_num_list) - 1
        pathways = []
        e_lists = []
        e_areas = []
        name_lists = []
        cage_size = self.q_cage.size

        # Get path info recursively
        _get_pathway_unit(idx=idx, pathway=[], cage_name=q_isomer_name, e_list=[], name_list=[])
        info = pd.DataFrame({'pathway': pathways, 'name': name_lists, 'e_list': e_lists, 'e_area': e_areas})
        sorted_info = info.sort_values(by='e_area')
        sorted_info.index = sorted(sorted_info.index)

        if self.ctr_path:
            sorted_info = sorted_info[sorted_info.index < self.max_path_num]
        if self.is_refine:
            sorted_info = _refine_e(old_info=sorted_info[sorted_info.index < self.refine_top_k])

        info_path = f'{q_isomer_name}_pathway_info.pickle'
        sorted_info.to_pickle(info_path)

        # Get relatively energy info
        e_lists = sorted_info['e_list']
        e_array = np.array(e_lists[0])
        for i in e_lists[1:]:
            e_array = np.vstack((e_array, np.array(i)))

        if len(e_array.shape) == 1:
            raise RuntimeError(f'There are only one root for rank {q_rank}.')


        rel_e_array = np.ones_like(e_array)
        for i in range(len(add_num_list) + 1):
            rel_e_array[:, i] = e_array[:, i] - min(e_array[:, i])
        rel_e_array = rel_e_array * hartree2kjmol

        # Dump path related xyz file
        for path_rank in range(min(self.q_path_num, len(e_lists))):
            name_list = sorted_info['name'][path_rank]
            one_path_workbase = os.path.join(q_rank_workbase, f'path_rank_{path_rank}')
            os.makedirs(exist_ok=True, name=one_path_workbase)
            os.chdir(one_path_workbase)
            log_path = f'path_rank_{path_rank}_root.log'
            write(filename=log_path, images=self.q_cage.atoms, format='xyz', append=True)
            for idx, add_num in enumerate(add_num_list):
                addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                name = name_list[idx + 1]
                new_path = f'{name}_addon_{add_num}.xyz'
                if self.optimizer.mode == 'xtb':
                    xyz_filename = 'xtbopt.xyz'
                elif self.optimizer.mode == 'ase':
                    xyz_filename = 'opt.xyz'
                elif self.optimizer.mode == 'gaussian':
                    xyz_filename = 'gau.xyz'

                opt_path = os.path.join(addon_path, f'opt_final', name, xyz_filename)
                if os.path.exists(opt_path):
                    shutil.copy(src=opt_path, dst=new_path)
                else:
                    opt_path = os.path.join(addon_path, f'opt_{self.optimizer.init_cycle}', name, xyz_filename)
                    shutil.copy(src=opt_path, dst=new_path)
                atoms = read(new_path, format='xyz')
                write(filename=log_path, images=atoms, format='xyz', append=True)

        np.save(file=os.path.join(q_rank_workbase, f'Path_relative_energy.npy'), arr=rel_e_array)
        # Plot relative energy for path rank
        rel_e_array = rel_e_array[:self.q_path_num, :]
        fig_1 = plt.figure(dpi=400, figsize=(len(add_num_list), min(self.q_path_num, len(e_lists) + 1, 100)))
        cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
        rel_e_df = pd.DataFrame(rel_e_array)
        rel_e_df.columns = [0, *add_num_list]
        sns.heatmap(rel_e_df, annot=True, cmap=cmap, linewidths=.5)
        plt.ylabel('Path rank.')
        plt.xlabel('Addon number.')
        plt.title('Path relative energy (kj/mol).')
        plt.savefig(os.path.join(q_rank_workbase, f'Path_relative_energy.png'))

        # Plot relative energy for isomers.
        isomer_rel_e = np.zeros_like(e_array)
        for idx, add_num in enumerate(add_num_list):
            deep_info_path = os.path.join(self.q_cage.workbase, f'{add_num}addons', 'deep_yes_info.pickle')
            deep_info = pd.read_pickle(deep_info_path)
            isomer_min = deep_info['energy'][0]
            isomer_rel_e[:, idx + 1] = e_array[:, idx + 1] - isomer_min
        isomer_rel_e = isomer_rel_e * hartree2kjmol

        np.save(file=os.path.join(q_rank_workbase, f'Isomer_relative_energy.npy'), arr=isomer_rel_e)
        isomer_rel_e = isomer_rel_e[:self.q_path_num, :]
        # Simple plot
        fig_2 = plt.figure(dpi=400, figsize=(len(add_num_list), min(self.q_path_num, len(e_lists) + 1, 100)))
        cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
        isomer_rel_e_df = pd.DataFrame(isomer_rel_e)
        isomer_rel_e_df.columns = [0, *add_num_list]
        sns.heatmap(isomer_rel_e_df, annot=True, cmap=cmap, linewidths=.5)
        plt.ylabel('Path rank.')
        plt.xlabel('Addon number.')
        plt.title('Isomer relative energy (kj/mol).')
        plt.savefig(os.path.join(q_rank_workbase, f'Isomer_relative_energy.png'))





    def get_path_info_mix(self):
        os.chdir(self.workbase)
        def _get_e_area_triangle(e_list: list):
            e_array = np.array(e_list)
            e_array = e_array - self.base_e
            e_area = e_array[1] / 2
            for i in range(2, len(e_array)):
                e_area += (e_array[i] + e_array[i - 1]) / 2
            return e_area * hartree2kjmol

        def _get_e_area(e_list: list):
            e_area = sum(e_list)
            return e_area

        def _get_pathway_unit(idx: int, pathway: list, cage_name: str, e_list: list, name_list: list):
            if idx == 0:
                _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                final_pathway = [init_pathway, addon_set, *pathway]
                q_isomer_e = flat_yes_list[idx][cage_name][0]
                final_e_list = [self.base_e, q_isomer_e, *e_list]
                final_name_list = [self.q_cage.name, cage_name, *name_list]

                e_areas.append(_get_e_area(final_e_list))
                pathways.append(final_pathway)
                e_lists.append(final_e_list)
                name_lists.append(final_name_list)
                return
            else:
                new_name_list = [cage_name, *name_list]
                _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                parent_name_list = flat_yes_list[idx][cage_name][0]
                q_isomer_e = flat_yes_list[idx][cage_name][1]
                new_e_list = [q_isomer_e, *e_list]
                for a_parent in parent_name_list:
                    _, parent_addon_set, _0 = name2seq(cage_size=cage_size, name=a_parent)
                    new_addon_set = addon_set - parent_addon_set
                    new_pathway = [new_addon_set, *pathway]
                    _get_pathway_unit(idx - 1, pathway=new_pathway, cage_name=a_parent, e_list=new_e_list,
                                      name_list=new_name_list)

        def _get_base_e():
            if self.optimizer.mode == 'ase':
                atoms = self.q_cage.atoms
                atoms.calc = self.optimizer.calc
                self.base_e = atoms.get_potential_energy()[0][0]
            elif self.optimizer.mode == 'xtb':
                opted_xyz = os.path.join('temp_opt', self.q_cage.name, 'xtbopt.xyz')
                if not os.path.exists(opted_xyz):
                    os.makedirs(name='temp_raw', exist_ok=True)
                    os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
                    os.makedirs(name='temp_opt', exist_ok=True)
                    self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
                                               cmd_list=self.optimizer.cmd_list)
                with open(opted_xyz, 'r') as f:
                    f.readline()
                    e_line = f.readline()
                    self.base_e = float(e_line.split()[1])
            elif self.optimizer.mode == 'gaussian':
                log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
                if not os.path.exists(log_path):
                    os.makedirs(name='temp_raw', exist_ok=True)
                    os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
                    os.makedirs(name='temp_opt', exist_ok=True)
                    self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
                                               cmd_list=self.optimizer.cmd_list, cycles=self.optimizer.init_cycle)
                log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
                self.last_image = read_gaussian_out(fd=log_path)[-1]
                self.base_e = self.last_image.get_potential_energy()

        def _refine_e(old_info: pd.DataFrame):
            if self.refine_optimizer.mode == 'ase':
                atoms = self.q_cage.atoms
                atoms.calc = self.optimizer.calc
                self.base_e = atoms.get_potential_energy()[0][0]
                return None
            elif self.refine_optimizer.mode == 'xtb':
                os.makedirs(name='refine_workbase', exist_ok=True)
                os.chdir('refine_workbase')
                os.makedirs(name='odd_raw', exist_ok=True)
                os.makedirs(name='even_raw', exist_ok=True)
                os.makedirs(name='opt_workbase', exist_ok=True)

                for add_num_idx, add_num in enumerate(add_num_list[:-1]):
                    addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                    for name_list in old_info['name']:
                        name = name_list[add_num_idx + 1]
                        if add_num % 2 == 1:
                            new_path = os.path.join('odd_raw', f'{name}.xyz')
                        else:
                            new_path = os.path.join('even_raw', f'{name}.xyz')
                        if os.path.exists(new_path):
                            continue
                        if self.optimizer.mode == 'xtb':
                            xyz_filename = 'xtbopt.xyz'
                        elif self.optimizer.mode == 'ase':
                            xyz_filename = 'opt.xyz'
                        elif self.optimizer.mode == 'gaussian':
                            xyz_filename = 'gau.xyz'
                        opt_path = os.path.join(addon_path, f'opt_final', name, xyz_filename)
                        if os.path.exists(opt_path):
                            shutil.copy(src=opt_path, dst=new_path)
                        else:
                            opt_path = os.path.join(addon_path, f'opt_{self.optimizer.init_cycle}', name, xyz_filename)
                            shutil.copy(src=opt_path, dst=new_path)
                self.refine_optimizer.run_a_batch(path_source='even_raw', path_destination='opt_workbase',
                                                  cmd_list=self.refine_optimizer.cmd_list)
                self.refine_optimizer.run_a_batch(path_source='odd_raw', path_destination='opt_workbase',
                                                  cmd_list=[*self.refine_optimizer.cmd_list, '--uhf 1'])

                if self.refine_optimizer.mode == 'xtb':
                    xyz_filename = 'xtbopt.xyz'
                elif self.refine_optimizer.mode == 'ase':
                    xyz_filename = 'opt.xyz'
                elif self.refine_optimizer.mode == 'gaussian':
                    xyz_filename = 'gau.xyz'

                os.chdir('opt_workbase')
                name_e_map = dict()
                for a_name in os.listdir('./'):
                    with open(os.path.join(a_name, xyz_filename), 'r') as f:
                        f.readline()
                        energy_line = f.readline()
                        energy = float(energy_line.split()[1])
                    name_e_map.update({a_name: energy})

                for add_num_idx, add_num in enumerate(add_num_list[:-1]):
                    addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                    for pathway_idx, name_list in enumerate(old_info['name']):
                        name = name_list[add_num_idx + 1]
                        e = name_e_map[name]
                        old_info['e_list'][pathway_idx][add_num_idx + 1] = e
                for idx, e_list in enumerate(old_info['e_list']):
                    e_area = _get_e_area(e_list=e_list)
                    old_info['e_area'][idx] = e_area

                os.chdir(self.workbase)
                refined_info = old_info.sort_values(by='e_area')
                refined_info.index = sorted(refined_info.index)
                return refined_info
            # elif self.optimizer.mode == 'gaussian':
            #     log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            #     if not os.path.exists(log_path):
            #         os.makedirs(name='temp_raw', exist_ok=True)
            #         os.symlink(src=self.q_cage.pristine_path, dst=os.path.join('temp_raw', self.q_cage.name + '.xyz'))
            #         os.makedirs(name='temp_opt', exist_ok=True)
            #         self.optimizer.run_a_batch(path_source='temp_raw', path_destination='temp_opt',
            #                                    cmd_list=self.optimizer.cmd_list, cycles=self.optimizer.init_cycle)
            #     log_path = os.path.join('temp_opt', self.q_cage.name, 'gau.log')
            #     self.last_image = read_gaussian_out(fd=log_path)[-1]
            #     self.base_e = self.last_image.get_potential_energy()

        self.base_e = 0

        # Dump low e isomers to a log
        deep_info_path = os.path.join(self.q_cage.workbase, f'{self.q_add_num}addons', 'deep_yes_info.pickle')
        deep_info = pd.read_pickle(deep_info_path)
        Max_q_rank = len(deep_info['energy'])
        self.log_low_e_num = min(self.log_low_e_num, Max_q_rank)

        log_path = f'top_{self.log_low_e_num}_isomers.log'
        for rank in range(self.log_low_e_num):
            src_path = deep_info['xyz_path'][rank]
            if rank == 0:
                shutil.copy(src=src_path, dst=log_path)
            else:
                with open(log_path, 'a') as w_file, open(src_path, 'r') as r_file:
                    for a_line in r_file.readlines():
                        w_file.write(a_line)

        add_num_list = list(range(self.start, self.q_add_num + self.step, self.step))
        flat_yes_list = []

        if self.ctr_path:
            assert len(add_num_list) > 2, print(
                'Current quary addon number is less than 3 steps away from the begining.\n'
                'Do not need control pathes.')
            for ii, i in enumerate(add_num_list):
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'flat_yes_info.pickle')
                a_flat_yes_info = pd.read_pickle(info_path)
                if ii < 2:
                    flat_yes_list.append(a_flat_yes_info)
                else:
                    parent_flat_info = pd.read_pickle(os.path.join(self.q_cage.workbase, f'{i-self.step}addons', 'flat_yes_info.pickle'))
                    sorted_flat_yes_info = self._sort_parent(flat_info=a_flat_yes_info, info_path=info_path, parent_flat_info=parent_flat_info)
                    flat_yes_list.append(sorted_flat_yes_info)
        else:
            for i in add_num_list:
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'flat_yes_info.pickle')
                a_flat_yes_info = pd.read_pickle(info_path)
                flat_yes_list.append(a_flat_yes_info)

        # Query the path of some low e isomers
        self.q_low_e_num = min(self.q_low_e_num, Max_q_rank)
        all_pathways = []
        all_e_lists = []
        all_e_areas = []
        all_name_lists = []
        isomer_dump = os.path.join(self.workbase, 'isomer_dump')
        os.makedirs(isomer_dump, exist_ok=True)
        os.chdir(isomer_dump)
        for q_rank in range(self.q_low_e_num):
            q_isomer_e = deep_info['energy'][q_rank]
            q_isomer_name = deep_info['name'][q_rank]
            q_isomer_atoms = read(deep_info['xyz_path'][q_rank])

            write(filename=q_isomer_name + f'_rank_{q_rank}.xyz', images=q_isomer_atoms, format='xyz', comment=str(q_isomer_e))
            init_pathway = {0}

            idx = len(add_num_list) - 1
            pathways = []
            e_lists = []
            e_areas = []
            name_lists = []

            cage_size = self.q_cage.size

            # Get path info recursively
            _get_pathway_unit(idx=idx, pathway=[], cage_name=q_isomer_name, e_list=[], name_list=[])
            all_pathways.extend(pathways)
            all_e_lists.extend(e_lists)
            all_e_areas.extend(e_areas)
            all_name_lists.extend(name_lists)

        os.chdir(self.workbase)
        info = pd.DataFrame({'pathway': all_pathways, 'name': all_name_lists, 'e_list': all_e_lists, 'e_area': all_e_areas})
        sorted_info = info.sort_values(by='e_area')
        sorted_info.index = sorted(sorted_info.index)
        sorted_info.to_pickle(f'all_info.pickle')
        if self.ctr_path:
            sorted_info = sorted_info[sorted_info.index < self.max_path_num]
        if self.is_refine:
            sorted_info = _refine_e(old_info=sorted_info[sorted_info.index < self.refine_top_k])

        sorted_info.to_pickle(f'refined_info.pickle')

        # Get relatively energy info
        e_lists = sorted_info['e_list']
        e_array = np.array(e_lists[0])
        for i in e_lists[1:]:
            e_array = np.vstack((e_array, np.array(i)))

        # if len(e_array.shape) == 1:
        #     print(f'There are only one root for rank {q_rank}.')
        #     continue

        rel_e_array = np.ones_like(e_array)
        for i in range(len(add_num_list) + 1):
            rel_e_array[:, i] = e_array[:, i] - min(e_array[:, i])
        rel_e_array = rel_e_array * hartree2kjmol

        # Dump path related xyz file
        for path_rank in range(min(self.q_path_num, len(e_lists))):
            name_list = sorted_info['name'][path_rank]
            one_path_workbase = os.path.join(self.workbase, 'path_dump', f'path_rank_{path_rank}')
            os.makedirs(exist_ok=True, name=one_path_workbase)
            os.chdir(one_path_workbase)
            log_path = f'path_rank_{path_rank}_root.log'
            write(filename=log_path, images=self.q_cage.atoms, format='xyz', append=True)
            for idx, add_num in enumerate(add_num_list):
                addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                name = name_list[idx + 1]
                new_path = f'{name}_addon_{add_num}.xyz'
                if self.optimizer.mode == 'xtb':
                    xyz_filename = 'xtbopt.xyz'
                elif self.optimizer.mode == 'ase':
                    xyz_filename = 'opt.xyz'
                elif self.optimizer.mode == 'gaussian':
                    xyz_filename = 'gau.xyz'

                opt_path = os.path.join(addon_path, f'opt_final', name, xyz_filename)
                if os.path.exists(opt_path):
                    shutil.copy(src=opt_path, dst=new_path)
                else:
                    opt_path = os.path.join(addon_path, f'opt_{self.optimizer.init_cycle}', name, xyz_filename)
                    shutil.copy(src=opt_path, dst=new_path)
                atoms = read(new_path, format='xyz')
                write(filename=log_path, images=atoms, format='xyz', append=True)
        os.chdir(self.workbase)
        np.save(file=f'Path_relative_energy.npy', arr=rel_e_array)
        # Plot relative energy for path rank
        rel_e_array = rel_e_array[:self.q_path_num, :]
        fig_1 = plt.figure(dpi=400, figsize=(len(add_num_list), min(self.q_path_num, len(e_lists) + 1, 20)))
        cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
        rel_e_df = pd.DataFrame(rel_e_array)
        rel_e_df.columns = [0, *add_num_list]
        sns.heatmap(rel_e_df, annot=True, cmap=cmap, linewidths=.5)
        plt.ylabel('Path rank.')
        plt.xlabel('Addon number.')
        plt.title('Path relative energy (kj/mol).')
        plt.savefig(f'Path_relative_energy.png')

        # Plot relative energy for isomers.
        isomer_rel_e = np.zeros_like(e_array)
        for idx, add_num in enumerate(add_num_list):
            deep_info_path = os.path.join(self.q_cage.workbase, f'{add_num}addons', 'deep_yes_info.pickle')
            deep_info = pd.read_pickle(deep_info_path)
            isomer_min = deep_info['energy'][0]
            isomer_rel_e[:, idx + 1] = e_array[:, idx + 1] - isomer_min
        isomer_rel_e = isomer_rel_e * hartree2kjmol

        np.save(file=f'Isomer_relative_energy.npy', arr=isomer_rel_e)
        isomer_rel_e = isomer_rel_e[:self.q_path_num, :]
        # Simple plot
        fig_2 = plt.figure(dpi=400, figsize=(len(add_num_list), min(self.q_path_num, len(e_lists) + 1, 20)))
        cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
        isomer_rel_e_df = pd.DataFrame(isomer_rel_e)
        isomer_rel_e_df.columns = [0, *add_num_list]
        sns.heatmap(isomer_rel_e_df, annot=True, cmap=cmap, linewidths=.5)
        plt.ylabel('Path rank.')
        plt.xlabel('Addon number.')
        plt.title('Isomer relative energy (kj/mol).')
        plt.savefig(f'Isomer_relative_energy.png')
