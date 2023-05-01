import os

import matplotlib.pyplot as plt
import networkx as nx
from ase import Atoms
from autosteper.cage import Cage, seq2name, name2seq
from autosteper.optimizers import *
from autosteper.tools import get_low_e_ranks, strip_extraFullerene, get_G, simple_dump_pathway, get_pathway_info, \
    get_connection, simple_parse_logs, match_swr_cages, lazy_file_mid_name, re_label_cages, plot_pathway_heatmap
from networkx import isomorphism
from tqdm import tqdm


def refine(old_workbase: str, new_workbase: str, ref_para: dict):
    refiner = switch_optimizers(mode=ref_para['opt_mode'], para=ref_para['opt_para'])
    add_num_list = range(ref_para['start'], ref_para['stop'] + ref_para['step'], ref_para['step'])
    for an_add_num in add_num_list:
        os.chdir(old_workbase)
        addon_path = str(an_add_num) + 'addons'
        os.chdir(addon_path)
        passed_info = pd.read_pickle('passed_info.pickle')
        xyz_paths = passed_info['xyz_path']
        e_arr = np.array(passed_info['energy'])
        os.chdir(new_workbase)
        os.makedirs(addon_path)
        os.chdir(addon_path)
        refiner.set_folders()
        os.chdir('raw')
        for a_rank in get_low_e_ranks(e_arr=e_arr, para=ref_para['cutoff']):
            an_old_path = xyz_paths[a_rank]
            a_name = os.path.basename(an_old_path)
            shutil.copy(src=an_old_path, dst=a_name)
        os.chdir('./..')
        ref_code = refiner.opt()
        if ref_code != 0:
            raise RuntimeError(f'Something wrong happened in {addon_path}.\n'
                               f'Refinement procedure terminated.')
    print('Normal termination of refinement procedure.')


def simple_test_iso(q_atoms: Atoms, passed_info_path: str, top_k: int = None):
    print('Please make sure the queried atoms has been optimized with at least semi-empirical level of theory.')
    target_G = get_G(q_atoms)
    passed_info = pd.read_pickle(passed_info_path)
    flag = 0
    for idx, a_xyz_path in enumerate(tqdm(passed_info['xyz_path'][:top_k])):
        an_atoms = read(a_xyz_path)
        pred_G = get_G(an_atoms)
        if nx.is_isomorphic(pred_G, target_G):
            rel_e = passed_info['energy'][idx] - passed_info['energy'][0]
            print('The queried atoms has been found in AutoStepper scan results.\n'
                  f'The target is in rank {str(idx + 1)}.\n'
                  f'Has a relative energy of {str(rel_e)} eV considering the lowest isomer in scan results.\n')
            flag = 1
            break
    if flag == 0:
        if top_k:
            print(f'The queried atoms has NOT been found in top {top_k} AutoStepper scan results.\n'
                  'Enlarge test scope may help.')
        else:
            print('The queried atoms has NOT been found in AutoStepper scan results.')
        return None
    else:
        return idx + 1


def _prep_relatives(q_atoms: Atoms = None, group: str = None, cage_size: int = None, q_seq: str = None,
                    q_cage: Cage = None):
    if q_atoms:
        print('Please make sure the queried atoms has an identical cage with the scanned cage.\n'
              'The atom sequence in two cages should be identical.\n'
              'The addons should be atomic addons.\n'
              'Complex senerio could be tackled with strip_extraFullerene function in tools module.')
        pristine_cage, q_add_set = strip_extraFullerene(atoms=q_atoms, group=group)
        cage_size = len(pristine_cage)
        mid_add_num = len(q_add_set)
    elif q_seq:
        _, q_add_set, _0 = seq2name(q_seq, q_cage)
        mid_add_num = len(q_add_set)
        cage_size = q_cage.size
    else:
        print('Please input query addon set.\n'
              'Currently only surpport ase.atoms input, sequence input and addon set input.')
    return q_add_set, cage_size, mid_add_num


def simple_log_relatives(workbase: str, dump_log_path: str,
                         step: int, fst_add_num: int = None, final_add_num: int = None,
                         q_atoms: Atoms = None, group: str = None, q_seq: str = None, q_cage: Cage = None):
    q_add_set, cage_size, mid_add_num = _prep_relatives(q_atoms=q_atoms, group=group, q_seq=q_seq, q_cage=q_cage)
    if fst_add_num:
        for an_add_num in range(fst_add_num, mid_add_num, step):
            addon_path = str(an_add_num) + 'addons'
            os.chdir(workbase)
            os.chdir(addon_path)
            with open(dump_log_path, 'a') as f:
                passed_info = pd.read_pickle('passed_info.pickle')
                min_e = passed_info['energy'][0]
                f.write(f'Start screen SubStructures in {addon_path}:\n')
                flag = 0
                for idx, a_name in enumerate(passed_info['name']):
                    _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
                    if len(pred_add_set - q_add_set) == 0:
                        f.write('      A SubStructure has been found:\n')
                        f.write(f'      The absent addon sites are:      {str(q_add_set - pred_add_set)}.\n')
                        f.write(f'      In rank:                         {str(idx + 1)}\n')
                        rel_e = passed_info['energy'][idx] - min_e
                        f.write(f'      Relative energy is:              {str(rel_e)}eV\n')
                        f.write(f'      Name is:                         {a_name}\n\n\n')
                        flag = 1
                if flag == 0:
                    f.write(f'{addon_path} do not have a SubStructure.\n')
                f.write(
                    '\n############################################################################################\n')

    mid_passed_info_path = os.path.join(workbase, f'{mid_add_num}addons', 'passed_info.pickle')
    mid_rank = simple_test_iso(q_atoms=q_atoms, passed_info_path=mid_passed_info_path)
    mid_rank = mid_rank - 1
    passed_info = pd.read_pickle(mid_passed_info_path)
    a_name = passed_info['name'][mid_rank]
    rel_e = passed_info['energy'][mid_rank] - passed_info['energy'][0]
    with open(dump_log_path, 'a') as f:
        f.write(f'\n\nStart screen Queried isomer:\n')
        if mid_rank:
            f.write('      Queried isomer has been found:\n')
            f.write(f'      In rank:                         {str(mid_rank + 1)}\n')
            f.write(f'      Relative energy is:              {str(rel_e)}eV\n')
            f.write(f'      Name is:                         {a_name}\n\n\n')
        else:
            f.write('Queried isomer has not been found:\n\n\n')
        f.write('\n############################################################################################\n')

    if final_add_num:
        for an_add_num in range(mid_add_num + step, final_add_num + step, step):
            addon_path = str(an_add_num) + 'addons'
            os.chdir(workbase)
            os.chdir(addon_path)
            with open(dump_log_path, 'a') as f:
                passed_info = pd.read_pickle('passed_info.pickle')
                min_e = passed_info['energy'][0]
                f.write(f'Start screen SuperStructures in {addon_path}:\n')
                flag = 0
                for idx, a_name in enumerate(passed_info['name']):
                    _, pred_add_set, _0 = name2seq(name=a_name, cage_size=cage_size)
                    if len(q_add_set - pred_add_set) == 0:
                        f.write('      A SuperStructure has been found:\n')
                        f.write(f'      The new addon sites are:      {str(pred_add_set - q_add_set)}.\n')
                        f.write(f'      In rank:                         {str(idx + 1)}\n')
                        rel_e = passed_info['energy'][idx] - min_e
                        f.write(f'      Relative energy is:              {str(rel_e)}eV\n')
                        f.write(f'      Name is:                         {a_name}\n\n\n')
                        flag = 1
                if flag == 0:
                    f.write(f'{addon_path} do not have a SubStructure.\n')
                f.write(
                    '\n############################################################################################\n')


def strict_scatter_relatives(workbase: str, dump_folder: str,
                             step: int, fst_add_num: int = None, final_add_num: int = None,
                             q_atoms: Atoms = None, q_seq: str = None, q_cage: Cage = None,
                             cutoff: dict = None, group: str = None):
    is_first = True

    def _prep_unit():
        addon_path = str(an_add_num) + 'addons'
        os.chdir(workbase)
        os.chdir(addon_path)
        passed_info = pd.read_pickle('passed_info.pickle')
        e_arr = np.array(passed_info['energy'])
        e_arr = e_arr - e_arr[0]
        names = passed_info['name']
        xyz_paths = passed_info['xyz_path']
        return e_arr, names, xyz_paths

    q_add_set, cage_size, mid_add_num = _prep_relatives(q_atoms=q_atoms, group=group, q_seq=q_seq, q_cage=q_cage)
    q_G = get_G(q_atoms)
    fig = plt.figure(dpi=400)
    fig_len = 0
    add_num_list = []
    if fst_add_num:
        for an_add_num in range(fst_add_num, mid_add_num + step, step):
            add_num_list.append(an_add_num)
            e_arr, names, xyz_paths = _prep_unit()
            non_rel_e_list = []
            rel_e_list = []
            for a_rank in get_low_e_ranks(e_arr=e_arr, para=cutoff):
                _, pred_add_set, _0 = name2seq(name=names[a_rank], cage_size=cage_size)
                if len(pred_add_set - q_add_set) != 0:
                    pred_atoms = read(xyz_paths[a_rank])
                    pred_G = get_G(pred_atoms)
                    GM = isomorphism.GraphMatcher(q_G, pred_G)
                    if not GM.subgraph_is_isomorphic():
                        non_rel_e_list.append(e_arr[a_rank])
                        continue
                rel_e_list.append(e_arr[a_rank])
            os.chdir(dump_folder)
            np.save(file=f'{an_add_num}_addons_rel_e.npy', arr=np.array(rel_e_list))
            np.save(file=f'{an_add_num}_addons_non_rel_e.npy', arr=np.array(non_rel_e_list))
            if is_first:
                plt.scatter(x=[an_add_num] * len(non_rel_e_list), y=non_rel_e_list, marker='+', c='blue',
                            label='non-rel')
                plt.scatter(x=[an_add_num] * len(rel_e_list), y=rel_e_list, marker='x', c='red', label='rel')
                is_first = False
            else:
                plt.scatter(x=[an_add_num] * len(non_rel_e_list), y=non_rel_e_list, marker='+', c='blue')
                plt.scatter(x=[an_add_num] * len(rel_e_list), y=rel_e_list, marker='x', c='red')

    if final_add_num:
        for an_add_num in range(mid_add_num + step, final_add_num + step, step):
            add_num_list.append(an_add_num)
            e_arr, names, xyz_paths = _prep_unit()
            non_rel_e_list = []
            rel_e_list = []
            for a_rank in get_low_e_ranks(e_arr=e_arr, para=cutoff):
                _, pred_add_set, _0 = name2seq(name=names[a_rank], cage_size=cage_size)
                if len(q_add_set - pred_add_set) != 0:
                    pred_atoms = read(xyz_paths[a_rank])
                    pred_G = get_G(pred_atoms)
                    GM = isomorphism.GraphMatcher(pred_G, q_G)
                    if not GM.subgraph_is_isomorphic():
                        non_rel_e_list.append(e_arr[a_rank])
                        continue
                rel_e_list.append(e_arr[a_rank])
            os.chdir(dump_folder)
            np.save(file=f'{an_add_num}_addons_rel_e.npy', arr=np.array(rel_e_list))
            np.save(file=f'{an_add_num}_addons_non_rel_e.npy', arr=np.array(non_rel_e_list))
            if is_first:
                plt.scatter(x=[an_add_num] * len(non_rel_e_list), y=non_rel_e_list, marker='+', c='blue',
                            label='non-rel')
                plt.scatter(x=[an_add_num] * len(rel_e_list), y=rel_e_list, marker='x', c='red', label='rel')
                is_first = False
            else:
                plt.scatter(x=[an_add_num] * len(non_rel_e_list), y=non_rel_e_list, marker='+', c='blue')
                plt.scatter(x=[an_add_num] * len(rel_e_list), y=rel_e_list, marker='x', c='red')

    plt.ylim(None, 1)
    plt.legend(loc='best')
    plt.xticks(ticks=add_num_list, labels=add_num_list)
    plt.xlabel('Number of addends')
    plt.ylabel('Relative energy (eV)')
    plt.savefig('rel_e.png', figsize=(len(add_num_list), 7))


def cook_disordered(disordered_root: str, dump_root: str, step: int, log_mode: str, keep_top_k_pathway: int,
                    group: str = None, file_mid_name: str = None, dpi: int = 400, has_pathway: bool=True):
    """
    Pipeline to cook disordered log files.
    Sort log files first, get connection and pathway infos, dump pathways.
    """
    disordered_root = os.path.abspath(disordered_root)
    dump_root = os.path.abspath(dump_root)
    sorted_root = os.path.join(dump_root, 'sorted')
    simple_parse_logs(dump_root=sorted_root, src_root=disordered_root,
                      mode=log_mode, group=group, file_mid_name=file_mid_name)
    get_connection(xyz_root=os.path.join(sorted_root, 'xyz'),
                   connection_dump=os.path.join(sorted_root, 'connection'),
                   step=step)
    if has_pathway:
        get_pathway_info(e_info_path=os.path.join(sorted_root, 'info', 'info.pickle'),
                         xyz_root=os.path.join(sorted_root, 'xyz'),
                         cnt_root=os.path.join(sorted_root, 'connection'),
                         dump_info_path=os.path.join(dump_root, r'pathway_info.pickle'))
        simple_dump_pathway(pathway_info_path=os.path.join(dump_root, r'pathway_info.pickle'),
                            dump_root=os.path.join(dump_root, r'pathways'),
                            src_xyz_root=os.path.join(sorted_root, 'xyz'),
                            top_K=keep_top_k_pathway, dpi=dpi)


def find_SWR(q_sorted_root: str, tgt_sorted_root: str, swr_dump_path: str,
             is_low_e: bool = True, is_unique: bool = True,
             group_symbol: str = None, file_mid_name: str = None):
    '''
    Find SWR from atoms in q_logs to atoms in tgt_logs.
    Pristine cage needs to be prepared in log root.

    is_low_e: set True to apply the energy criterion,
    is_unique: set True to keep one SWR pair for each queried atoms if there are many target atoms met requirements,
    '''
    cwd_swr = os.getcwd()
    q_sorted_root = os.path.abspath(q_sorted_root)
    tgt_sorted_root = os.path.abspath(tgt_sorted_root)
    rel_file_mid_name = lazy_file_mid_name(file_mid_name)
    swr_dump_path = os.path.abspath(swr_dump_path)
    os.makedirs(swr_dump_path, exist_ok=True)

    q_xyz_root = os.path.join(q_sorted_root, 'xyz')
    q_cage = read(os.path.join(q_xyz_root, f'0_addons{rel_file_mid_name}1.xyz'))
    q_G = get_G(q_cage)
    q_cnt_root = os.path.join(q_sorted_root, 'connection')

    tgt_xyz_root = os.path.join(tgt_sorted_root, 'xyz')
    tgt_cage = read(os.path.join(tgt_xyz_root, f'0_addons{rel_file_mid_name}1.xyz'))
    tgt_G = get_G(tgt_cage)
    tgt_cnt_root = os.path.join(tgt_sorted_root, 'connection')

    match_map_list, swr_pair_list, q_swr_adj_nodes_set_list = match_swr_cages(query_cage=q_cage, target_cage=tgt_cage)

    chked_q_add = {}
    max_rank = 0

    for a_file in os.listdir(q_xyz_root):
        a_name = os.path.splitext(a_file)[0]
        a_rank = int(a_name.split('_')[-1])
        max_rank = max(max_rank, a_rank)
        an_add_num = int(a_name.split('_')[0])
        if an_add_num not in chked_q_add.keys():
            chked_q_add.update({an_add_num: 0})
    for an_add_num in chked_q_add.keys():
        for a_rank in range(max_rank):
            a_file = f'{an_add_num}_addons{rel_file_mid_name}{a_rank + 1}.xyz'
            a_file_path = os.path.join(q_xyz_root, a_file)
            if not os.path.exists(a_file_path):
                continue
            a_q_atoms = read(a_file_path)
            _, q_addon_set = strip_extraFullerene(atoms=a_q_atoms, group=group_symbol)
            a_q_G = get_G(a_q_atoms)
            tgt_atoms_list = []
            rank_list = []
            tgt_sites_list = []
            tgt_name_list = []
            swr_pair_sites_list = []
            # Match target atoms
            for idx, a_swr_pair in enumerate(swr_pair_list):
                if len(set(a_swr_pair[0]) & q_addon_set) > 0:
                    continue
                a_q_swr_adj_nodes_set = q_swr_adj_nodes_set_list[idx]
                ept_nodes = a_q_swr_adj_nodes_set - q_addon_set
                if not len(ept_nodes) > 0:
                    continue
                a_q_G_for_chk = a_q_G.copy()
                a_q_G_for_chk.remove_nodes_from(a_swr_pair[0])
                for a_tgt_file in os.listdir(tgt_xyz_root):
                    a_tgt_rank = int(a_tgt_file.split('.')[0].split('_')[-1]) - 1
                    if a_tgt_rank in rank_list:
                        continue
                    tgt_add_num = int(a_tgt_file.split('_')[0])
                    if tgt_add_num == len(q_addon_set) + 2:
                        a_tgt_atoms = read(os.path.join(tgt_xyz_root, a_tgt_file))
                        _, tgt_addon_set = strip_extraFullerene(atoms=a_tgt_atoms, group=group_symbol)
                        if len(tgt_addon_set & set(a_swr_pair[1])) > 0:
                            continue
                        a_tgt_G_for_chk = get_G(a_tgt_atoms)
                        a_tgt_G_for_chk.remove_nodes_from(a_swr_pair[1])
                        a_GM = isomorphism.GraphMatcher(a_tgt_G_for_chk, a_q_G_for_chk)
                        if not a_GM.subgraph_is_isomorphic():
                            continue
                        a_match_map = match_map_list[idx]
                        mapped_tgt_addon_set = set()
                        for a_tgt_addon in tgt_addon_set:
                            mapped_tgt_addon_set.add(a_match_map[a_tgt_addon])
                        if len(mapped_tgt_addon_set & ept_nodes) == 0:
                            continue
                        # A SWR pair has been found
                        rank_list.append(a_tgt_rank)
                        new_target_atoms = re_label_cages(old_atoms=a_tgt_atoms, is_reverse_map=True,
                                                          match_map=a_match_map, swr_pair=a_swr_pair)
                        tgt_atoms_list.append(new_target_atoms)
                        tgt_sites_list.append(mapped_tgt_addon_set)
                        tgt_name_list.append(a_tgt_file)
                        swr_pair_sites_list.append(a_swr_pair[0])
            if len(rank_list) == 0:
                continue
            rank_atoms_map = dict(zip(rank_list, tgt_atoms_list))
            rank_sites_map = dict(zip(rank_list, tgt_sites_list))
            rank_swr_map = dict(zip(rank_list, swr_pair_sites_list))
            rank_name_map = dict(zip(rank_list, tgt_name_list))
            rank_list = sorted(rank_list)
            # Applicate the low e criterion
            if is_low_e:
                tgt_info = pd.read_pickle(os.path.join(tgt_sorted_root, 'info', 'info.pickle'))
                q_info = pd.read_pickle(os.path.join(q_sorted_root, 'info', 'info.pickle'))
                tgt_e_list = []
                for a_tgt_rank in rank_list:
                    tgt_e_list.append(tgt_info[len(q_addon_set) + 2][a_tgt_rank])
                cnt = np.load(os.path.join(q_cnt_root, a_file[:-4] + '.npy'))
                new_q_info = []
                for idx, a_q_e in enumerate(q_info[len(q_addon_set) + 2].dropna()):
                    if not a_q_e == None:
                        if cnt[idx] == 1:
                            new_q_info.append(a_q_e)
                if len(new_q_info) != 0:
                    new_rank_list = []
                    for idx, a_tgt_e in enumerate(tgt_e_list):
                        if a_tgt_e < min(new_q_info):
                            new_rank_list.append(rank_list[idx])
                    if len(new_rank_list) == 0:
                        continue
                    rank_list = new_rank_list
            ###################################################
            # dump
            os.chdir(swr_dump_path)
            os.makedirs(f'{len(q_addon_set)}_to_{len(tgt_addon_set)}_swr_{chked_q_add[len(q_addon_set)] + 1}',
                        exist_ok=True)
            os.chdir(f'{len(q_addon_set)}_to_{len(tgt_addon_set)}_swr_{chked_q_add[len(q_addon_set)] + 1}')
            chked_q_add[len(q_addon_set)] += 1
            for idx, a_rank in enumerate(rank_list):
                if is_unique and idx > 0:
                    break
                write(filename='q_atoms.xyz', images=a_q_atoms, format='xyz')
                write(filename=f'tgt_atoms_rank_{idx + 1}.xyz', images=rank_atoms_map[a_rank], format='xyz')
                with open('sites_info.txt', 'a') as f_w:
                    f_w.write('=========A SWR pair has been found.======\n\n'
                              f'The query atoms are in file q_atoms.xyz.\n'
                              f'The target atoms are in file tgt_atoms_rank_{idx + 1}.xyz.\n\n'
                              f'The addon sites in query atoms are:\n{str(sorted([x + 1 for x in q_addon_set]))}\n\n'
                              f'The addon sites in target atoms are:\n'
                              f'{str(sorted([x + 1 for x in rank_sites_map[a_rank]]))}\n\n'
                              f'The SWR bond in query atoms are: {str([x + 1 for x in rank_swr_map[a_rank]])}\n\n'
                              f'The SWR bond in target atoms are: {str([x + 1 for x in rank_swr_map[a_rank]])}\n\n'
                              f'========================================\n\n\n')
    os.chdir(cwd_swr)


def get_binding_e(sorted_root: str, addends_e: float, cage_e: float, file_mid_name: str = None):
    cwd_binding = os.getcwd()
    sorted_root = os.path.abspath(sorted_root)
    os.chdir(sorted_root)
    isomer_e_info = pd.read_pickle(os.path.join('info', 'info.pickle'))
    first_add = isomer_e_info.columns[0]
    step = first_add
    max_rank = len(isomer_e_info[first_add])
    get_connection(xyz_root='./xyz',
                   connection_dump='./connection',
                   step=step,
                   file_mid_name=file_mid_name)
    bind_e_info = {}
    a_e_list = []
    for a_e in isomer_e_info[first_add].dropna():
        a_e_list.append(a_e - cage_e - addends_e)
    bind_e_info[first_add] = a_e_list
    os.chdir('connection')
    for an_add in isomer_e_info.columns[1:]:
        a_e_list = []
        for a_rank, a_e in enumerate(isomer_e_info[an_add].dropna()):
            rank_list = []
            prev_add = an_add - step
            for a_cnt_f in os.listdir('./'):
                if a_cnt_f.startswith(f'{prev_add}_'):
                    a_cnt = np.load(a_cnt_f)
                    if a_cnt[a_rank] == 1:
                        a_prev_rank = int(os.path.splitext(a_cnt_f)[0].split('_')[-1]) - 1
                        rank_list.append(a_prev_rank)
            if len(rank_list) == 0:
                continue
            else:
                lowest_e_rank = sorted(rank_list)[0]
                a_prev_e = isomer_e_info[prev_add][lowest_e_rank]
                a_e_list.append(a_e - a_prev_e - addends_e)
        bind_e_info[an_add] = a_e_list
    for an_add in isomer_e_info.columns:
        a_e_list = bind_e_info[an_add]
        for i in range(max_rank - len(a_e_list)):
            a_e_list.append(None)
        isomer_e_info[an_add] = a_e_list
    os.chdir(sorted_root)
    os.chdir('info')
    isomer_e_info.to_pickle('binding_e.pickle')
    isomer_e_info.to_excel('binding_e.xlsx')


def count_swr_unit(swr_out_root: str):
    if len(os.listdir(swr_out_root)) == 0:
        return 0
    cwd_swr_count = os.getcwd()
    swr_counter = {}
    for a_swr in os.listdir(swr_out_root):
        a_q_add_num = int(a_swr.split('_')[0])
        a_swr_root = os.path.join(swr_out_root, a_swr)
        xyz_counter = 0
        for a_file in os.listdir(a_swr_root):
            if a_file.endswith('xyz'):
                xyz_counter = xyz_counter + 1
        tgt_xyz_counter = xyz_counter - 1
        if a_q_add_num in swr_counter.keys():
            tgt_xyz_counter = tgt_xyz_counter + swr_counter[a_q_add_num]
        swr_counter.update({a_q_add_num: tgt_xyz_counter})
    os.chdir(cwd_swr_count)
    return swr_counter


def count_SWR(swr_1_workbase: str, swr_2_workbase: str, swr_1_legend: str, swr_2_legend: str,
              dump_pic_path: str = None):
    swr_1 = count_swr_unit(swr_1_workbase)
    swr_2 = count_swr_unit(swr_2_workbase)
    if swr_1 == 0 or swr_2 == 0:
        print('No SWR has been found in one of the swr folders.')
        return
    # get formateed x labels
    add_list_set = set(swr_1.keys()) | set(swr_2.keys())
    add_list = list(add_list_set)
    add_list.sort()
    # get formatted swr data
    swr_clctor = {swr_1_legend: [], swr_2_legend: []}
    for an_add in add_list:
        prev_swr_1 = swr_clctor[swr_1_legend]
        if not an_add in swr_1.keys():
            prev_swr_1.append(0)
        else:
            prev_swr_1.append(swr_1[an_add])

        prev_swr_2 = swr_clctor[swr_2_legend]
        if not an_add in swr_2.keys():
            prev_swr_2.append(0)
        else:
            prev_swr_2.append(swr_2[an_add])
        swr_clctor.update({swr_1_legend: prev_swr_1, swr_2_legend: prev_swr_2})
    # plot
    x = np.array(add_list)  # the label locations
    width = 0.5  # the width of the bars
    multiplier = 0
    fig, ax = plt.subplots(layout='constrained')
    for a_legend, counts in swr_clctor.items():
        offset = width * multiplier
        a_count = ax.bar(x + offset, counts, width, label=a_legend)
        ax.bar_label(a_count, padding=3)
        multiplier += 1
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Counts')
    ax.set_xlabel('Number of addends')
    ax.set_title('Count of SWR scenarios')
    ax.set_xticks(x + width, add_list)
    ax.legend(loc='best')
    if dump_pic_path:
        plt.savefig(dump_pic_path, dpi=400)
        plt.close()
    else:
        plt.show()
    return swr_clctor, add_list


def swr_pipeline(dump_root: str, sorted_root_1: str, sorted_root_2: str,
                 label_1: str, label_2: str, is_unique: bool, is_low_e: bool,
                 proj_ring_seq_1: list, proj_ring_seq_2: list, legend_swr_1_to_2: str, legend_swr_2_to_1: str,
                 group_symbol: str, file_mid_name: str = None,
                 show_C_label=False, C_label_color="black", C_label_transparency=1,
                 pentagon_color="orange", pentagon_transparency=0.5,
                 sphere_ratio=0.8, parr_ratio=0.4, addon_color: str = None, addon_label_size: int = 200,
                 show_addon_nums: bool = False, addon_nums_rotation: int = 0,
                 dpi=400, fontsize: int = 10):
    """

    Args:
        # For dump convenience
        dump_root: The root to dump All the SWR results.
        label_1: The label of system 1, spiral code recommended.
        label_2: The label of system 2, spiral code recommended.
        # For find swr
        sorted_root_1: The highly structured information source. Cook disordered result for system 1.
        sorted_root_2: Cook disordered result for system 2.
        is_low_e: set True to apply the energy criterion,
        is_unique: set True to keep one SWR pair for each queried atoms if there are many target atoms met requirements,
        # For plot swr
        proj_ring_seq_1: proj_ring_seq for system 1.
        proj_ring_seq_2: proj_ring_seq for system 2.
        group_symbol: symbol of groups.
        # For count swr
        legend_swr_1_to_2: the legend on count-SWR-picture, swr from system 1 to system 2.
        legend_swr_2_to_1: the legend for swr from system 1 to system 2.

        Parameters above MUST be provide.
        Rest of the parameters are optional to tailor pictures, etc.

    """
    from autosteper.plotter import FullereneDataParser_Plotter

    a_plotter = FullereneDataParser_Plotter()

    find_SWR(q_sorted_root=sorted_root_1, tgt_sorted_root=sorted_root_2, is_unique=is_unique, is_low_e=is_low_e,
             swr_dump_path=fr'{dump_root}/q_{label_1}_to_tgt_{label_2}',
             group_symbol=group_symbol, file_mid_name=file_mid_name
             )

    a_plotter.plot_swr(src_swr_root=fr'{dump_root}/q_{label_1}_to_tgt_{label_2}',
                       proj_ring_seq=proj_ring_seq_1, C_label_transparency=C_label_transparency,
                       show_C_label=show_C_label, C_label_color=C_label_color,
                       pentagon_transparency=pentagon_transparency, pentagon_color=pentagon_color,
                       sphere_ratio=sphere_ratio, parr_ratio=parr_ratio,
                       addon_color=addon_color, addon_label_size=addon_label_size,
                       show_addon_nums=show_addon_nums, addon_nums_rotation=addon_nums_rotation,
                       fontsize=fontsize, dpi=dpi, group_symbol=group_symbol
                       )

    find_SWR(q_sorted_root=sorted_root_2, tgt_sorted_root=sorted_root_1, is_unique=is_unique, is_low_e=is_low_e,
             swr_dump_path=fr'{dump_root}/q_{label_2}_to_tgt_{label_1}',
             group_symbol=group_symbol, file_mid_name=file_mid_name
             )

    a_plotter.plot_swr(src_swr_root=fr'{dump_root}/q_{label_2}_to_tgt_{label_1}',
                       proj_ring_seq=proj_ring_seq_2, C_label_transparency=C_label_transparency,
                       show_C_label=show_C_label, C_label_color=C_label_color,
                       pentagon_transparency=pentagon_transparency, pentagon_color=pentagon_color,
                       sphere_ratio=sphere_ratio, parr_ratio=parr_ratio,
                       addon_color=addon_color, addon_label_size=addon_label_size,
                       show_addon_nums=show_addon_nums, addon_nums_rotation=addon_nums_rotation,
                       fontsize=fontsize, dpi=dpi, group_symbol=group_symbol
                       )

    count_SWR(swr_1_workbase=fr'{dump_root}/q_{label_1}_to_tgt_{label_2}',
              swr_2_workbase=fr'{dump_root}/q_{label_2}_to_tgt_{label_1}',
              swr_1_legend=legend_swr_1_to_2,
              swr_2_legend=legend_swr_2_to_1,
              dump_pic_path=rf'{dump_root}/swr_count_{label_2}_and_{label_1}.png')


def clc_status_unit(info_path: str):
    status_codes = ['1', '2', '3', '4', '5', '6', '7']
    status_counter = {}
    for a_code in status_codes:
        status_counter.update({a_code: 0})
    with open(info_path, 'r') as f_r:
        for a_line in f_r.readlines():
            a_code = a_line.split()[0]
            status_counter[a_code] = status_counter[a_code] + 1
    return status_counter


def clc_failed(workbase: str, dump_pic_path: str = None, ylim: int = None):
    failed_types = ["radical", "transfer", "3ring", "broken", "bridge", "cluster", "inner"]
    failed_counts = {}
    for a_type in failed_types:
        failed_counts.update({a_type: []})
    add_list = []
    for a_folder in os.listdir(workbase):
        an_add = int(a_folder.split('add')[0])
        add_list.append(an_add)
    add_list.sort()
    for an_add in add_list:
        a_folder = f'{an_add}addons'
        a_failed_path = os.path.join(workbase, a_folder, 'failed_job_paths')
        a_status_counter = clc_status_unit(a_failed_path)
        for a_code, a_count in a_status_counter.items():
            a_failed_type = failed_types[int(a_code) - 1]
            previous_count = failed_counts[a_failed_type]
            previous_count.append(a_count)
            failed_counts.update({a_failed_type: previous_count})
    width = 0.5
    fig, ax = plt.subplots()
    bottom = np.zeros(len(add_list))
    for failed_type, count in failed_counts.items():
        p = ax.bar(add_list, count, width, label=failed_type, bottom=bottom)
        bottom += count
    ax.bar_label(p, padding=3)
    ax.set_title("Count of different failed types")
    ax.legend(loc="upper left")
    ax.set_xticks(add_list)
    ax.set_xticklabels(add_list)
    ax.set_xlabel('Addon number')
    ax.set_ylabel('Count')
    if ylim:
        ax.set_ylim([0, max(ylim, bottom[-1] + 5)])
    if dump_pic_path:
        plt.savefig(dump_pic_path, dpi=400)
        plt.close()
    else:
        plt.show()
    return failed_counts, add_list


class Path_Parser():
    def __init__(self):
        pass

    def get_path_from_sim(self, dump_root: str, pristine_cage_path: str, sim_workbase: str,
                          sim_step: int, sim_start: int, q_add_num: int, q_path_rank: int, q_isomer_rank: int,
                          is_ctr_path: bool = False, ctl_parent_num: int = None, max_path_num: int = None,
                          is_mix: bool = False):
        cwd_path = os.getcwd()
        self.add_num_list = list(range(sim_start, q_add_num + sim_step, sim_step))
        pristine_cage = read(pristine_cage_path)
        pristine_cage_name = os.path.basename(pristine_cage_path)
        cage_size = len(pristine_cage)
        dump_root = os.path.abspath(dump_root)
        os.makedirs(dump_root, exist_ok=True)
        os.chdir(dump_root)
        low_e_dump = os.path.abspath(f'{q_add_num}_addons_low_e_isomers')
        os.makedirs(low_e_dump, exist_ok=True)

        def _control_parent_info(parent_info: pd.DataFrame, info_path: str, pre_parent_info: pd.DataFrame):
            for a_cage in parent_info.keys():
                parent_names = parent_info[a_cage][0]
                name_e_map = {}
                for a_parent in parent_names:
                    name_e_map.update({pre_parent_info[a_parent][1]: a_parent})
                new_parent_names = []
                for a_e in sorted(name_e_map.keys())[:ctl_parent_num]:
                    new_parent_names.append(name_e_map[a_e])
                parent_info[a_cage][0] = new_parent_names
            parent_info.to_pickle(info_path)
            return parent_info

        def _query_isomer_unit():
            def _get_pathway_unit(idx: int, pathway: list, cage_name: str, e_list: list, name_list: list):
                if idx == 0:
                    _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                    final_pathway = [init_pathway, addon_set, *pathway]
                    q_isomer_e = parent_info_list[idx][cage_name][0]
                    final_e_list = [0, q_isomer_e, *e_list]
                    final_name_list = [pristine_cage_name, cage_name, *name_list]
                    e_areas.append(sum(final_e_list))
                    pathways.append(final_pathway)
                    e_lists.append(final_e_list)
                    name_lists.append(final_name_list)
                    return
                else:
                    new_name_list = [cage_name, *name_list]
                    _, addon_set, _0 = name2seq(cage_size=cage_size, name=cage_name)
                    parent_name_list = parent_info_list[idx][cage_name][0]
                    q_isomer_e = parent_info_list[idx][cage_name][1]
                    new_e_list = [q_isomer_e, *e_list]
                    for a_parent in parent_name_list:
                        _, parent_addon_set, _0 = name2seq(cage_size=cage_size, name=a_parent)
                        new_addon_set = addon_set - parent_addon_set
                        new_pathway = [new_addon_set, *pathway]
                        _get_pathway_unit(idx - 1, pathway=new_pathway, cage_name=a_parent, e_list=new_e_list,
                                          name_list=new_name_list)

            q_isomer_e = a_passed_info['energy'][q_rank]
            q_isomer_name = a_passed_info['name'][q_rank]
            q_isomer_atoms = read(a_passed_info['xyz_path'][q_rank])

            _, add_set, _0 = name2seq(name=q_isomer_name, cage_size=cage_size)
            add_for_dump = '_'.join([str(x + 1) for x in sorted(list(add_set))])
            write(filename=os.path.join(low_e_dump, f'{add_for_dump}_rank_{q_rank + 1}.xyz'),
                  images=q_isomer_atoms, format='xyz', comment=str(q_isomer_e))
            init_pathway = {0}
            pathways = []
            e_lists = []
            e_areas = []
            name_lists = []
            # Get path info recursively
            idx = len(self.add_num_list) - 1
            _get_pathway_unit(idx=idx, pathway=[], cage_name=q_isomer_name, e_list=[], name_list=[])
            all_pathways.extend(pathways)
            all_e_lists.extend(e_lists)
            all_e_areas.extend(e_areas)
            all_name_lists.extend(name_lists)

        # get infos
        parent_info_list = []
        if is_ctr_path:
            assert len(self.add_num_list) > 2, print(
                'Current quary addon number is less than 3 steps away from the begining.\n'
                'Do not need control pathes.')
            for ii, i in enumerate(self.add_num_list):
                info_path = os.path.join(sim_workbase, f'{i}addons', 'parent_info.pickle')
                a_parent_info = pd.read_pickle(info_path)
                if ii < 2:
                    parent_info_list.append(a_parent_info)
                else:
                    pre_parent_info = pd.read_pickle(
                        os.path.join(sim_workbase, f'{self.add_num_list[ii - 1]}addons', 'parent_info.pickle'))
                    a_controled_info = _control_parent_info(parent_info=a_parent_info, info_path=info_path,
                                                            pre_parent_info=pre_parent_info)
                    parent_info_list.append(a_controled_info)
        else:
            for i in self.add_num_list:
                info_path = os.path.join(sim_workbase, f'{i}addons', 'parent_info.pickle')
                a_parent_info = pd.read_pickle(info_path)
                parent_info_list.append(a_parent_info)

        passed_info_maps = []
        passed_info_list = []
        for i in self.add_num_list:
            a_passed_info = pd.read_pickle(os.path.join(sim_workbase, f'{i}addons', 'passed_info.pickle'))
            passed_info_list.append(a_passed_info)
            a_info_map = dict(zip(a_passed_info['name'], a_passed_info['xyz_path']))
            passed_info_maps.append(a_info_map)

        # Query the path of low e isomers
        Max_q_rank = len(a_passed_info['energy'])
        all_pathways = []
        all_e_lists = []
        all_e_areas = []
        all_name_lists = []
        q_isomer_rank = min(q_isomer_rank, Max_q_rank)
        if is_mix:
            for q_rank in range(q_isomer_rank):
                _query_isomer_unit()
        else:
            q_rank = q_isomer_rank - 1
            _query_isomer_unit()

        info = pd.DataFrame(
            {'pathway': all_pathways, 'name': all_name_lists, 'e_list': all_e_lists, 'e_area': all_e_areas})
        sorted_info = info.sort_values(by='e_area')
        sorted_info.index = sorted(sorted_info.index)
        if is_ctr_path:
            sorted_info = sorted_info[sorted_info.index < max_path_num]

        sorted_info.to_pickle(f'raw_pathway_info.pickle')
        sorted_info[:1000].to_excel(f'raw_pathway_info.xlsx')

        # Get relative energy info
        e_lists = sorted_info['e_list']
        e_array = np.array(e_lists[0])
        for i in e_lists[1:]:
            e_array = np.vstack((e_array, np.array(i)))

        if len(e_array.shape) == 1:
            raise RuntimeError(f'There are only one root. Please check parameters.')

        rel_e_array = np.ones_like(e_array)
        for i in range(len(self.add_num_list) + 1):
            rel_e_array[:, i] = e_array[:, i] - min(e_array[:, i])

        pathway_xyz_path_dict = {}
        # Dump path related xyz file
        for path_rank in range(min(q_path_rank, len(e_lists))):
            xyz_path_list = []
            name_list = sorted_info['name'][path_rank]
            one_path_workbase = os.path.join(dump_root, f'path_rank_{path_rank + 1}')
            os.makedirs(exist_ok=True, name=one_path_workbase)
            os.chdir(one_path_workbase)
            log_path = f'path_rank_{path_rank + 1}_traj.log'
            write(filename=log_path, images=pristine_cage, format='xyz', append=True)
            for idx, add_num in enumerate(self.add_num_list):
                a_name_path_map = passed_info_maps[idx]
                name = name_list[idx + 1]
                _, add_set, _0 = name2seq(name=name, cage_size=cage_size)
                add_for_dump = '_'.join([str(x + 1) for x in sorted(list(add_set))])
                new_filename = f'{add_num}_addons_{add_for_dump}.xyz'
                shutil.copy(src=a_name_path_map[name], dst=new_filename)
                xyz_path_list.append(os.path.abspath(new_filename))
                atoms = read(new_filename, format='xyz')
                write(filename=log_path, images=atoms, format='xyz', append=True)
            shutil.copy(src=pristine_cage_path, dst='pristine_cage.xyz')
            pathway_xyz_path_dict.update({path_rank + 1: xyz_path_list})
        pathway_xyz_path_dict.update({'pristine_cage_path': os.path.abspath(r'pristine_cage.xyz')})
        os.chdir(dump_root)
        pathway_xyz_path_df = pd.DataFrame(pathway_xyz_path_dict)
        pathway_xyz_path_df.to_pickle('pathway_xyz_path_map.pickle')  # For refine
        # Plot relative energy for path rank
        rel_e_array = rel_e_array[:q_path_rank, :]
        rel_e_df = pd.DataFrame(rel_e_array)
        rel_e_df.columns = [0, *self.add_num_list]
        rel_e_df.index = np.array(rel_e_df.index) + 1
        # rel_e_df.to_pickle('test.pickle')
        plot_pathway_heatmap(rel_e_df)
        os.chdir(cwd_path)

    def update_path_info(self, old_workbase: str, new_workbase: str,
                         refine_top_k: int, refine_mode: str, refiner_para: dict):
        refiner = switch_optimizers(mode=refine_mode, para=refiner_para)
        if refiner.has_parity == True:
            refiner.group = refiner_para['group']
        cwd_update = os.getcwd()
        old_xyz_map = pd.read_pickle(os.path.join(old_workbase, 'pathway_xyz_path_map.pickle'))
        old_rel_e_info = pd.read_pickle(os.path.join(old_workbase, 'pathway_rel_e.pickle'))

        new_workbase = os.path.abspath(new_workbase)
        refine_workbase = os.path.join(new_workbase, 'refine_workbase')
        odd_workbase = os.path.join(refine_workbase, 'odd')
        odd_raw = os.path.join(odd_workbase, 'raw')
        os.makedirs(odd_raw, exist_ok=True)
        even_workbase = os.path.join(refine_workbase, 'even')
        even_raw = os.path.join(even_workbase, 'raw')
        os.makedirs(even_raw, exist_ok=True)

        for an_old_path_rank in range(refine_top_k):
            for an_old_xyz_path in old_xyz_map[an_old_path_rank + 1]:
                old_filename = os.path.basename(an_old_xyz_path)
                add_num = int(old_filename.split('_')[0])
                if add_num % 2 == 1:
                    shutil.copy(src=an_old_xyz_path, dst=os.path.join(odd_raw, old_filename))
                else:
                    shutil.copy(src=an_old_xyz_path, dst=os.path.join(even_raw, old_filename))
        new_name_xyz_path_map = {}
        new_name_e_map = {}
        if len(os.listdir(odd_raw)) > 0:
            os.chdir(odd_workbase)
            refiner.set_folders()
            refiner.opt()
            passed_info = pd.read_pickle(r'passed_info.pickle')
            for idx, a_name in enumerate(passed_info['name']):
                new_name_e_map.update({a_name + '.xyz': passed_info['energy'][idx]})
                new_name_xyz_path_map.update({a_name + '.xyz': passed_info['xyz_path'][idx]})
        if len(os.listdir(even_raw)) > 0:
            os.chdir(even_workbase)
            refiner.set_folders()
            refiner.opt()
            passed_info = pd.read_pickle(r'passed_info.pickle')
            for idx, a_name in enumerate(passed_info['name']):
                new_name_e_map.update({a_name + '.xyz': passed_info['energy'][idx]})
                new_name_xyz_path_map.update({a_name + '.xyz': passed_info['xyz_path'][idx]})
        for an_old_path_rank in range(refine_top_k):
            for an_old_xyz_path in old_xyz_map[an_old_path_rank + 1]:
                old_filename = os.path.basename(an_old_xyz_path)
                add_num = int(old_filename.split('_')[0])
                old_rel_e_info[add_num][an_old_path_rank + 1] = new_name_e_map[old_filename]
        old_rel_e_info = old_rel_e_info[old_rel_e_info.index < refine_top_k + 1]
        old_rel_e_info['pathway'] = old_rel_e_info.iloc[:, 0:len(old_rel_e_info.columns) - 1].sum(axis=1)
        for a_column in old_rel_e_info.columns:
            old_rel_e_info[a_column] = np.array(old_rel_e_info[a_column]) - min(old_rel_e_info[a_column])
        old_rel_e_info = old_rel_e_info.sort_values(by='pathway')
        old_ranks = old_rel_e_info.index.copy()
        old_rel_e_info.index = sorted(old_rel_e_info.index)

        pathway_xyz_path_dict = {}
        pristine_cage_path = old_xyz_map['pristine_cage_path'][0]
        pristine_cage = read(pristine_cage_path)
        for idx, an_old_rank in enumerate(old_ranks):
            if idx == refine_top_k:
                break
            xyz_path_list = []
            os.chdir(new_workbase)
            os.makedirs(f'path_rank_{idx + 1}', exist_ok=True)
            os.chdir(f'path_rank_{idx + 1}')
            write(filename=r'pristine_cage.xyz', images=pristine_cage)
            write(filename=f'path_rank_{idx + 1}_traj.log', images=pristine_cage, format='xyz')
            for an_old_xyz_path in old_xyz_map[an_old_rank]:
                old_filename = os.path.basename(an_old_xyz_path)
                new_xyz_path = new_name_xyz_path_map[old_filename]
                new_atoms = read(new_xyz_path)
                write(filename=old_filename, images=new_atoms)
                write(filename=f'path_rank_{idx + 1}_traj.log', images=new_atoms, append=True, format='xyz')
                xyz_path_list.append(os.path.abspath(old_filename))
            pathway_xyz_path_dict.update({idx + 1: xyz_path_list})
        pathway_xyz_path_dict.update({'pristine_cage_path': os.path.abspath(r'pristine_cage.xyz')})
        os.chdir(new_workbase)
        pathway_xyz_path_df = pd.DataFrame(pathway_xyz_path_dict)
        pathway_xyz_path_df.to_pickle('pathway_xyz_path_map.pickle')  # For refine
        old_rel_e_info = old_rel_e_info.drop('pathway', axis=1)
        plot_pathway_heatmap(old_rel_e_info)
        os.chdir(cwd_update)

    def re_generate_pathways(self, old_workbase: str, new_workbase: str,
                             step: int, last_log_mode: str, keep_top_k_pathway: int,
                             group: str = None, file_mid_name: str = None, dpi: int = 400):
        cwd_re_gen = os.getcwd()

        new_workbase = os.path.abspath(new_workbase)
        os.makedirs(new_workbase, exist_ok=True)
        os.chdir(new_workbase)
        os.makedirs('old_logs', exist_ok=True)
        abs_old_logs = os.path.abspath('old_logs')
        if last_log_mode == 'xtb':
            chk_name = 'xtbopt.log'
        elif last_log_mode == 'gauss':
            chk_name = 'gau.log'
        old_even_cooking = os.path.join(old_workbase, 'refine_workbase', 'even', 'cooking')
        for a_folder in os.listdir(old_even_cooking):
            os.chdir(old_even_cooking)
            os.chdir(a_folder)
            shutil.copy(src=chk_name, dst=os.path.join(abs_old_logs, f'{a_folder}.log'))
        old_odd_cooking = os.path.join(old_workbase, 'refine_workbase', 'odd', 'cooking')
        for a_folder in os.listdir(old_odd_cooking):
            os.chdir(old_odd_cooking)
            os.chdir(a_folder)
            shutil.copy(src=chk_name, dst=os.path.join(abs_old_logs, f'{a_folder}.log'))
        os.chdir(new_workbase)
        cook_disordered(disordered_root=abs_old_logs, dump_root=new_workbase, group=group, file_mid_name=file_mid_name,
                        step=step, log_mode=last_log_mode, keep_top_k_pathway=keep_top_k_pathway, dpi=dpi)
