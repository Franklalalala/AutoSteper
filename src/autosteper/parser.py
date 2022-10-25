import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from ase import Atoms
from autosteper.cage import Cage, seq2name
from autosteper.optimizers import *
from autosteper.tools import get_low_e_ranks, strip_extraFullerene, get_G, simple_dump_pathway, get_pathway_info, \
    get_connection, simple_parse_logs, map_SWR, lazy_file_mid_name
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
        pristine_cage, q_add_set = strip_extraFullerene(atoms=q_atoms, group=group, cage_size=cage_size)
    elif q_seq:
        _, q_add_set, _0 = seq2name(q_seq, q_cage)
    else:
        print('Please input query addon set.\n'
              'Currently only surpport ase.atoms input, sequence input and addon set input.')

    mid_add_num = len(q_add_set)
    if q_cage:
        cage_size = q_cage.size
    else:
        print('Please input cage size.')
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
                             cutoff: dict = None):
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

    q_add_set, cage_size, mid_add_num = _prep_relatives(q_atoms=q_atoms, q_seq=q_seq, q_cage=q_cage)
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
                    group: str = None, cage_size: int = None, file_mid_name: str = None):
    """
    Pipeline to cook disordered log files.
    Sort log files first, get connection and pathway infos, dump pathways.
    """
    disordered_root = os.path.abspath(disordered_root)
    dump_root = os.path.abspath(dump_root)
    sorted_root = os.path.join(dump_root, 'sorted')
    simple_parse_logs(dump_root=sorted_root, src_root=disordered_root,
                      mode=log_mode, group=group, cage_size=cage_size, file_mid_name=file_mid_name)
    get_connection(xyz_root=os.path.join(sorted_root, 'xyz'),
                   connection_dump=os.path.join(sorted_root, 'connection'),
                   step=step)
    get_pathway_info(e_info_path=os.path.join(sorted_root, 'info', 'info.pickle'),
                     xyz_root=os.path.join(sorted_root, 'xyz'),
                     cnt_root=os.path.join(sorted_root, 'connection'),
                     dump_info_path=os.path.join(dump_root, r'pathway_info.pickle'))
    simple_dump_pathway(pathway_info_path=os.path.join(dump_root, r'pathway_info.pickle'),
                        dump_root=os.path.join(dump_root, r'pathways'),
                        src_xyz_root=os.path.join(sorted_root, 'xyz'),
                        top_K=keep_top_k_pathway)


def find_SWR(q_sorted_root: str, tgt_sorted_root: str, swr_dump_path: str,
             step: int, is_low_e: bool = True, is_unique: bool = True,
             group: str = None, cage_size: int = None, file_mid_name: str = None):
    '''
    Find SWR from atoms in q_logs to atoms in tgt_logs.
    Pristine cage needs to be prepared in log root.
    if is_unique is true, for every atoms in q_root, only one SWR target is outputed,
    typically for the lowest energy isomer, here we take the rank info in the name as criteria.
    if is_low_e is true, for every atoms in q_root, only one SWR target is outputed,
    and it should have a lower energy than the 'ought to be' parents.
    '''
    rel_file_mid_name = lazy_file_mid_name(file_mid_name)
    swr_dump_path = os.path.abspath(swr_dump_path)
    os.makedirs(swr_dump_path, exist_ok=True)

    q_xyz_root = os.path.join(q_sorted_root, 'xyz')
    q_atoms = read(os.path.join(q_xyz_root, f'0_addons{rel_file_mid_name}1.xyz'))
    q_G = get_G(q_atoms)
    q_cnt_root = os.path.join(q_sorted_root, 'connection')

    tgt_xyz_root = os.path.join(tgt_sorted_root, 'xyz')
    tgt_atoms = read(os.path.join(tgt_xyz_root, f'0_addons{rel_file_mid_name}1.xyz'))
    tgt_G = get_G(tgt_atoms)
    tgt_cnt_root = os.path.join(tgt_sorted_root, 'connection')

    swr = map_SWR(q_atoms, tgt_atoms)
    # swr = [[22,24], [22,24]]

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
    chked_q_add = {}
    add_num_list = []
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
            _, q_addon_set = strip_extraFullerene(group=group, cage_size=cage_size, atoms=a_q_atoms)
            ept_nodes = q_swr_adj_nodes_set - q_addon_set

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
                tgt_atoms_list = []
                rank_list = []
                for a_tgt_file in os.listdir(tgt_xyz_root):
                    tgt_add_num = int(a_tgt_file.split('_')[0])
                    if tgt_add_num == len(q_addon_set) + step:
                        a_tgt_atoms = read(os.path.join(tgt_xyz_root, a_tgt_file))
                        _, tgt_addon_set = strip_extraFullerene(group=group, cage_size=cage_size, atoms=a_tgt_atoms)
                        if len(done_nodes_for_chk - tgt_addon_set) == 0 and \
                                len(ept_nodes_for_chk - tgt_addon_set) < len(ept_nodes_for_chk):
                            a_tgt_G_for_chk = get_G(a_tgt_atoms)
                            a_GM = isomorphism.GraphMatcher(a_tgt_G_for_chk, a_q_G_for_chk)
                            if a_GM.subgraph_is_isomorphic():
                                tgt_atoms_list.append(a_tgt_atoms)
                                rank_list.append(int(a_tgt_file.split('.')[0].split('_')[-1]) - 1)
                if len(rank_list) > 0:
                    if is_low_e:
                        tgt_info = pd.read_pickle(os.path.join(tgt_sorted_root, 'info', 'info.pickle'))
                        q_info = pd.read_pickle(os.path.join(q_sorted_root, 'info', 'info.pickle'))
                        tgt_e = tgt_info[len(q_addon_set) + step][sorted(rank_list)[0]]
                        cnt = np.load(os.path.join(q_cnt_root, a_file[:-4] + '.npy'))
                        new_q_info = []
                        for idx, a_q_e in enumerate(q_info[len(q_addon_set) + step].dropna()):
                            if not a_q_e == None:
                                if cnt[idx] == 1:
                                    new_q_info.append(a_q_e)
                        if len(new_q_info) != 0:
                            if tgt_e > min(new_q_info):
                                continue
                    ###################################################
                    # dump
                    rank_list = sorted(rank_list)
                    os.chdir(swr_dump_path)
                    os.makedirs(f'{len(q_addon_set)}_to_{len(tgt_addon_set)}_swr_{chked_q_add[len(q_addon_set)] + 1}',
                                exist_ok=True)
                    os.chdir(f'{len(q_addon_set)}_to_{len(tgt_addon_set)}_swr_{chked_q_add[len(q_addon_set)] + 1}')
                    chked_q_add[len(q_addon_set)] += 1
                    for idx, a_rank in enumerate(rank_list):
                        if is_unique and idx > 0:
                            break
                        rank_atoms_map = dict(zip(rank_list, tgt_atoms_list))
                        write(filename='q_atoms.xyz', images=a_q_atoms, format='xyz')
                        write(filename=f'tgt_atoms_rank_{idx + 1}.xyz', images=rank_atoms_map[a_rank], format='xyz')


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


class Path_Parser():
    def __init__(self, path_para: dict, dump_root: str,
                 q_cage: Cage, refiner_para: dict = None):
        self.q_cage = q_cage
        self.workbase = os.path.join(dump_root, f'path_info_{self.q_cage.name}')
        os.makedirs(self.workbase, exist_ok=True)
        #################################################################################
        # Receive key path parameter
        if 'is_mix' in path_para.keys():
            self.is_mix = path_para['is_mix']
        else:
            self.is_mix = False

        self.step = path_para['step']
        self.start = path_para['start']
        self.q_add_num = path_para['q_add_num']
        self.add_num_list = list(range(self.start, self.q_add_num + self.step, self.step))

        self.q_path_rank = path_para['q_path_rank']
        self.q_isomer_rank = path_para['q_isomer_rank']
        self.log_low_e_num = path_para['log_low_e_num']

        if 'ctl_path_para' in path_para.keys():
            self.ctr_path = True
            self.ctl_parent_num = path_para['ctl_path_para']['ctl_parent_num']
            self.max_path_num = path_para['ctl_path_para']['max_path_num']
        else:
            self.ctr_path = False
        #################################################################################
        if refiner_para:
            self.is_refine = True
            self.refine_top_k = refiner_para['refine_top_k']
            self.refine_mode = refiner_para['opt_mode']
            self.refiner = switch_optimizers(mode=self.refine_mode, para=refiner_para)
        else:
            self.is_refine = False

    def get_path_info(self):
        os.chdir(self.workbase)

        def _control_parent_info(parent_info: pd.DataFrame, info_path: str, pre_parent_info: pd.DataFrame):
            for a_cage in parent_info.keys():
                parent_names = parent_info[a_cage][0]
                name_e_map = {}
                for a_parent in parent_names:
                    name_e_map.update({pre_parent_info[a_parent][1]: a_parent})
                new_parent_names = []
                for a_e in sorted(name_e_map.keys())[:self.ctl_parent_num]:
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
                    final_e_list = [self.base_e, q_isomer_e, *e_list]
                    final_name_list = [self.q_cage.name, cage_name, *name_list]
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

            write(filename=q_isomer_name + f'_rank_{q_rank + 1}.xyz', images=q_isomer_atoms, format='xyz',
                  comment=str(q_isomer_e))
            init_pathway = {0}

            pathways = []
            e_lists = []
            e_areas = []
            name_lists = []

            cage_size = self.q_cage.size
            # Get path info recursively
            idx = len(self.add_num_list) - 1
            _get_pathway_unit(idx=idx, pathway=[], cage_name=q_isomer_name, e_list=[], name_list=[])
            all_pathways.extend(pathways)
            all_e_lists.extend(e_lists)
            all_e_areas.extend(e_areas)
            all_name_lists.extend(name_lists)

        def _refine_e(old_info: pd.DataFrame):
            def _refine_unit(is_even: Union[bool, int, float]):
                self.refiner.set_folders()
                os.chdir('raw')
                for add_num_idx, add_num in enumerate(self.add_num_list):
                    if add_num % 2 == is_even:
                        continue
                    for name_list in old_info['name']:
                        name = name_list[add_num_idx + 1]
                        old_path = passed_info_maps[add_num_idx][name]
                        shutil.copy(src=old_path, dst=f'{name}.xyz')
                os.chdir('./..')
                self.refiner.opt()
                new_name_e_map.update(dict(zip(self.refiner.passed_names, self.refiner.e_list)))

            os.makedirs(name='refine_workbase', exist_ok=True)
            os.chdir('refine_workbase')
            new_name_e_map = {}
            if self.refiner.has_parity:
                os.makedirs(name='odd_workbase', exist_ok=True)
                os.makedirs(name='even_workbase', exist_ok=True)
                os.chdir('even_workbase')
                _refine_unit(is_even=True)
                os.chdir('./..')
                os.chdir('odd_workbase')
                _refine_unit(is_even=False)
                os.chdir('./..')
            else:
                _refine_unit(is_even=0.5)
            os.chdir(self.workbase)
            old_e_list = list(old_info['e_list'])
            for add_num_idx, add_num in enumerate(self.add_num_list):
                addon_path = os.path.join(self.q_cage.workbase, f'{add_num}addons')
                for pathway_idx, name_list in enumerate(old_info['name']):
                    name = name_list[add_num_idx + 1]
                    e = new_name_e_map[name]
                    old_e_list[pathway_idx][add_num_idx + 1] = e
            old_info['e_list'] = old_e_list
            new_e_area = []
            for e_list in old_e_list:
                e_area = sum(e_list)
                new_e_area.append(e_area)
            old_info['e_area'] = new_e_area
            os.chdir(self.workbase)
            refined_info = old_info.sort_values(by='e_area')
            refined_info.index = sorted(refined_info.index)
            return refined_info

        self.base_e = 0

        # get infos
        parent_info_list = []
        if self.ctr_path:
            assert len(self.add_num_list) > 2, print(
                'Current quary addon number is less than 3 steps away from the begining.\n'
                'Do not need control pathes.')
            for ii, i in enumerate(self.add_num_list):
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'parent_info.pickle')
                a_parent_info = pd.read_pickle(info_path)
                if ii < 2:
                    parent_info_list.append(a_parent_info)
                else:
                    pre_parent_info = pd.read_pickle(
                        os.path.join(self.q_cage.workbase, f'{self.add_num_list[ii - 1]}addons', 'parent_info.pickle'))
                    a_controled_info = _control_parent_info(parent_info=a_parent_info, info_path=info_path,
                                                            pre_parent_info=pre_parent_info)
                    parent_info_list.append(a_controled_info)
        else:
            for i in self.add_num_list:
                info_path = os.path.join(self.q_cage.workbase, f'{i}addons', 'parent_info.pickle')
                a_parent_info = pd.read_pickle(info_path)
                parent_info_list.append(a_parent_info)

        passed_info_maps = []
        passed_info_list = []
        for i in self.add_num_list:
            a_passed_info = pd.read_pickle(os.path.join(self.q_cage.workbase, f'{i}addons', 'passed_info.pickle'))
            passed_info_list.append(a_passed_info)
            a_info_map = dict(zip(a_passed_info['name'], a_passed_info['xyz_path']))
            passed_info_maps.append(a_info_map)

        # Dump low e isomers to a log
        Max_q_rank = len(a_passed_info['energy'])
        self.log_low_e_num = min(self.log_low_e_num, Max_q_rank)
        log_path = f'top_{self.log_low_e_num}_isomers.log'
        for rank in range(self.log_low_e_num):
            src_path = a_passed_info['xyz_path'][rank]
            if rank == 0:
                shutil.copy(src=src_path, dst=log_path)
            else:
                with open(log_path, 'a') as w_file, open(src_path, 'r') as r_file:
                    for a_line in r_file.readlines():
                        w_file.write(a_line)

        # Query the path of low e isomers
        all_pathways = []
        all_e_lists = []
        all_e_areas = []
        all_name_lists = []
        self.q_isomer_rank = min(self.q_isomer_rank, Max_q_rank)
        if self.is_mix:
            for q_rank in range(self.q_isomer_rank):
                _query_isomer_unit()
        else:
            q_rank = self.q_isomer_rank - 1
            _query_isomer_unit()

        info = pd.DataFrame(
            {'pathway': all_pathways, 'name': all_name_lists, 'e_list': all_e_lists, 'e_area': all_e_areas})
        sorted_info = info.sort_values(by='e_area')
        sorted_info.index = sorted(sorted_info.index)
        sorted_info.to_pickle(f'all_info.pickle')
        if self.ctr_path:
            sorted_info = sorted_info[sorted_info.index < self.max_path_num]
        if self.is_refine:
            sorted_info = _refine_e(old_info=sorted_info[sorted_info.index < self.refine_top_k])

        sorted_info.to_pickle(f'sorted_info.pickle')

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
        rel_e_array = rel_e_array

        # Dump path related xyz file
        for path_rank in range(min(self.q_path_rank, len(e_lists))):
            name_list = sorted_info['name'][path_rank]
            one_path_workbase = os.path.join(self.workbase, f'path_rank_{path_rank + 1}')
            os.makedirs(exist_ok=True, name=one_path_workbase)
            os.chdir(one_path_workbase)
            log_path = f'path_rank_{path_rank + 1}_traj.log'
            write(filename=log_path, images=self.q_cage.atoms, format='xyz', append=True)
            for idx, add_num in enumerate(self.add_num_list):
                a_name_path_map = passed_info_maps[idx]
                name = name_list[idx + 1]
                new_path = f'{name}_addon_{add_num}.xyz'
                shutil.copy(src=a_name_path_map[name], dst=new_path)
                atoms = read(new_path, format='xyz')
                write(filename=log_path, images=atoms, format='xyz', append=True)
        os.chdir(self.workbase)
        np.save(file=f'Path_relative_energy.npy', arr=rel_e_array)

        # Plot relative energy for path rank
        rel_e_array = rel_e_array[:self.q_path_rank, :]
        fig_1 = plt.figure(dpi=400, figsize=(len(self.add_num_list), min(self.q_path_rank, len(e_lists) + 1, 20)))
        cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
        rel_e_df = pd.DataFrame(rel_e_array)
        rel_e_df.columns = [0, *self.add_num_list]
        rel_e_df.index = np.array(rel_e_df.index) + 1
        rel_pathway = np.array(rel_e_df.sum(axis=1))
        rel_e_df['Pathway'] = rel_pathway - rel_pathway[0]
        sns.heatmap(rel_e_df, annot=True, cmap=cmap, linewidths=.5)
        plt.ylabel('Path rank')
        plt.xlabel('Number of addends')
        plt.title('Path relative energy (eV)')
        plt.savefig(f'Path_relative_energy.png')
        plt.close(fig=fig_1)

        # Plot relative energy for isomers.
        isomer_rel_e = np.zeros_like(e_array)
        for idx, add_num in enumerate(self.add_num_list):
            passed_info = passed_info_list[idx]
            isomer_min = passed_info['energy'][0]
            isomer_rel_e[:, idx + 1] = e_array[:, idx + 1] - isomer_min
        isomer_rel_e = isomer_rel_e
        np.save(file=f'Isomer_relative_energy.npy', arr=isomer_rel_e)
        isomer_rel_e = isomer_rel_e[:self.q_path_rank, :]
        # Simple plot
        fig_2 = plt.figure(dpi=400, figsize=(len(self.add_num_list), min(self.q_path_rank, len(e_lists) + 1, 20)))
        cmap = sns.light_palette((260, 75, 60), input="husl", as_cmap=True)
        isomer_rel_e_df = pd.DataFrame(isomer_rel_e)
        isomer_rel_e_df.columns = [0, *self.add_num_list]
        isomer_rel_e_df.index = np.array(isomer_rel_e_df.index) + 1
        sns.heatmap(isomer_rel_e_df, annot=True, cmap=cmap, linewidths=.5)
        plt.ylabel('Path rank')
        plt.xlabel('Number of addends')
        plt.title('Isomer relative energy (eV)')
        plt.savefig(f'Isomer_relative_energy.png')
        plt.close(fig=fig_2)
