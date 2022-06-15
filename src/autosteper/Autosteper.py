import shutil
import os
import pandas as pd
from ase.units import Hartree, eV
from autosteper.cage import name2seq, blk_list
from autosteper.generator import Generator
from autosteper.optimizers import *
from autosteper.path_parser import Path_Parser


Hartree2eV = Hartree/eV
class AutoSteper():
    def __init__(self, para: dict):
        self.cage = Cage(pristine_path=para['pristine_cage'], workbase=para['workbase'])
        self.generator = Generator(para['gen_para'])
        self.checker = Checker(chk_skin=self.generator.skin, group=self.generator.group, cage=self.cage)
        if para['opt_mode'] == 'xtb':
            self.optimizer = XTB_Optimizer(para['opt_para'], checker=self.checker, cage= self.cage)
        elif para['opt_mode'] == 'ase':
            self.optimizer = ASE_Optimizer(para['opt_para'], checker=self.checker, cage= self.cage)
        elif para['opt_mode'] == 'gaussian':
            self.optimizer = Gaussian_Optimizer(para['opt_para'], checker=self.checker, cage= self.cage)
        if 'run_para' in para.keys():
            self.start = para['run_para']['start']
            self.stop = para['run_para']['stop']
            self.step = para['run_para']['step']
            self.wht_list_para = para['run_para']['wht_list_para']

        if 'path_para' in para.keys():
            self.path_parser = Path_Parser(path_para=para['path_para'], step=self.step,
                                           workbase=para['workbase'], start=self.start,
                                           q_cage=self.cage, optimizer=self.optimizer)

        if 'random_para' in para.keys():
            self.addon_list = para['random_para']['addon_list']
            self.random_num = para['random_para']['random_num']
            self.try_times = para['random_para']['try_times']

        if 'pre_scan_para' in para.keys():
            self.is_pre_scan = True
            self.start_ps_num = para['pre_scan_para']['start_ps_num']
            self.ps_calc = para['pre_scan_para']['calculator']
            self.ps_cut_para = para['pre_scan_para']['ps_cut_para']
        else:
            self.is_pre_scan = False

        if 'blk_para' in para.keys():
            self.cage.blk_list = blk_list(size=self.cage.size, blk_para=para['blk_para'])

    def _get_yes_info(self, all_parent_info: dict=None):
        cwd_ = os.getcwd()
        deep_yes_path = os.path.join(cwd_, 'deep_yes_info.pickle')
        if os.path.exists(deep_yes_path):
            return deep_yes_path, os.path.abspath('flat_yes_info.pickle')
        flat_yes_info = {}
        name_list = []
        energy_list = []
        xyz_path_list = []
        for yes_paths_name in ['init_yes_paths', 'yes_paths']:
            if not os.path.exists(yes_paths_name):
                continue
            with open(yes_paths_name, 'r') as f:
                for a_log in f.readlines():
                    a_log = a_log.strip()
                    a_name = os.path.basename(os.path.split(a_log)[0])
                    name_list.append(a_name)
                    a_xyz_path = os.path.splitext(a_log)[0] + '.xyz'
                    xyz_path_list.append(a_xyz_path)
                    if self.optimizer.mode in ['xtb', 'gaussian']:
                        with open(a_xyz_path, 'r') as xyz_f:
                            xyz_f.readline()
                            energy_line = xyz_f.readline()
                            energy = float(energy_line.split()[1])
                    elif self.optimizer.mode == 'ase':
                        with open(a_log, 'r') as xyz_f:
                            energy_line = xyz_f.readlines()[-1]
                            energy = float(energy_line.split()[-2].split('*')[0])
                    energy_list.append(energy)
                    if all_parent_info == None:
                        flat_yes_info.update({a_name: [energy]})
                    else:
                        flat_yes_info.update({a_name: [all_parent_info[a_name][0], energy]})

        deep_yes_info = pd.DataFrame({'name': name_list, 'energy': energy_list, 'xyz_path': xyz_path_list})
        sorted_deep_yes = deep_yes_info.sort_values(by='energy')
        sorted_deep_yes.index = sorted(sorted_deep_yes.index)
        sorted_deep_yes.to_pickle(path=deep_yes_path)
        flat_yes_info_df = pd.DataFrame(flat_yes_info)
        flat_yes_info_df.to_pickle(path='flat_yes_info.pickle')
        return deep_yes_path, os.path.abspath('flat_yes_info.pickle')

    def _post_pre_scan(self):
        def _copy_file():
            name = self.pre_scan_map[a_e]
            os.symlink(src=os.path.abspath(os.path.join(f'{self.optimizer.path_raw_init}', f'{name}.xyz')),
                       dst=os.path.abspath(os.path.join('post_pre_scan_raw', f'{name}.xyz')))

        os.makedirs('post_pre_scan_raw', exist_ok=True)
        e_list = sorted(self.pre_scan_map.keys())
        if self.ps_cut_para['mode'] == 'rank':
            for a_rank, a_e in enumerate(e_list):
                if a_rank == self.ps_cut_para['rank']:
                    break
                _copy_file()
        elif self.ps_cut_para['mode'] == 'value':
            for a_e in e_list:
                if a_e >= e_list[0] + self.ps_cut_para['value']/Hartree2eV:
                    break
                _copy_file()
        elif self.ps_cut_para['mode'] == 'value_rank':
            for a_rank, a_e in enumerate(e_list):
                if a_e >= e_list[0] + self.ps_cut_para['value']/Hartree2eV:
                    break
                if a_rank == self.ps_cut_para['rank']:
                    break
                _copy_file()
        elif self.ps_cut_para['mode'] == None:
            for a_e in e_list:
                _copy_file()
        else:
            raise RuntimeError('Please check your pre-scan cutoff mode keyword.\nCurrently only support: None, rank, value and value_rank.')
        self.optimizer.path_raw_init = 'post_pre_scan_raw'


    def random(self):
        cwd_ = str(os.getcwd())
        for an_add_num in self.addon_list:
            if an_add_num < 5:
                self.cage.set_add_num(an_add_num)
                gen_out_path = f'{self.cage.name}_{self.cage.add_num}_addons.out'
                self.generator.gen_seq(cage=self.cage,
                                       gen_out_path=gen_out_path,
                                       mode='base')
                with open(gen_out_path, 'r') as file:
                    gen_num = len(file.readlines())
                if gen_num > self.random_num:
                    os.remove(gen_out_path)
                    os.chdir(cwd_)
                else:
                    print(f'The total number of {self.cage.add_num} addon isomers for '
                          f'cage in {self.cage.pristine_path} is less than random opt number, '
                          f'random opt procedure is terminated for this system.')
                    os.chdir(cwd_)
                    shutil.rmtree(self.cage.addon_path)
                    continue

            opt_times = 0
            while True:
                self.cage.set_add_num(an_add_num)
                gen_out_path = f'{self.cage.name}_{self.cage.add_num}_addons.out'
                self.generator.gen_seq(cage=self.cage,
                                       gen_out_path=gen_out_path,
                                       mode='random',
                                       random_num=self.random_num)
                self.optimizer.set_init_folders()
                self.generator.build(is_first=True,
                                     gen_out_path=gen_out_path,
                                     dump_folder=self.optimizer.path_raw_init,
                                     cage=self.cage,
                                     prev_xyz_path=self.cage.pristine_path)
                random_status = self.optimizer.opt()
                if random_status == 2:
                    os.chdir(cwd_)
                    shutil.rmtree(self.cage.addon_path)
                else:
                    _, _ = self._get_yes_info()
                    break
                if opt_times == self.try_times - 1:
                    print(f'Random opt procedure has performed {self.try_times} times, '
                          f'still get wrong status, it\'s been terminated, please check this system.')
                    break
                else:
                    opt_times = opt_times + 1


    def _first_step(self):
        self.cage.set_add_num(self.start)
        if self.is_pre_scan and self.cage.add_num >= self.start_ps_num:
            is_pre_scan = True
            self.pre_scan_map = dict()
        else:
            is_pre_scan = False
        gen_out_path = f'{self.cage.name}_{self.cage.add_num}_addons.out'
        self.generator.gen_seq(cage=self.cage,
                               gen_out_path=gen_out_path,
                               mode='base')
        self.optimizer.set_init_folders()
        if is_pre_scan:
            self.pre_scan_map = self.generator.build(is_first=True,
                                 gen_out_path=gen_out_path,
                                 dump_folder=self.optimizer.path_raw_init,
                                 cage=self.cage,
                                 prev_xyz_path=self.cage.pristine_path,
                                 calc=self.ps_calc,
                                 pre_scan_map=self.pre_scan_map
                                 )
            self._post_pre_scan()
        else:
            _ = self.generator.build(is_first=True,
                                 gen_out_path=gen_out_path,
                                 dump_folder=self.optimizer.path_raw_init,
                                 cage=self.cage,
                                 prev_xyz_path=self.cage.pristine_path)
        step_status = self.optimizer.opt()
        if step_status == 0:
            return 0
        else:
            self.prev_deep_yes, self.prev_flat_yes = self._get_yes_info()
        return step_status

    def _take_a_step(self):
        def _build_unit():
            prev_xyz_path = prev_deep_yes['xyz_path'][idx]
            prev_name = prev_deep_yes['name'][idx]
            prev_seq, prev_addon_set, _ = name2seq(name=prev_name, cage_size=self.cage.size)
            sub_nauty_o = os.path.join(sub_nauty_path, f'{prev_name}.out')
            self.generator.gen_seq(mode='step', cage=self.cage, gen_out_path=sub_nauty_o, prev_seq=prev_seq)
            if is_pre_scan:
                self.all_parent_info, self.pre_scan_map = self.generator.build(is_first=False,
                                                                               gen_out_path=sub_nauty_o,
                                                                               dump_folder=self.optimizer.path_raw_init,
                                                                               cage=self.cage,
                                                                               prev_xyz_path=prev_xyz_path,
                                                                               parent_info=self.all_parent_info,
                                                                               prev_addon_set=prev_addon_set,
                                                                               parent_name=prev_name,
                                                                               calc=self.ps_calc,
                                                                               pre_scan_map=self.pre_scan_map
                                                                               )
            else:
                self.all_parent_info, _ = self.generator.build(is_first=False,
                                                               gen_out_path=sub_nauty_o,
                                                               dump_folder=self.optimizer.path_raw_init,
                                                               cage=self.cage,
                                                               prev_xyz_path=prev_xyz_path,
                                                               parent_info=self.all_parent_info,
                                                               prev_addon_set=prev_addon_set,
                                                               parent_name=prev_name
                                                               )

        self.cage.set_add_num(self.new_add_num)

        if self.is_pre_scan and self.cage.add_num >= self.start_ps_num:
            is_pre_scan = True
            self.pre_scan_map = dict()
        else:
            is_pre_scan = False

        self.optimizer.set_init_folders()
        self.all_parent_info = {}
        sub_nauty_path = 'sub_nauty'
        os.makedirs(sub_nauty_path, exist_ok=True)

        prev_deep_yes = pd.read_pickle(self.prev_deep_yes)
        if self.wht_list_para['mode'] == 'rank':
            for idx in range(min(len(prev_deep_yes['energy']), self.wht_list_para['rank'])):
                _build_unit()
        elif self.wht_list_para['mode'] == 'value_rank':
            for idx, a_prev_energy in enumerate(prev_deep_yes['energy']):
                if a_prev_energy > prev_deep_yes['energy'][0] + self.wht_list_para['value']/Hartree2eV:
                    idx = idx - 1
                    break
                if idx == self.wht_list_para['rank']:
                    idx = idx - 1
                    break
                _build_unit()
        elif self.wht_list_para['mode'] == 'value':
            for idx, a_prev_energy in enumerate(prev_deep_yes['energy']):
                if a_prev_energy > prev_deep_yes['energy'][0] + self.wht_list_para['value']/Hartree2eV:
                    idx = idx - 1
                    break
                _build_unit()
        elif self.wht_list_para['mode'] == None:
            for idx, a_prev_energy in enumerate(prev_deep_yes['energy']):
                _build_unit()
        else:
            raise RuntimeError('Please check your run cutoff mode keyword.\nCurrently only support: None, rank, value and value_rank.')
        prev_flat_yes = pd.read_pickle(self.prev_flat_yes)
        wht_list_names = prev_deep_yes['name'][: idx+1]
        new_prev_flat = prev_flat_yes[wht_list_names]
        new_prev_flat.to_pickle(self.prev_flat_yes)


        all_parent_info_df = pd.DataFrame(self.all_parent_info)
        all_parent_info_df.to_pickle('all_parent_info.pickle')

        if is_pre_scan:
            self._post_pre_scan()

        if self.cage.blk_list:
            if (self.cage.blk_list.end_chk_num - self.step) > self.cage.add_num >= self.cage.blk_list.start_clct_num:
                self.cage.has_blk_list = True
            else:
                self.cage.has_blk_list = False

        step_status = self.optimizer.opt()
        self.prev_deep_yes, self.prev_flat_yes = self._get_yes_info(all_parent_info=self.all_parent_info)

        if self.cage.has_blk_list:
            if self.cage.blk_list.clct_unstb:
                self.cage.blk_list.clct_unstable(info_path=self.prev_deep_yes)
            else:
                self.cage.blk_list.blk_list_arr = self.cage.blk_list.failed_arr

        return step_status

    def run(self):
        step_status = self._first_step()
        if step_status == 0:
            print("AutoSteper failed on first step.")
        else:
            while step_status == 1:
                self.new_add_num = self.cage.add_num + self.step
                if self.stop != None:
                    if self.new_add_num > self.stop:
                        print("Normal Termination of AutoSteper.")
                        break
                step_status = self._take_a_step()
                if step_status == 0:
                    print(f"AutoSteper failed after optimizing {self.cage.add_num} addon system.")
                    break




