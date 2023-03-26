import warnings

from autosteper.cage import Cage
from autosteper.cage import blk_list
from autosteper.checker import Checker
from autosteper.generator import Generator
from autosteper.optimizers import *
from autosteper.tools import get_low_e_ranks


class AutoSteper():
    def __init__(self, para: dict):
        self.cage = Cage(pristine_path=para['pristine_path'])
        self.cage.set_workbase(root=para['root'])
        self.generator = Generator(gen_para=para['gen_para'], pst_cage=self.cage)
        self.optimizer = switch_optimizers(mode=para['opt_mode'], para=para['opt_para'])
        self.optimizer.checker = Checker(group=self.generator.group, pst_cage=self.cage)

        if 'run_para' in para.keys():
            self.start = para['run_para']['start']
            self.stop = para['run_para']['stop']

            if 'step' not in para['run_para'].keys():
                self.step = 1
            else:
                self.step = para['run_para']['step']

            if 'wht_list_para' not in para['run_para'].keys():
                self.wht_list_para = {
                    'mode': 'rank_and_value',
                    'rank': 200,
                    'value': 1
                }
            else:
                self.wht_list_para = para['run_para']['wht_list_para']

        if 'random_para' in para.keys():
            self.addon_list = para['random_para']['addon_list']
            self.random_num = para['random_para']['random_num']
            self.try_times = para['random_para']['try_times']

        if 'pre_scan_para' in para.keys():
            ps_para = para['pre_scan_para']
            self.cage.ps_num_list = range(ps_para['start_ps_para'], ps_para['final_ps_para'] + 1)
            self.cage.calc = ps_para['calculator']
            self.ps_cut_para = ps_para['ps_cut_para']

        if 'blk_para' in para.keys():
            self.cage.blk_list = blk_list(size=self.cage.size, blk_para=para['blk_para'])

    def get_parent_info(self, all_parent_info: dict = None):
        parent_info_path = os.path.abspath('parent_info.pickle')
        passed_info_path = os.path.abspath('passed_info.pickle')
        passed_info = pd.read_pickle(passed_info_path)
        passed_e = passed_info['energy']
        passed_names = passed_info['name']
        parent_info = dict()
        if all_parent_info == None:
            for idx, a_name in enumerate(passed_names):
                parent_info.update({a_name: [passed_e[idx]]})
        else:
            for idx, a_name in enumerate(passed_names):
                parent_info.update({a_name: [all_parent_info[a_name][0], passed_e[idx]]})
        parent_info_df = pd.DataFrame(parent_info)
        parent_info_df.to_pickle(parent_info_path)
        return passed_info_path, parent_info_path

    def _post_pre_scan(self):
        def _copy_file():
            name = name_list[e_old_rank_map[a_rank]]
            os.symlink(src=os.path.abspath(os.path.join(f'{self.optimizer.path_raw}', f'{name}.xyz')),
                       dst=os.path.abspath(os.path.join('post_pre_scan_raw', f'{name}.xyz')))

        os.makedirs('post_pre_scan_raw', exist_ok=True)
        name_list = list(self.generator.ps_map.keys())
        e_series = pd.Series(self.generator.ps_map.values()).sort_values()
        self.generator.ps_map = dict()
        e_old_rank_map = dict(zip(range(len(e_series)), e_series.index))
        for a_rank in get_low_e_ranks(e_arr=np.array(e_series), para=self.ps_cut_para):
            _copy_file()
        self.optimizer.path_raw = os.path.abspath(r'post_pre_scan_raw')

    def random(self):
        for an_add_num in self.addon_list:
            if an_add_num < 5:
                self.cage.set_add_num(an_add_num)
                gen_out_path = f'{self.cage.name}_{self.cage.add_num}_addons.out'
                self.generator.gen_seq(gen_out_path=gen_out_path,
                                       mode='base')
                with open(gen_out_path, 'r') as file:
                    gen_num = len(file.readlines())
                if gen_num > self.random_num:
                    os.remove(gen_out_path)
                    os.chdir(self.cage.workbase)
                else:
                    print(f'The total number of {self.cage.add_num} addon isomers for '
                          f'cage in {self.cage.pristine_path} is less than random opt number, '
                          f'random opt procedure is terminated for this system.')
                    os.chdir(self.cage.workbase)
                    shutil.rmtree(self.cage.addon_path)
                    continue
            opt_times = 0
            while True:
                self.cage.set_add_num(an_add_num)
                gen_out_path = f'{self.cage.name}_{self.cage.add_num}_addons.out'
                self.generator.gen_seq(gen_out_path=gen_out_path,
                                       mode='random',
                                       random_num=self.random_num)
                self.optimizer.set_folders()
                self.generator.build(is_first=True,
                                     gen_out_path=gen_out_path,
                                     dump_folder=self.optimizer.path_raw,
                                     dyn_cage=self.cage,
                                     prev_xyz_path=self.cage.pristine_path)
                random_status = self.optimizer.opt()
                if random_status != 0:
                    os.chdir(self.cage.workbase)
                    shutil.rmtree(self.cage.addon_path)
                    if random_status == -2:
                        warnings.warn(
                            f'Something wrong happened while optimizing isomers in {self.optimizer.path_cooking}.\n'
                            f'The whole batch will be discarded.')
                    elif random_status == -1:
                        warnings.warn(f'All jobs failed while optimizing isomers in {self.optimizer.path_cooking}.\n'
                                      f'The whole batch will be discarded.')
                else:
                    break
                if opt_times == self.try_times - 1:
                    print(f'Random opt procedure has performed {self.try_times} times, '
                          f'still get wrong or failed status. It\'s been terminated, please check this system.')
                    break
                else:
                    opt_times = opt_times + 1

    def _first_step(self):
        self.cage.set_add_num(self.start)

        gen_out_path = f'{self.cage.name}_{self.cage.add_num}_addons.out'
        self.generator.gen_seq(gen_out_path=gen_out_path, mode='base')
        self.optimizer.set_folders()

        self.generator.build(is_first=True,
                             gen_out_path=gen_out_path,
                             dump_folder=self.optimizer.path_raw,
                             dyn_cage=self.cage,
                             prev_xyz_path=self.cage.pristine_path,
                             )
        if self.cage.is_pre_scan:
            self._post_pre_scan()

        step_status = self.optimizer.opt()
        if step_status == -1:
            return -1
        elif step_status == -2:
            raise RuntimeError(f'AutoSteper {self.mode} mode do not support '
                               f'{self.optimizer.deal_wrong_mode} deal wrong mode.')
        else:
            if self.cage.has_blk_list:
                self.cage.blk_list.clct_failed(r'status_info.pickle')
                if self.cage.blk_list.clct_unstb:
                    self.cage.blk_list.clct_unstable(info_path=self.prev_passed_info)
                else:
                    self.cage.blk_list.blk_list_arr = self.cage.blk_list.failed_arr

            self.prev_passed_info, self.prev_parent_info = self.get_parent_info()
            return 0

    def restart_first_step(self, restart_add_num):
        for i in range(restart_add_num, self.stop + 1):
            a_prev_addon_path = os.path.join(self.cage.workbase, f'{i}addons')
            if os.path.exists(a_prev_addon_path):
                shutil.rmtree(a_prev_addon_path)
        self.cage.set_add_num(restart_add_num - 1)
        self.prev_passed_info, self.prev_parent_info = \
            os.path.abspath('passed_info.pickle'), os.path.abspath('parent_info.pickle')
        return 0

    def base(self, add_num: int):
        self.start = add_num
        step_status = self._first_step()
        if step_status == 0:
            print("Normal Termination of AutoSteper.")
        else:
            print("All jobs failed.")

    def _take_a_step(self):
        def _step_build_unit():
            prev_xyz_path = prev_passed_info['xyz_path'][a_rank]
            prev_name = prev_passed_info['name'][a_rank]
            prev_seq, prev_addon_set, _ = name2seq(name=prev_name, cage_size=self.cage.size)
            sub_nauty_o = os.path.join(sub_nauty_path, f'{prev_name}.out')
            self.generator.gen_seq(mode='step', gen_out_path=sub_nauty_o, prev_seq=prev_seq)
            self.all_parent_info = self.generator.build(is_first=False,
                                                        gen_out_path=sub_nauty_o,
                                                        dump_folder=self.optimizer.path_raw,
                                                        dyn_cage=self.cage,
                                                        prev_xyz_path=prev_xyz_path,
                                                        parent_info=self.all_parent_info,
                                                        prev_addon_set=prev_addon_set,
                                                        parent_name=prev_name,
                                                        )

        self.cage.set_add_num(self.new_add_num)
        self.optimizer.set_folders()
        self.all_parent_info = {}
        sub_nauty_path = 'sub_nauty'
        os.makedirs(sub_nauty_path, exist_ok=True)

        prev_passed_info = pd.read_pickle(self.prev_passed_info)
        for a_rank in get_low_e_ranks(e_arr=np.array(prev_passed_info['energy']), para=self.wht_list_para):
            _step_build_unit()
        prev_parent_info = pd.read_pickle(self.prev_parent_info)
        wht_list_names = prev_passed_info['name'][: a_rank + 1]
        new_prev_flat = prev_parent_info[wht_list_names]
        new_prev_flat.to_pickle(self.prev_parent_info)

        all_parent_info_df = pd.DataFrame(self.all_parent_info)
        all_parent_info_df.to_pickle('all_parent_info.pickle')

        if self.cage.is_pre_scan:
            self._post_pre_scan()

        step_status = self.optimizer.opt()
        if step_status == -2:
            raise RuntimeError(f'AutoSteper {self.mode} mode do not support '
                               f'{self.optimizer.deal_wrong_mode} deal wrong mode.')
        elif step_status == -1:
            return step_status
        self.prev_passed_info, self.prev_parent_info = self.get_parent_info(all_parent_info=self.all_parent_info)
        if self.cage.has_blk_list:
            self.cage.blk_list.clct_failed(r'status_info.pickle')
            if self.cage.blk_list.clct_unstb:
                self.cage.blk_list.clct_unstable(info_path=self.prev_passed_info)
            else:
                self.cage.blk_list.blk_list_arr = self.cage.blk_list.failed_arr
        return step_status

    # re-opt jobs in restart_add_num
    def restart(self, restart_add_num: int):
        step_status = self.restart_first_step(restart_add_num)
        while step_status == 0:
            self.new_add_num = self.cage.add_num + self.step
            if self.stop != None:
                if self.new_add_num > self.stop:
                    print("Normal Termination of AutoSteper.")
                    break
            step_status = self._take_a_step()
            if step_status == -1:
                print(f"All jobs failed when optimizing {self.cage.add_num} addon system.")
                break

    def run(self):
        step_status = self._first_step()
        if step_status == -1:
            print("AutoSteper failed on first step.")
        else:
            while step_status == 0:
                self.new_add_num = self.cage.add_num + self.step
                if self.stop != None:
                    if self.new_add_num > self.stop:
                        print("Normal Termination of AutoSteper.")
                        break
                step_status = self._take_a_step()
                if step_status == -1:
                    print(f"All jobs failed when optimizing {self.cage.add_num} addon system.")
                    break
