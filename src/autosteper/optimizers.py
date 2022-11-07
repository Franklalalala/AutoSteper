import math
import os
import shutil
from typing import Union

import numpy as np
import pandas as pd
from ase.io import read, write
from ase.io.gaussian import read_gaussian_out, write_gaussian_in
from ase.units import Hartree
from autosteper.cage import name2seq
from autosteper.tools import get_low_e_ranks
from dpdispatcher import Task, Submission, Machine, Resources


def switch_optimizers(mode: str, para: Union[list, dict]):
    if mode == 'xtb':
        return XTB_Optimizer(para)
    elif mode == 'ase':
        return ASE_Optimizer(para)
    elif mode == 'gaussian':
        return Gaussian_Optimizer(para)
    elif mode == 'multi':
        return Multi_Optimizer(para)


# job status code:
# 0: Normal termination.
# -1: Failed.
# -2: Wrong.

# step status code:
# 0: Normal termination.
# -1: All jobs failed.
# -2: At least one job end wrong.

# Energy should be in unit eV

class Optimizer():
    def __init__(self, opt_para: dict):
        # self.deal_wrong_mode = opt_para['deal_wrong_mode']
        self.mode = None
        self.cmd_list = opt_para['cmd_list']
        self.out_list = opt_para['out_list']
        self.deal_wrong_mode = opt_para['deal_wrong_mode']
        self.resrc_para = opt_para['resrc_para']
        self.mach_para = opt_para['mach_para']
        self.checker = None
        self.is_odd = False
        if 'has_parity' in opt_para.keys():
            self.cage_size = opt_para['cage_size']
            self.has_parity = opt_para['has_parity']
        else:
            self.has_parity = False

    def set_folders(self):
        self.path_raw = os.path.abspath('raw')
        os.makedirs(exist_ok=True, name=self.path_raw)
        self.path_cooking = os.path.abspath('cooking')
        os.makedirs(exist_ok=True, name=self.path_cooking)
        self.path_cooked = os.path.abspath('cooked')
        os.makedirs(exist_ok=True, name=self.path_cooked)

    def simple_chg_parity(self):
        a_name = os.path.splitext(os.listdir(self.path_raw)[0])[0]
        _, _0, an_arr = name2seq(name=a_name, cage_size=self.cage_size)
        add_num = an_arr.sum()
        if add_num % 2 == 1:
            self.is_odd = True
        else:
            self.is_odd = False

    def pack(self):
        pass

    def unpack(self):
        pass

    def run_a_batch(self):
        run_cwd_ = os.getcwd()

        def _run_unit():
            try:
                submission.run_submission()
            except Exception as e:
                print(os.getcwd())
                if os.path.exists('0/xyz/4do2a5oxwe0w.xyz'):
                    print('11111111111111111')
                with open(r'./opt_error.log', 'a') as f:
                    f.write(str(e) + '\n' + str(self.task_list) + '\n')

        self.mach = Machine(batch_type=self.mach_para['batch_type'],
                            context_type=self.mach_para['context_type'],
                            remote_root=self.mach_para['remote_root'],
                            remote_profile=self.mach_para['remote_profile'],
                            local_root=self.path_cooking)
        self.resrc = Resources(number_node=self.resrc_para['number_node'],
                               cpu_per_node=self.resrc_para['cpu_per_node'],
                               gpu_per_node=self.resrc_para['gpu_per_node'],
                               group_size=self.resrc_para['group_size'],
                               queue_name=self.resrc_para['queue_name'],
                               envs=self.resrc_para['envs']
                               )
        if self.resrc_para['sub_batch_size'] == None:
            submission = Submission(work_base=self.path_cooking,
                                    machine=mach,
                                    resources=resrc,
                                    task_list=self.task_list,
                                    forward_common_files=[],
                                    backward_common_files=[]
                                    )
            _run_unit()
        else:
            sub_batch_size = self.resrc_para['sub_batch_size']
            num_groups = math.ceil(len(self.task_list) / sub_batch_size)
            for i in range(num_groups):
                cursor = i * sub_batch_size
                a_task_list = self.task_list[cursor:cursor + sub_batch_size]
                submission = Submission(work_base=self.path_cooking,
                                        machine=self.mach,
                                        resources=self.resrc,
                                        task_list=a_task_list,
                                        forward_common_files=[],
                                        backward_common_files=[]
                                        )
                _run_unit()
        os.chdir(run_cwd_)

    def prep_unpack(self):
        self.names = []
        self.passed_names = []
        self.e_list = []
        self.path_list = []
        self.nimages_list = []
        self.status_codes = []

    def post_unpack(self):
        info = pd.DataFrame({'name': self.passed_names, 'energy': self.e_list, 'xyz_path': self.path_list,
                             'nimages': self.nimages_list})
        passed_info = info.sort_values(by='energy')
        passed_info.index = sorted(passed_info.index)
        passed_info.to_pickle(path='passed_info.pickle')
        status_info = pd.DataFrame(dict(zip(self.names, [[i] for i in self.status_codes])))
        status_info.to_pickle('status_info.pickle')

    def deal_wrong_complete(self):
        deal_wrong_cwd_ = os.getcwd()
        for self.resrc_para['group_size'], self.resrc_para['sub_batch_size'] in self.resrc_para['wrong_rcs_para']:
            self.new_raws = os.path.join(deal_wrong_cwd_, 'deal_wrong', 'raw')
            os.makedirs(self.new_raws)
            for i_code, a_code in enumerate(self.status_codes):
                if a_code == -2:
                    a_wrong_name = self.names[i_code]
                    shutil.copy(src=os.path.join(deal_wrong_cwd_, 'raw', a_wrong_name + '.xyz'),
                                dst=os.path.join(self.new_raws, a_wrong_name + '.xyz'))
            os.chdir('deal_wrong')
            self.set_folders()
            self.pack()
            self.run_a_batch()
            self.unpack()
            for a_job in os.listdir('cooking'):
                shutil.copytree(src=os.path.join('cooking', a_job), dst=os.path.join(deal_wrong_cwd_, 'cooking', a_job),
                                dirs_exist_ok=True)
            os.chdir(deal_wrong_cwd_)
            new_path_list = []
            for a_xyz_path in self.path_list:
                a_new_path = os.path.abspath(os.path.join('cooked', os.path.basename(a_xyz_path)))
                shutil.copy(src=a_xyz_path, dst=a_new_path)
                new_path_list.append(a_new_path)
            old_info = pd.read_pickle('passed_info.pickle')
            new_info = pd.read_pickle(os.path.join(deal_wrong_cwd_, 'deal_wrong', 'passed_info.pickle'))
            new_info['xyz_path'] = new_path_list
            info = pd.concat([old_info, new_info], ignore_index=True)
            sorted_info = info.sort_values(by='energy')
            sorted_info.index = sorted(sorted_info.index)
            sorted_info.to_pickle(path='passed_info.pickle')
            old_status_info = pd.read_pickle('status_info.pickle')
            for i_code, a_code in enumerate(self.status_codes):
                if a_code != -2:
                    old_status_info[self.names[i_code]] = a_code
            old_status_info.to_pickle('status_info.pickle')
            shutil.rmtree('deal_wrong')

    def opt(self):
        if len(os.listdir(self.path_raw)) == 0:
            raise RuntimeError(f'There are no isomers in {self.path_raw}.\n'
                               f'Please check this system.')
        self.pack()
        self.run_a_batch()
        self.unpack()
        if -2 in self.status_codes:
            if self.deal_wrong_mode == 'Complete':
                self.deal_wrong_complete()
                raise RuntimeError(f'Something wrong happened while optimizing isomers in {self.path_cooking}.\n'
                                   f'Re-opted recursively to get a minimize the wrong jobs.')
            elif self.deal_wrong_mode == 'Report':
                raise RuntimeError(f'Something wrong happened while optimizing isomers in {self.path_cooking}.')
            elif self.deal_wrong_mode == 'Tough':
                return -2
        if self.checker:
            self.status_codes = self.checker.check(passed_info_path='passed_info.pickle',
                                                   status_info_path='status_info.pickle')
        if 0 not in self.status_codes:
            return -1
        else:
            return 0


class XTB_Optimizer(Optimizer):
    def __init__(self, opt_para: dict):
        super(XTB_Optimizer, self).__init__(opt_para=opt_para)
        self.mode = 'xtb'
        self.opted_xyz_name = 'xtbopt.xyz'

    def pack(self):
        if self.has_parity:
            self.simple_chg_parity()
        if self.is_odd:
            self.cmd_list.append('--uhf 1')
        self.task_list = []
        for item in os.listdir(self.path_raw):
            isomer_name = os.path.splitext(item)[0]
            os.makedirs(os.path.join(self.path_cooking, isomer_name), exist_ok=True)
            shutil.copy(os.path.join(self.path_raw, item), os.path.join(self.path_cooking, isomer_name, item))
            self.cmd_list.append(item)
            a_task = Task(command=' '.join(self.cmd_list),
                          task_work_path=f"{isomer_name}/",
                          forward_files=[item],
                          backward_files=self.out_list
                          )
            self.task_list.append(a_task)
            del self.cmd_list[-1]
        if self.is_odd:
            del self.cmd_list[-1]

    def unpack(self):
        cwd_ = os.getcwd()
        self.prep_unpack()
        for a_job in os.listdir(self.path_cooking):
            os.chdir(os.path.join(self.path_cooking, a_job))
            self.names.append(a_job)
            if not os.path.exists(self.opted_xyz_name):
                self.status_codes.append(-2)
                os.chdir(cwd_)
                continue
            self.status_codes.append(0)
            self.passed_names.append(a_job)
            new_path = os.path.join(self.path_cooked, a_job + '.xyz')
            shutil.copy(src='xtbopt.xyz', dst=new_path)
            with open(self.opted_xyz_name, 'r') as f:
                _ = f.readline()
                e_line = f.readline()
                energy = float(e_line.split()[1]) * Hartree
            with open(r'xtbopt.log', "r") as file:
                lines = file.readlines()
                natoms = int(lines[0])
                nimages = len(lines) // (natoms + 2)
            self.e_list.append(energy)
            self.nimages_list.append(nimages)
            self.path_list.append(new_path)
            os.chdir(cwd_)
        self.post_unpack()


class ASE_Optimizer(Optimizer):
    def __init__(self, opt_para: dict):
        super(ASE_Optimizer, self).__init__(opt_para=opt_para)
        self.mode = 'ase'
        self.opted_xyz_name = 'opt.xyz'
        self.ase_para = opt_para['ase_para']

    def pack(self):
        cwd_ = os.getcwd()
        self.task_list = []
        num_raws = len(os.listdir(self.path_raw))
        num_per_worker = math.ceil(num_raws / self.ase_para['num_pll'])
        model_path = os.path.abspath(self.ase_para['model_path'])

        if isinstance(self.ase_para['py_script'], str):
            py_script = os.path.abspath(self.ase_para['py_script'])
        else:
            new_list = []
            for a_script in self.ase_para['py_script']:
                new_list.append(os.path.abspath(a_script))
            py_script = new_list

        for i in range(min(self.ase_para['num_pll'], num_raws)):
            os.chdir(self.path_cooking)
            cursor = i * num_per_worker
            sub_raw = str(i)
            os.makedirs(sub_raw, exist_ok=True)
            os.chdir(sub_raw)

            if isinstance(py_script, list):
                py_script = py_script[i % len(py_script)]
            else:
                py_script = py_script
            pll_file_name = os.path.basename(py_script)
            shutil.copy(src=py_script, dst=pll_file_name)
            model_name = os.path.basename(model_path)
            shutil.copy(src=model_path, dst=model_name)

            sub_raw_xyz = 'xyz'
            os.makedirs(sub_raw_xyz, exist_ok=True)
            for a_raw_xyz in os.listdir(self.path_raw)[cursor: cursor + num_per_worker]:
                shutil.copy(src=os.path.join(self.path_raw, a_raw_xyz), dst=os.path.join(sub_raw_xyz, a_raw_xyz))

            start_node = self.ase_para['base_node'] + self.ase_para['cpu_per_worker'] * i
            end_node = self.ase_para['base_node'] + self.ase_para['cpu_per_worker'] * (i + 1) - 1
            real_cmd_list = [f'taskset -c {start_node}-{end_node}', *self.cmd_list, pll_file_name]
            a_task = Task(command=' '.join(real_cmd_list),
                          task_work_path=f"{i}/",
                          forward_files=[sub_raw_xyz, pll_file_name, model_name],
                          backward_files=self.out_list
                          )
            self.task_list.append(a_task)
            os.chdir(cwd_)

    def unpack(self):
        cwd_ = os.getcwd()
        self.prep_unpack()
        for a_folder in os.listdir(self.path_cooking):
            os.chdir(os.path.join(self.path_cooking, a_folder, 'cooked'))
            for a_job in os.listdir('./'):
                os.chdir(a_job)
                self.names.append(a_job)
                if not os.path.exists(self.opted_xyz_name):
                    self.status_codes.append(-2)
                    os.chdir(cwd_)
                    continue

                self.status_codes.append(0)
                self.passed_names.append(a_job)
                new_path = os.path.join(self.path_cooked, a_job + '.xyz')
                shutil.copy(src=self.opted_xyz_name, dst=new_path)

                self.path_list.append(new_path)
                with open(r'opt.log', 'r') as xyz_f:
                    lines = xyz_f.readlines()
                    energy_line = lines[-1]
                    energy = float(energy_line.split()[-2].split('*')[0]) * Hartree
                    nimages = len(lines) - 2

                self.e_list.append(energy)
                self.nimages_list.append(nimages)
                os.chdir('./..')
            os.chdir(cwd_)
        self.post_unpack()


class Gaussian_Optimizer(Optimizer):
    def __init__(self, opt_para: dict):
        super(Gaussian_Optimizer, self).__init__(opt_para=opt_para)
        self.mode = 'gaussian'
        self.gau_para = opt_para['gau_para']

    def pack(self):
        if self.has_parity:
            self.simple_chg_parity()
            if self.is_odd:
                multi = 1
            else:
                multi = 0
        else:
            multi = self.gau_para['multi']

        self.task_list = []
        self.cmd_list.append('gau.gjf')
        for item in os.listdir(self.path_raw):
            isomer_name = os.path.splitext(item)[0]
            atoms = read(os.path.join(self.path_raw, item))
            os.makedirs(os.path.join(self.path_cooking, isomer_name), exist_ok=True)
            with open(file=os.path.join(self.path_cooking, isomer_name, 'gau.gjf'), mode='w') as f:
                write_gaussian_in(fd=f,
                                  atoms=atoms,
                                  properties=[' '],
                                  method='',
                                  basis=self.gau_para['cmd_line'],
                                  nprocshared=str(self.gau_para['cpu_per_worker']),
                                  mem=self.gau_para['mem'],  # 10GB by default
                                  mult=multi,
                                  charge=self.gau_para['charge'],
                                  chk='gau.chk'
                                  )
            a_task = Task(command=' '.join(self.cmd_list),
                          task_work_path=f"{isomer_name}/",
                          forward_files=['gau.gjf'],
                          backward_files=self.out_list
                          )
            self.task_list.append(a_task)

    def unpack(self):
        cwd_ = os.getcwd()
        self.prep_unpack()
        for a_job in os.listdir(self.path_cooking):
            os.chdir(os.path.join(self.path_cooking, a_job))
            self.names.append(a_job)
            if not os.path.exists(r'gau.log'):
                self.status_codes.append(-2)
                os.chdir(cwd_)
                continue
            else:
                with open(file='gau.log', mode='r') as fd:
                    end_line = fd.readlines()[-1]
                if not end_line.startswith(r' Normal termination'):
                    self.status_codes.append(-1)
                    os.chdir(cwd_)
                    continue
            self.status_codes.append(0)
            self.passed_names.append(a_job)
            traj = read_gaussian_out(
                fd='gau.log')  # The read_gaussian_out function is slightly changed to get all traj.

            write(filename='gau_traj.log', images=traj, format='xyz')
            nimages = len(traj)

            last_image = traj[-1]
            new_path = os.path.join(self.path_cooked, a_job + '.xyz')
            write(filename=new_path, images=last_image, format='xyz')
            energy = last_image.get_potential_energy()

            self.e_list.append(energy)
            self.nimages_list.append(nimages)
            self.path_list.append(new_path)
            os.chdir(cwd_)
        self.post_unpack()


class Multi_Optimizer(Optimizer):
    def __init__(self, opt_para: list):
        self.mode = 'multi'
        self.opt_para = opt_para

    def opt(self):
        cwd_multi = os.getcwd()
        for idx, a_para in enumerate(self.opt_para):
            a_workbase = f'opt_{idx}'
            os.makedirs(a_workbase, exist_ok=True)
            os.chdir(a_workbase)
            an_optimizer = switch_optimizers(mode=a_para['opt_mode'], para=a_para)
            an_optimizer.set_folders()
            if idx == 0:
                shutil.copytree(src=os.path.join(cwd_multi, 'raw'), dst=an_optimizer.path_raw, dirs_exist_ok=True)
            else:
                wht_list_para = a_para['wht_list_para']
                pre_passed_info = pd.read_pickle(os.path.join(cwd_multi, f'opt_{idx - 1}', 'passed_info.pickle'))
                if 'nimg_th' in wht_list_para.keys():
                    wht_list_para.update({'nimages': pre_passed_info['nimages']})

                for a_rank in get_low_e_ranks(e_arr=np.array(pre_passed_info['energy']), para=wht_list_para):
                    a_xyz_path = pre_passed_info['xyz_path'][a_rank]
                    a_xyz_name = os.path.basename(a_xyz_path)
                    shutil.copy(src=a_xyz_path, dst=os.path.join(an_optimizer.path_raw, a_xyz_name))
            status = an_optimizer.opt()
            if status == -1:
                return -1
            os.chdir(cwd_multi)
        for file in os.listdir(a_workbase):
            file_path = os.path.join(a_workbase, file)
            if os.path.isfile(file_path):
                shutil.copy(src=file_path, dst=file)
        return 0
