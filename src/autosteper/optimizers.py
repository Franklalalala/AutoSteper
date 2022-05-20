import os
import shutil
import math
from autosteper.cage import Cage
from autosteper.checker import Checker
from dpdispatcher import Task, Submission, Machine, Resources
from ase.io import read, write
from ase.optimize import *
from ase.io import read
from ase.io.gaussian import write_gaussian_in
import subprocess


class Optimizer():
    def __init__(self, opt_para: dict, checker: Checker, cage: Cage):
        self.init_cycle = opt_para['init_cycle']
        self.is_Opt_Twice = opt_para['is_Opt_Twice']
        self.deal_wrong_mode = opt_para['deal_wrong_mode']
        if self.is_Opt_Twice:
            self.final_cycle = opt_para['final_cycle']
        self.checker = checker
        self.cage = cage
        self.mode = None
        self.cmd_list = None
        self.calc = None

    def set_init_folders(self):
        self.path_raw_init = f'raw_{self.init_cycle}'
        os.makedirs(exist_ok=True, name=self.path_raw_init)
        self.path_opt_init = f'opt_{self.init_cycle}'
        os.makedirs(exist_ok=True, name=self.path_opt_init)

    def set_final_folders(self):
        self.path_raw_final = 'raw_final'
        os.makedirs(exist_ok=True, name=self.path_raw_final)
        self.path_opt_final = 'opt_final'
        os.makedirs(exist_ok=True, name=self.path_opt_final)

    def run_a_batch(self):
        pass

    def opt_twice(self):
        pass

    def opt_once(self):
        pass

    def opt(self):
        if self.is_Opt_Twice:
            status = self.opt_twice()
        else:
            status = self.opt_once()
        return status


class XTB_Optimizer(Optimizer):
    def __init__(self, opt_para: dict, checker: Checker, cage: Cage):
        super(XTB_Optimizer, self).__init__(opt_para=opt_para, checker=checker, cage=cage)
        self.xtb_path = opt_para['xtb_path']
        self.cmd_list = [self.xtb_path, *opt_para['cmd_list']]
        self.out_list = opt_para['out_list']
        self.resrc_para = opt_para['resrc_para']
        self.mach_para = opt_para['mach_para']
        self.sub_batch_size = opt_para['sub_batch_size']
        self.mode = 'xtb'

    # group_size is how many jobs in a single .sub file uploaded.
    # the old version of this wrapper upload all raw files in path_source in a single submission,
    # which will be divided into multiple groups, in each group, there will be group size jobs.
    # if one file failed, the whole submission will terminated, and the results of already calculated jobs will not be transfered to workbase.
    # It's wasted, so I propose the parameter, sub_batch_size, which will divide the one submission to several small submissions
    def run_a_batch(self, path_source: str, path_destination: str, cmd_list: list):
        path_destination = os.path.abspath(path_destination)
        self.mach = Machine(batch_type=self.mach_para['batch_type'],
                            context_type=self.mach_para['context_type'],
                            remote_root=self.mach_para['remote_root'],
                            remote_profile=self.mach_para['remote_profile'],
                            local_root=path_destination)
        self.resrc = Resources(number_node=self.resrc_para['number_node'],
                               cpu_per_node = self.resrc_para['cpu_per_node'],
                               gpu_per_node = self.resrc_para['gpu_per_node'],
                               group_size = self.resrc_para['group_size'],
                               queue_name = self.resrc_para['queue_name'],
                               envs = self.resrc_para['envs']
                               )
        task_list = []
        for item in os.listdir(path_source):
            isomer_name = os.path.splitext(item)[0]
            os.makedirs(os.path.join(path_destination, isomer_name), exist_ok=True)
            path_item = os.path.join(path_source, item)
            shutil.copy(path_item, os.path.join(path_destination, isomer_name, item))
            cmd_list.append(item)
            a_task = Task(command=' '.join(cmd_list),
                          task_work_path=f"{isomer_name}/",
                          forward_files=[item],
                          backward_files=self.out_list
                          )
            task_list.append(a_task)
            del cmd_list[-1]

        if self.sub_batch_size == None:
            submission = Submission(work_base=path_destination,
                                    machine=self.mach,
                                    resources=self.resrc,
                                    task_list=task_list,
                                    forward_common_files=[],
                                    backward_common_files=[]
                                    )
            try:
                submission.run_submission()
            except Exception as e:
                with open('opt_error.log', 'a') as f:
                    f.write(str(e) + '\n' + str(task_list) + '\n')

        else:
            num_groups = math.ceil(len(task_list) / self.sub_batch_size)
            for i in range(num_groups):
                cursor = i * self.sub_batch_size
                a_task_list = task_list[cursor:cursor + self.sub_batch_size]
                submission = Submission(work_base=path_destination,
                                        machine=self.mach,
                                        resources=self.resrc,
                                        task_list=a_task_list,
                                        forward_common_files=[],
                                        backward_common_files=[]
                                        )
                try:
                    submission.run_submission()
                except Exception as e:
                    with open('opt_error.log', 'a') as f:
                        f.write(str(e) + '\n' + str(task_list) + '\n')

    def opt_twice(self):
        cmd_list = self.cmd_list.copy()
        if self.cage.add_num % 2 == 1:
            cmd_list.append('--uhf 1')
        cmd_list.extend([f'--cycles {self.init_cycle}'])
        self.run_a_batch(path_source=self.path_raw_init, path_destination=self.path_opt_init, cmd_list=cmd_list)
        del cmd_list[-1]

        self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_init, is_init=True, init_cycle=self.init_cycle)

        if os.stat('wrong_paths').st_size != 0:
            if self.deal_wrong_mode == 'complete':
                # re-opt wrong_path with small sub batch size
                self.deal_wrong_complete()
            elif self.deal_wrong_mode == 'report':
                print(f'Something wrong happened while optimizing isomers in {self.path_raw_init}.')
                return 0
            elif self.deal_wrong_mode == 'tough':
                return 2

        # if all jobs are failed or wronged, stop opt
        if os.stat('yes_paths').st_size == 0:
            self.combine_log()
            if os.stat('init_yes_paths').st_size == 0:
                return 0
            else:
                return 1
        else:
            self.set_final_folders()

            with open('yes_paths', 'r') as f:
                for a_log in f.readlines():
                    src = os.path.join(os.path.split(a_log)[0], 'xtbopt.xyz')
                    dst = os.path.join(self.path_raw_final, os.path.basename(os.path.split(a_log)[0]) + '.xyz')
                    os.symlink(src=src, dst=dst)

            cmd_list.extend([f'--cycles {self.final_cycle}'])
            self.run_a_batch(path_source=self.path_raw_final, path_destination=self.path_opt_final, cmd_list=cmd_list)
            # can deal with wrong here, but not meaningful.
            self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_final, is_init=False)
            self.combine_log()
            if os.stat('yes_paths').st_size == 0 and os.stat('init_yes_paths').st_size == 0:
                return 0
            else:
                return 1

    def opt_once(self):
        cmd_list = self.cmd_list.copy()
        if self.cage.add_num % 2 == 1:
            cmd_list.append('--uhf 1')
        cmd_list.extend([f'--cycles {self.init_cycle}'])
        self.run_a_batch(path_source=self.path_raw_init, path_destination=self.path_opt_init, cmd_list=cmd_list)

        self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_init, is_init=True, init_cycle=self.init_cycle)

        if os.stat('wrong_paths').st_size != 0:
            if self.deal_wrong_mode == 'complete':
                # re-opt wrong_path with small sub batch size
                self.deal_wrong_complete()
            elif self.deal_wrong_mode == 'report':
                print(f'Something wrong happened while optimizing isomers in {self.path_raw_init}.')
                return 0
            elif self.deal_wrong_mode == 'tough':
                return 2


        # if all jobs are failed or wronged, stop opt
        if (os.stat('yes_paths').st_size != 0) or (os.stat('init_yes_paths').st_size != 0):
            return 1
        else:
            return 0

    def combine_log(self):
        cwd_ = os.getcwd()
        os.makedirs(name='combined_log', exist_ok=True)
        if os.path.exists('opt_final'):
            for an_init in os.listdir(f'opt_{self.init_cycle}'):
                if an_init not in os.listdir('opt_final'):
                    os.symlink(
                        src=os.path.join(cwd_, f'opt_{self.init_cycle}', an_init, 'xtbopt.log'),
                        dst=os.path.join('combined_log', f'{an_init}.log')
                    )
                else:
                    init_path = os.path.join('combined_log', f'{an_init}.log')
                    shutil.copy(
                        src=os.path.join(f'opt_{self.init_cycle}', an_init, 'xtbopt.log'),
                        dst=init_path
                    )
                    final_path = os.path.join('opt_final', an_init, 'xtbopt.log')

                    with open(init_path, 'a') as init_file, open(final_path, 'r') as final_file:
                        atom_num = int(final_file.readline().strip())
                        for idx, a_line in enumerate(final_file.readlines()):
                            if idx < atom_num + 1:
                                continue
                            else:
                                init_file.write(a_line)

        else:
            for an_init in os.listdir(f'opt_{self.init_cycle}'):
                os.symlink(
                    src=os.path.join(cwd_, f'opt_{self.init_cycle}', an_init, 'xtbopt.log'),
                    dst=os.path.join('combined_log', f'{an_init}.log')
                )


    def deal_wrong_complete(self):
        def _wrong_unit():
            os.makedirs(path_wrong_opt, exist_ok=True)
            with open('wrong_paths', 'r') as f:
                for a_wrong in f.readlines():
                    a_wrong = a_wrong.rstrip()
                    shutil.copy(src=a_wrong, dst=os.path.join(path_wrong_opt, os.path.split(a_wrong)[1]))
            self.run_a_batch(path_source=path_wrong_opt, path_destination=self.path_opt_init)
            shutil.rmtree(path_wrong_opt)
            self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_init, is_init=True,
                               init_cycle=self.init_cycle)

        # here is the situation of initial.
        path_wrong_opt = f'wrong_{self.init_cycle}'
        old_para = [self.resrc_para['group_size'], self.sub_batch_size]
        for self.resrc_para['group_size'], self.sub_batch_size in [[3, 15], [2, 8], [1, 4], [1, 2], [1, 1]]:
            _wrong_unit()
        self.resrc_para['group_size'], self.sub_batch_size = old_para




class ASE_Optimizer(Optimizer):
    def __init__(self, opt_para: dict, checker: Checker, cage: Cage):
        super(ASE_Optimizer, self).__init__(opt_para=opt_para, checker=checker, cage=cage)
        self.is_pll = opt_para['is_pll']
        if self.is_pll:
            if self.is_Opt_Twice:
                self.pll_unit_init = opt_para['pll_para']['pll_unit_init']
                self.pll_unit_final = opt_para['pll_para']['pll_unit_final']
            else:
                self.pll_unit_file = opt_para['pll_para']['pll_unit_file']
            self.num_worker = opt_para['pll_para']['num_worker']
            self.cpu_per_worker = opt_para['pll_para']['cpu_per_worker']
            self.base_node = opt_para['pll_para']['base_node']
        self.calc = opt_para['calculator']
        self.fmax = opt_para['fmax']
        self.ase_optimizer = opt_para['ase_optimizer']
        self.mode = 'ase'

    def run_a_batch(self, path_source: str, path_destination: str, steps: int):
        cwd_ = os.getcwd()
        path_source = os.path.abspath(path_source)
        path_destination = os.path.abspath(path_destination)
        for a_file in os.listdir(path_source):
            name = a_file.split('.')[0]
            sub_opt = os.path.join(path_destination, name)

            # if os.path.exists(sub_opt):
            #     continue

            os.makedirs(sub_opt, exist_ok=True)
            path_a_file = os.path.join(path_source, a_file)
            shutil.copy(path_a_file, os.path.join(sub_opt, a_file))
            os.chdir(sub_opt)
            atoms = read(a_file, format='xyz')
            atoms.calc = self.calc
            opt = self.ase_optimizer(atoms, logfile='opt.log')
            opt.run(steps=steps, fmax=self.fmax)
            write(filename='opt.xyz', images=atoms, format='xyz')
        os.chdir(cwd_)

    def run_a_batch_parallel(self, path_source: str, raw_num: int):
        path_source = os.path.abspath(path_source)
        cwd_ = os.getcwd()
        num_per_worker = math.ceil(raw_num / self.num_worker)
        res_list = [None] * self.num_worker
        path_temp = os.path.join(cwd_, 'temp')
        os.makedirs(path_temp, exist_ok=True)
        for i in range(self.num_worker):
            os.chdir(path_temp)
            cursor = i * num_per_worker
            sub_raw = str(i)
            os.makedirs(sub_raw, exist_ok=True)
            os.chdir(sub_raw)
            pll_file_name = os.path.basename(self.pll_unit_file)
            shutil.copy(src=self.pll_unit_file, dst=pll_file_name)
            sub_raw_xyz = 'xyz'
            os.makedirs(sub_raw_xyz, exist_ok=True)
            for a_raw_xyz in os.listdir(path_source)[cursor: cursor + num_per_worker]:
                shutil.copy(src=os.path.join(path_source, a_raw_xyz), dst=os.path.join(sub_raw_xyz, a_raw_xyz))
            start_node = self.base_node + self.cpu_per_worker * i
            end_node = self.base_node + self.cpu_per_worker * (i + 1) - 1
            res_list[i] = subprocess.Popen(f'taskset -c {start_node}-{end_node} python {pll_file_name}', shell=True)
        for a_proc in res_list:
            a_proc.wait()
        os.chdir(cwd_)
        shutil.rmtree(path_temp)

    def opt_once(self):
        raw_num = len(os.listdir(self.path_raw_init))
        if self.is_pll and raw_num >= self.num_worker:
            self.run_a_batch_parallel(path_source=self.path_raw_init, raw_num=raw_num)
        else:
            self.run_a_batch(path_source=self.path_raw_init, path_destination=self.path_opt_init, steps=self.init_cycle)
        self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_init, is_init=True, init_cycle=self.init_cycle)
        if os.stat('wrong_paths').st_size != 0:
            if self.deal_wrong_mode == 'report':
                print(f'Something wrong happened while optimizing isomers in {self.path_raw_init}.')
                return 0
            elif self.deal_wrong_mode == 'tough':
                return 2
        # if all jobs are failed, stop opt
        if os.stat('yes_paths').st_size == 0 and os.stat('init_yes_paths').st_size == 0:
            return 0
        else:
            return 1

    def opt_twice(self):
        self.run_a_batch(path_source=self.path_raw_init, path_destination=self.path_opt_init, steps=self.init_cycle)
        self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_init, is_init=True, init_cycle=self.init_cycle)
        if os.stat('wrong_paths').st_size != 0:
            if self.deal_wrong_mode == 'report':
                print(f'Something wrong happened while optimizing isomers in {self.path_raw_init}.')
                return 0
            elif self.deal_wrong_mode == 'tough':
                return 2
        # if all jobs are failed or wronged, stop opt
        if os.stat('yes_paths').st_size == 0:
            if os.stat('init_yes_paths').st_size == 0:
                return 0
            else:
                return 1
        else:
            self.set_final_folders()

            with open('yes_paths', 'r') as f:
                for a_log in f.readlines():
                    src = os.path.join(os.path.split(a_log)[0], 'opt.xyz')
                    dst = os.path.join(self.path_raw_final, os.path.basename(os.path.split(a_log)[0]) + '.xyz')
                    os.symlink(src=src, dst=dst)

            self.run_a_batch(path_source=self.path_raw_final, path_destination=self.path_opt_final, steps=self.final_cycle)
            # can deal with wrong here, but not meaningful.
            self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_final, is_init=False)
            if os.stat('yes_paths').st_size == 0 and os.stat('init_yes_paths').st_size == 0:
                return 0
            else:
                return 1


class Gaussian_Optimizer(Optimizer):
    def __init__(self, opt_para: dict, checker: Checker, cage: Cage):
        super(Gaussian_Optimizer, self).__init__(opt_para=opt_para, checker=checker, cage=cage)
        self.gau_path = opt_para['gau_path']
        self.cmd_list = [self.gau_path]
        self.out_list = opt_para['out_list']
        self.resrc_para = opt_para['resrc_para']
        self.mach_para = opt_para['mach_para']
        self.sub_batch_size = opt_para['sub_batch_size']
        self.mode = 'gaussian'
        self.mem = opt_para['mem']
        self.charge = opt_para['charge']
        self.job_type = opt_para['job_type']
        self.mult = None
        self.method = opt_para['method']
        self.other_para = opt_para['other_para']

    def run_a_batch(self, path_source: str, path_destination: str, cmd_list: list, cycles:int):
        path_destination = os.path.abspath(path_destination)
        self.mach = Machine(batch_type=self.mach_para['batch_type'],
                            context_type=self.mach_para['context_type'],
                            remote_root=self.mach_para['remote_root'],
                            remote_profile=self.mach_para['remote_profile'],
                            local_root=path_destination)
        self.resrc = Resources(number_node=self.resrc_para['number_node'],
                               cpu_per_node = self.resrc_para['cpu_per_node'],
                               gpu_per_node = self.resrc_para['gpu_per_node'],
                               group_size = self.resrc_para['group_size'],
                               queue_name = self.resrc_para['queue_name'],
                               envs = self.resrc_para['envs']
                               )
        task_list = []
        for item in os.listdir(path_source):
            isomer_name = os.path.splitext(item)[0]
            os.makedirs(os.path.join(path_destination, isomer_name), exist_ok=True)
            path_item = os.path.join(path_source, item)
            atoms = read(path_item)
            gjf_name = 'gau.gjf'
            gjf_path = os.path.join(path_destination, isomer_name, gjf_name)
            with open(gjf_path, 'w') as fd:
                write_gaussian_in(fd=fd, atoms=atoms, method=f'{self.job_type}(MaxCycles={cycles}) {self.method} {self.other_para}',
                                  charge=self.charge, mult=self.mult, mem=self.mem, nprocshared=self.resrc_para['cpu_per_node'], chk='gau.chk')
            cmd_list.append(gjf_name)
            a_task = Task(command=' '.join(cmd_list),
                          task_work_path=f"{isomer_name}/",
                          forward_files=[gjf_name],
                          backward_files=self.out_list
                          )
            task_list.append(a_task)
            del cmd_list[-1]

        if self.sub_batch_size == None:
            submission = Submission(work_base=path_destination,
                                    machine=self.mach,
                                    resources=self.resrc,
                                    task_list=task_list,
                                    forward_common_files=[],
                                    backward_common_files=[]
                                    )
            try:
                submission.run_submission()
            except Exception as e:
                with open('opt_error.log', 'a') as f:
                    f.write(str(e) + '\n' + str(task_list) + '\n')

        else:
            num_groups = math.ceil(len(task_list) / self.sub_batch_size)
            for i in range(num_groups):
                cursor = i * self.sub_batch_size
                a_task_list = task_list[cursor:cursor + self.sub_batch_size]
                submission = Submission(work_base=path_destination,
                                        machine=self.mach,
                                        resources=self.resrc,
                                        task_list=a_task_list,
                                        forward_common_files=[],
                                        backward_common_files=[]
                                        )
                try:
                    submission.run_submission()
                except Exception as e:
                    with open('opt_error.log', 'a') as f:
                        f.write(str(e) + '\n' + str(task_list) + '\n')


    def opt_once(self):
        if self.cage.add_num % 2 == 1:
            self.mult = 2
        else:
            self.mult = 1
        self.run_a_batch(path_source=self.path_raw_init, path_destination=self.path_opt_init, cmd_list=self.cmd_list.copy(), cycles=self.init_cycle)
        self.checker.check(opt_mood=self.mode, opt_root=self.path_opt_init, is_init=True, init_cycle=self.init_cycle)
        if os.stat('wrong_paths').st_size != 0:
            if self.deal_wrong_mode == 'report':
                print(f'Something wrong happened while optimizing isomers in {self.path_raw_init}.')
                return 0
            elif self.deal_wrong_mode == 'tough':
                return 2

        # if all jobs are failed, stop opt
        if (os.stat('yes_paths').st_size != 0) or (os.stat('init_yes_paths').st_size != 0):
            return 1
        else:
            return 0

