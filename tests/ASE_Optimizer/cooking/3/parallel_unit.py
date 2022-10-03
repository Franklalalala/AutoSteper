import os
import shutil
import math
from ase.io import read, write
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from torchlightmolnet.caculator import torchCaculator
from torchlightmolnet.lightning.molnet import LightMolNet
from torchlightmolnet.dataset.atomref import refatoms_xTB, get_refatoms
from torchlightmolnet import Properties
import torch



def run_a_batch(steps: int, fmax: float, path_raw: str, path_opt: str):
    for a_file in os.listdir(path_raw).copy():
        os.chdir(cwd_)
        name = a_file.split('.')[0]
        sub_opt = os.path.join(path_opt, name)
        os.makedirs(sub_opt)
        path_a_file = os.path.join(path_raw, a_file)
        shutil.copy(path_a_file, os.path.join(sub_opt, a_file))
        os.chdir(sub_opt)
        atoms = read(a_file, format='xyz')
        atoms.calc = calculator
        opt = FIRE(atoms, logfile='opt.log')
        opt.run(steps=steps, fmax=fmax)
        write(filename='opt.xyz', images=atoms, format='xyz')



model_path = r'/home/mkliu/schnet_opt/paper_4_27/Cl_final.ckpt'
net = LightMolNet(atomref=get_refatoms(refatoms_xTB)[Properties.energy_U0])
state_dict = torch.load(model_path)
net.load_state_dict(state_dict["state_dict"])
calculator = torchCaculator(net=net)

cwd_ = os.getcwd()
os.makedirs(name=r'cooked', exist_ok=True)
path_opt = os.path.abspath(r'cooked')
run_a_batch(steps=100, fmax=0.0005, path_raw='xyz', path_opt=path_opt)
