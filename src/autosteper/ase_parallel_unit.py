import os
import shutil
import math
from ase.io import read, write
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from ase.calculators.emt import EMT


def run_a_batch(steps: int, fmax: float, path_raw: str, path_opt: str):
    for a_file in os.listdir(path_raw):
        name = a_file.split('.')[0]
        sub_opt = os.path.join(path_opt, name)
        os.makedirs(sub_opt, exist_ok=True)
        path_a_file = os.path.join(path_raw, a_file)
        shutil.copy(path_a_file, os.path.join(sub_opt, a_file))
        os.chdir(sub_opt)
        atoms = read(a_file, format='xyz')
        atoms.calc = calculator
        opt = FIRE(atoms, logfile='opt.log')
        opt.run(steps=steps, fmax=fmax)
        write(filename='opt.xyz', images=atoms, format='xyz')

# Please check arguments: opt folder, opt steps, fmax, optimizer...
calculator = EMT()
path_raw = str(os.path.abspath('./xyz/'))
print(path_raw)
path_opt = str(os.path.abspath('../../opt_100'))
print(path_opt)
run_a_batch(steps=100, fmax=0.0005, path_raw=path_raw, path_opt=path_opt)


