import numpy as np
from numpy import sin, cos
import pandas as pd
import os
from ase.atoms import Atoms



def rotate_around_axis(norm_axis: np.array, old_vec: np.array, theta: float):
    a, b, c = norm_axis
    rotate_mat = np.array([[(b ** 2 + c ** 2) * cos(theta) + a ** 2, a * b * (1 - cos(theta)) - c * sin(theta),
                            a * c * (1 - cos(theta)) + b * sin(theta)],
                           [a * b * (1 - cos(theta)) + c * sin(theta), b ** 2 + (1 - b ** 2) * cos(theta),
                            b * c * (1 - cos(theta)) - a * sin(theta)],
                           [a * c * (1 - cos(theta)) - b * sin(theta), b * c * (1 - cos(theta)) + a * sin(theta),
                            c ** 2 + (1 - c ** 2) * cos(theta)]])
    new_vec = old_vec@rotate_mat
    return new_vec




def get_yes_info(opt_mode: str, all_parent_info: dict=None):
    cwd_ = os.getcwd()
    deep_yes_path = os.path.join(cwd_, 'deep_yes_info.pickle')
    if os.path.exists(deep_yes_path):
        return deep_yes_path
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
                if opt_mode in ['xtb', 'gaussian']:
                    with open(a_xyz_path, 'r') as xyz_f:
                        xyz_f.readline()
                        energy_line = xyz_f.readline()
                        energy = float(energy_line.split()[1])
                elif opt_mode == 'ase':
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
    return deep_yes_path





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
