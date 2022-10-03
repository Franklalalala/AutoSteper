import numpy as np
import pandas as pd
from ase.atoms import Atoms
from ase.neighborlist import build_neighbor_list
from numpy import sin, cos


def rotate_around_axis(norm_axis: np.array, old_vec: np.array, theta: float):
    a, b, c = norm_axis
    rotate_mat = np.array([[(b ** 2 + c ** 2) * cos(theta) + a ** 2, a * b * (1 - cos(theta)) - c * sin(theta),
                            a * c * (1 - cos(theta)) + b * sin(theta)],
                           [a * b * (1 - cos(theta)) + c * sin(theta), b ** 2 + (1 - b ** 2) * cos(theta),
                            b * c * (1 - cos(theta)) - a * sin(theta)],
                           [a * c * (1 - cos(theta)) - b * sin(theta), b * c * (1 - cos(theta)) + a * sin(theta),
                            c ** 2 + (1 - c ** 2) * cos(theta)]])
    new_vec = old_vec @ rotate_mat
    return new_vec


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


def sort_atomic(atoms: Atoms):
    cage = Atoms()
    addon = Atoms()
    for an_atom in atoms:
        if an_atom.symbol == 'C':
            cage.append(an_atom)
        else:
            addon.append(an_atom)
    pristine_cage = cage.copy()
    cage.extend(addon)
    return pristine_cage, cage


def strip_extraFullerene(atoms: Atoms, is_group: bool = None, group: str = None, cage_size: int = None):
    if is_group:
        pristine_cage, atoms = deal_group(atoms=atoms, group=group, cage_size=cage_size)
    else:
        pristine_cage, atoms = sort_atomic(atoms=atoms)
    if group in ['CH3', 'CF3']:
        cutoffs = [natural_cutoffs(pristine_cage)[0]] * len(atoms)
        neighborlist = NeighborList(cutoffs, skin=0.15, bothways=True, self_interaction=False)
    else:
        neighborlist = build_neighbor_list(atoms=atoms, skin=0.15, bothways=True, self_interaction=False)
    dok_matrix = neighborlist.get_connectivity_matrix()
    bonds = np.array([*dok_matrix.keys()])
    map_dict = dict()
    for i in range(bonds.shape[0]):
        bond = bonds[i]
        symbol_0 = atoms[bond[0]].symbol
        symbol_1 = atoms[bond[1]].symbol
        if not symbol_0 == symbol_1:
            if bond[0] < bond[1]:
                if symbol_0 == 'C':
                    cage_element_idx = bond[0]
                    addon_idx = bond[1]
                else:
                    cage_element_idx = bond[1]
                    addon_idx = bond[0]
                map_dict.update({addon_idx: cage_element_idx})
    addon_set = set(map_dict.values())
    return pristine_cage, addon_set


def get_low_e_ranks(e_arr: np.ndarray, para: dict, is_reverse: bool=False):
    assert para['mode'] in ['None', 'rank', 'value',
                            'value_and_rank', 'rank_and_value',
                            'value_or_rank', 'rank_or_value'], f'Please check your run cutoff mode keyword.'
    if para['mode'] == 'None':
        rank_list = range(len(e_arr))
    if is_reverse:
        if para['mode'] == 'rank':
            rank_list = range(len(e_arr))[::-1][:para['rank']]
        else:
            e_arr = e_arr - min(e_arr)
            max_e = max(e_arr)
            for idx, a_e in enumerate(e_arr[::-1]):
                if a_e <= max_e - para['value']:
                    break
            if para['mode'] == 'value_and_rank' or para['mode'] == 'rank_and_value':
                idx = range(min(idx, para['rank']))
            if para['mode'] == 'value_or_rank' or para['mode'] == 'rank_or_value':
                idx = range(max(idx, para['rank']))
            rank_list = range(len(e_arr))[::-1][:idx]
    else:
        if para['mode'] == 'rank':
            rank_list = range(len(e_arr))[:para['rank']]
        else:
            e_arr = e_arr - min(e_arr)
            for idx, a_e in enumerate(e_arr):
                if a_e > para['value']:
                    break
            if para['mode'] == 'value_and_rank' or para['mode'] == 'rank_and_value':
                idx = range(min(idx, para['rank']))
            if para['mode'] == 'value_or_rank' or para['mode'] == 'rank_or_value':
                idx = range(max(idx, para['rank']))
            rank_list = range(len(e_arr))[:idx]

    if 'nimg_th' in para.keys():
        new_rank_list = []
        for idx, a_nimg in enumerate(para['nimages']):
            if a_nimg >= para['nimg_th']:
                new_rank_list.append(rank_list[idx])
        rank_list = new_rank_list
    for a_rank in rank_list:
        yield a_rank
