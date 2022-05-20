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
