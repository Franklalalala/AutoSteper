import os

import networkx as nx
import numpy as np
import pandas as pd
from ase.io import read
from autosteper.tools import get_low_e_ranks
from fullerenedatapraser.molecular.fullerene import FullereneFamily


class Cage():
    def __init__(self, pristine_path):
        self.pristine_path = os.path.abspath(pristine_path)
        self.name, fmt = os.path.splitext(os.path.basename(pristine_path))
        # parameter for standard name
        self.atoms = read(pristine_path, format=fmt[1:])
        symbol_set = self.atoms.symbols.species()
        assert len(symbol_set) == 1, 'Currently only support pure cages that contain only one element.'
        self.symbol = self.atoms[0].symbol
        self.centre = self.atoms.get_positions().mean(axis=0)
        self.size = len(self.atoms)
        self.max_add_36_size = len(to_36_base(int('1' * self.size, 2)))
        # graph6str
        self._get_graph6str()
        self.blk_list = None
        self.has_blk_list = False
        self.is_pre_scan = False
        self.ps_num_list = []
        self.calc = None

    def _get_graph6str(self):
        f = FullereneFamily(spiral=0, atoms=self.atoms)
        self.graph6str = nx.to_graph6_bytes(f.graph, header=False).decode().split()[0]
        if os.name == 'posix':
            self.graph6str = self.graph6str.replace('`', '\\`')

    def set_workbase(self, root: str):
        self.workbase = os.path.abspath(os.path.join(root, self.name))
        os.makedirs(self.workbase, exist_ok=True)

    def set_add_num(self, add_num: int = None):
        self.add_num = add_num
        self.addon_path = os.path.join(self.workbase, f'{self.add_num}addons')
        os.makedirs(name=self.addon_path, exist_ok=True)
        os.chdir(self.addon_path)

        # blk is relatively late to apply (one step)
        if self.blk_list:
            if add_num in self.blk_list.num_list:
                self.has_blk_list = True
            else:
                self.has_blk_list = False
        if add_num in self.ps_num_list:
            self.is_pre_scan = True
        else:
            self.is_pre_scan = False


def to_36_base(num):
    return ((num == 0) and "0") or (
                to_36_base(num // 36).lstrip("0") + "0123456789abcdefghijklmnopqrstuvwxyz"[num % 36])


def name2seq(name: str, cage_size: int):
    seq = str()
    addon_set = set()
    bin_str = format(int(name, 36), 'b')
    bin_str = '0' * (cage_size - len(bin_str)) + bin_str
    bin_list = []
    for idx, a_position in enumerate(bin_str):
        if a_position == '1':
            seq += str(idx) + ' '
            addon_set.add(idx)
            bin_list.append(1)
        else:
            bin_list.append(0)
    return seq, addon_set, np.array(bin_list)


def seq2name(seq: str, cage: Cage):
    addon_list = []
    bin_name = ['0'] * cage.size
    bin_arr = np.zeros(cage.size)
    for i in seq.split():
        int_i = int(i)
        bin_name[int_i] = '1'
        addon_list.append(int_i)
        bin_arr[int_i] = 1
    bin_name_str = ''.join(bin_name)
    base_36_name = to_36_base(int(bin_name_str, 2))
    for i in range(cage.max_add_36_size - len(base_36_name)):
        base_36_name = '0' + base_36_name
    addon_set = set(addon_list)
    return base_36_name, addon_set, bin_arr


class blk_list():
    def __init__(self, size: int, blk_para: dict = None):
        self.num_list = range(blk_para['start_clct_num'], blk_para['final_chk_num'] + 1)
        self.clct_unstb = blk_para['clct_unstb']
        self.failed_arr = np.ones(size)
        self.blk_list_arr = np.ones(size)
        if self.clct_unstb:
            self.unstb_para = blk_para['unstb_para']
            self.collector = None
            self.container = []
            self.container_size = blk_para['container_size']
            self.cage_size = size

    def clct_failed(self, status_info_path: str):
        status_info = pd.read_pickle(status_info_path)
        status = status_info.T[0]
        names = status_info.columns
        for idx, a_status in enumerate(status):
            if a_status == 0:
                continue
            a_name = names[idx]
            _, _0, a_failed_arr = name2seq(name=a_name, cage_size=self.cage_size)
            self.failed_arr = np.vstack((self.failed_arr, a_failed_arr))

    def clct_unstable(self, info_path: str):
        self.collector = np.ones(self.cage_size)
        deep_yes = pd.read_pickle(info_path)
        passed_info = pd.read_pickle(info_path)
        e_arr = np.array(passed_info['energy'])
        names = passed_info['name']
        for a_rank in get_low_e_ranks(e_arr=e_arr, para=self.unstb_para, is_reverse=True):
            a_name = names[a_rank]
            _, _0, bin_arr = name2seq(name=a_name, cage_size=self.cage_size)
            self.collector = np.vstack((self.collector, bin_arr))

        self.container.append(self.collector.copy())
        if len(self.container) > self.container_size:
            del self.container[0]
        self.blk_list_arr = self.failed_arr
        for an_unstable_arr in self.container:
            self.blk_list_arr = np.vstack((self.blk_list_arr, an_unstable_arr))

    def chk_blk(self, q_bin_arr_list: list):
        blk_bin_arr_T = self.blk_list_arr.T
        res_arrs = np.array(q_bin_arr_list) @ blk_bin_arr_T - blk_bin_arr_T.sum(axis=0)
        uncutted_idx_list = []
        for an_idx in range(len(q_bin_arr_list)):
            if not 0 in res_arrs[an_idx]:
                uncutted_idx_list.append(an_idx)
        return uncutted_idx_list
