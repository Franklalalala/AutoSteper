# This file is trying to visualize isomers, pathways, and SWR pairs in terms of Schlegel diagram with assistance of FullereneDataParser.
# FullereneDataParser(https://github.com/XJTU-ICP/FullereneDataParser) is an excellent tool designed and maintained by Yanbo Han, a XJTU-ICP member.
# A vivid description of FullereneDataParser could be found in:
# XCSI: "Effect of orbital angles on the modeling of the conjugate system with curvature". DOI:10.1039/D2CP03549A
# If any function or class in this file is utilized in your work, please cite the article above.
import os
import shutil
from typing import Union

import matplotlib.pyplot as plt
from ase.atoms import Atoms
from ase.io import read, write
from autosteper.tools import re_label_a_pathway, strip_extraFullerene, get_max_adduct_filename, match_max_adduct, \
    match_swr_cages
from fullerenedataparser.graph.visualize.cage import planarity_graph_draw
from fullerenedataparser.molecular.fullerene import FullereneCage


def get_proj_point(ring_seq: list, circle_vertex_list):
    new_ring_seq = set()
    for a_vertex in ring_seq:
        new_ring_seq.add(a_vertex - 1)
    proj_point = 0
    if len(new_ring_seq) == 5:
        proj_point = proj_point - 1
        print('We strongly recommend to project from a hexagon.')
        print('Please set the pentage_alpha to 0 to get a pleasant view.')
        for a_circle in circle_vertex_list:
            if len(a_circle) != 5:
                continue
            a_circle_set = set(a_circle)
            if a_circle_set - new_ring_seq == set():
                break
            proj_point = proj_point - 1

    if len(new_ring_seq) == 6:
        for a_circle in circle_vertex_list:
            if len(a_circle) != 6:
                continue
            a_circle_set = set(a_circle)
            if a_circle_set - new_ring_seq == set():
                break
            proj_point = proj_point + 1
    return proj_point


color_map = {
    'Cl': 'c',
    'F': 'tan',
    'Br': 'gold',
    'I': 'm',
    'H': 'deepskyblue',
    'OH': 'red',
    'CH3': 'y',
    'CF3': 'g',
}


class FullereneDataParser_Plotter():
    """
    Plotter for Schlegel diagrams, powered by FullereneDataParser.
    FullereneDataParser(https://github.com/XJTU-ICP/FullereneDataParser) is an excellent tool
    designed and maintained by Yanbo Han, a XJTU-ICP member.
    An vivid description of FullereneDataParser could be found in:
    XCSI: "Effect of orbital angles on the modeling of the conjugate system with curvature". DOI:10.1039/D2CP03549A
    If this class is utilized in your work, please cite the article above.

    """

    def __init__(self):
        print(r'Note: if autosteper.plotter.plotter is utilized in your work, '
              r'please cite this article: DOI:10.1039/D2CP03549A')

    def load_cage(self, coord_file_path: str = None, atoms: Atoms = None,
                  show_C_label=False, C_label_color="black", C_label_transparency=1,
                  proj_ring_seq: Union[list, set] = None,
                  pentagon_color="orange", pentagon_transparency=0.5,
                  sphere_ratio=0.8, parr_ratio=0.4, ax=None):
        """

        Args:
            coord_file_path: the file that contains coordinates
            atoms: an ASE Atoms object
            show_C_label: set true to see Carbon labels
            C_label_color: color for carbon labels
            C_label_transparency: transparency of carbon labels
            proj_ring_seq: which ring to project from
            pentagon_color: color of pentagon fill
            pentagon_transparency: the transparency of pentagon fill color
            sphere_ratio, parr_ratio: ratio to control graph deformation between projection of platform and hemi-sphere.
            ax: plotter handle, default is none to create a new handle

        """
        if coord_file_path:
            atoms = read(coord_file_path)
        elif not atoms:
            raise RuntimeError('Please input geometry-containing file, or input an ASE Atoms object.')
        f = FullereneCage(spiral=0, atoms=atoms, nospiralflag=True)
        if proj_ring_seq:
            circle_vertex_list = f.circle_vertex_list
            projection_point = get_proj_point(proj_ring_seq, circle_vertex_list)
        else:
            projection_point = 0
        self.ax, self.pos = planarity_graph_draw(f,
                                                 sphere_ratio=sphere_ratio,
                                                 parr_ratio=parr_ratio,
                                                 atom_label=show_C_label,
                                                 line_color=C_label_color,
                                                 line_alpha=C_label_transparency,
                                                 projection_point=projection_point,
                                                 path=None,
                                                 pentage_alpha=pentagon_transparency,
                                                 pentage_color=pentagon_color,
                                                 ax=ax
                                                 )

    def load_addons(self, addon_set: Union[list, set], replace_addon_map: dict=None, addon_color: str = None, addon_label_size: int = 200,
                    group_symbol: str = None, show_addon_nums: bool = False, addon_nums_rotation: int = 0, fontsize: int=10):
        """

        Args:
            addon_set: cage sites that are functionalied by groups. (Caution: addon set start from 0)
            addon_color: color of addons on diagram
            addon_label_size: size of labels
            group_symbol: symbol of groups
            show_addon_nums: set true to see numbers of addons
            addon_nums_rotation: set true to rotate these numbers
            replace_addon_map: a map to replace addon numbers

        """

        if not group_symbol:
            group_symbol = 'H'
        if not addon_color:
            addon_color = color_map[group_symbol]
        for an_addon in addon_set:
            self.ax.scatter(self.pos[:, 0][an_addon], self.pos[:, 1][an_addon], c=addon_color, alpha=1,
                            s=addon_label_size)
            if show_addon_nums:
                if replace_addon_map:
                    self.ax.text(self.pos[:, 0][an_addon], self.pos[:, 1][an_addon], str(replace_addon_map[an_addon + 1]), fontsize=fontsize,
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 rotation=addon_nums_rotation)
                else:
                    self.ax.text(self.pos[:, 0][an_addon], self.pos[:, 1][an_addon], str(an_addon + 1), fontsize=fontsize,
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 rotation=addon_nums_rotation)

    def plot_pathway_unit(self, src_pathway_root: str, new_pathway_workbase: str, is_re_label: bool = True,
                          is_match_max_adduct: bool = False, max_adduct_path: str = None, diff_len: int = None,
                          proj_ring_seq: Union[list, set] = None,
                          show_C_label=False, C_label_color="black", C_label_transparency=1,
                          pentagon_color="orange", pentagon_transparency=0.5,
                          sphere_ratio=0.8, parr_ratio=0.4, addon_color: str = None, addon_label_size: int = 200,
                          group_symbol: str = None, show_addon_nums: bool = False, addon_nums_rotation: int = 0,
                          dpi=400, fontsize: int=10):
        """

        Args:
            src_pathway_root: the original pathway root
            new_pathway_workbase: new pathway workbase
            is_match_max_adduct: set true to match max adduct to a specific isomer
            max_adduct_path: the path to the specific isomer
            diff_len: how much shifts between two isomers
            is_re_label: set true to plot after re-label
            dpi: the quality of dumped pictures

            others: see above

        """

        cwd_plot_unit = os.getcwd()
        src_pathway_root = os.path.abspath(src_pathway_root)
        pathway_id = os.path.basename(src_pathway_root)
        new_pathway_workbase = os.path.abspath(new_pathway_workbase)
        os.chdir(new_pathway_workbase)

        if is_match_max_adduct:
            max_adduct_path = os.path.abspath(max_adduct_path)
            pristine_cage, _ = strip_extraFullerene(coord_file_path=max_adduct_path, group=group_symbol)
            new_src_pathway_root = os.path.join(new_pathway_workbase, 'matched_src_pathway', pathway_id)
            shutil.copytree(src=src_pathway_root, dst=new_src_pathway_root)
            src_pathway_root = new_src_pathway_root
            os.chdir(src_pathway_root)

            max_adduct_filename = os.path.basename(max_adduct_path)
            if not max_adduct_filename in os.listdir('./'):
                old_max_adduct_filename = get_max_adduct_filename(r'./')
                half_matched_atoms, re_matched_atoms = match_max_adduct(query_atoms_path=old_max_adduct_filename,
                                                                        tgt_atoms_path=max_adduct_path,
                                                                        diff_len=diff_len,
                                                                        group_symbol=group_symbol)
                write(filename=old_max_adduct_filename, images=re_matched_atoms, format='xyz')

        os.chdir(new_pathway_workbase)
        if is_re_label:
            os.makedirs(name='re_label', exist_ok=True)
            os.chdir('re_label')
            os.makedirs(pathway_id, exist_ok=True)
            re_label_a_pathway(src_root=src_pathway_root, re_label_root=pathway_id, group_symbol=group_symbol,
                               has_header=True)
            src_pathway_root = os.path.abspath(pathway_id)
        os.chdir(new_pathway_workbase)
        os.makedirs(name='pics', exist_ok=True)
        os.chdir('pics')
        os.makedirs(pathway_id, exist_ok=True)
        os.chdir(pathway_id)

        for a_file in os.listdir(src_pathway_root):
            if not a_file.endswith('xyz'):
                continue
            if is_match_max_adduct:
                _, addon_set = strip_extraFullerene(coord_file_path=os.path.join(src_pathway_root, a_file),
                                                    group=group_symbol)
            else:
                pristine_cage, addon_set = strip_extraFullerene(coord_file_path=os.path.join(src_pathway_root, a_file),
                                                                group=group_symbol)
            self.load_cage(atoms=pristine_cage, proj_ring_seq=proj_ring_seq,
                           show_C_label=show_C_label, C_label_color=C_label_color,
                           C_label_transparency=C_label_transparency,
                           pentagon_transparency=pentagon_transparency, pentagon_color=pentagon_color,
                           sphere_ratio=sphere_ratio, parr_ratio=parr_ratio)

            self.load_addons(addon_set=addon_set, group_symbol=group_symbol,
                             addon_color=addon_color, addon_label_size=addon_label_size,
                             show_addon_nums=show_addon_nums, addon_nums_rotation=addon_nums_rotation, fontsize=fontsize)
            plt.savefig(os.path.splitext(a_file)[0] + '.png', dpi=dpi)
            plt.close()
        os.chdir(cwd_plot_unit)

    def plot_swr_unit(self, src_swr_root: str, proj_ring_seq: Union[list, set] = None,
                      show_C_label=False, C_label_color="black", C_label_transparency=1,
                      pentagon_color="orange", pentagon_transparency=0.5,
                      sphere_ratio=0.8, parr_ratio=0.4, addon_color: str = None, addon_label_size: int = 200,
                      group_symbol: str = None, show_addon_nums: bool = False, addon_nums_rotation: int = 0,
                      dpi=400, fontsize: int=10):
        """

        Args:
            src_swr_root: root to the AutoSteper generated SWR pairs
            others: see above

        """
        cwd_swr = os.getcwd()
        src_swr_root = os.path.abspath(src_swr_root)
        if len(os.listdir(src_swr_root)) == 0:
            print(f'There are no swr reported in {os.getcwd()}')
            return
        a_folder = os.listdir(src_swr_root)[0]
        os.chdir(src_swr_root)
        os.chdir(a_folder)
        for a_file in os.listdir('./'):
            if a_file.endswith('xyz') and a_file.startswith('tgt'):
                target_cage, _ = strip_extraFullerene(coord_file_path=a_file, group=group_symbol)
            if a_file.endswith('xyz') and a_file.startswith('q_'):
                query_cage, _ = strip_extraFullerene(coord_file_path=a_file, group=group_symbol)
        new_query_cage, match_map, swr_sites = match_swr_cages(query_cage=query_cage, target_cage=target_cage)

        for a_folder in os.listdir(src_swr_root):
            os.chdir(src_swr_root)
            os.chdir(a_folder)
            for a_file in os.listdir('./'):
                if not a_file.endswith('xyz'):
                    continue
                _, addon_set = strip_extraFullerene(coord_file_path=a_file)
                if a_file.startswith('tgt'):
                    pristine_cage = target_cage
                else:
                    pristine_cage = new_query_cage
                    new_addon_set = set()
                    for an_old_site in addon_set:
                        new_addon_set.add(match_map[an_old_site])
                    addon_set = new_addon_set
                self.load_cage(atoms=pristine_cage, proj_ring_seq=proj_ring_seq,
                               show_C_label=show_C_label, C_label_color=C_label_color,
                               C_label_transparency=C_label_transparency,
                               pentagon_transparency=pentagon_transparency, pentagon_color=pentagon_color,
                               sphere_ratio=sphere_ratio, parr_ratio=parr_ratio)

                self.load_addons(addon_set=addon_set, group_symbol=group_symbol,
                                 addon_color=addon_color, addon_label_size=addon_label_size,
                                 show_addon_nums=show_addon_nums, addon_nums_rotation=addon_nums_rotation, fontsize=fontsize)

                self.load_addons(addon_set=swr_sites, addon_color='black', addon_label_size=addon_label_size - 50,
                                 show_addon_nums=show_addon_nums, addon_nums_rotation=addon_nums_rotation, fontsize=fontsize)
                plt.savefig(os.path.splitext(a_file)[0] + '.png', dpi=dpi)
                plt.close()
        os.chdir(cwd_swr)

