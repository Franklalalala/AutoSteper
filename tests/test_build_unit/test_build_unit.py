import os

from ase.io import write
from autosteper.cage import Cage
from autosteper.generator import Generator


pristine_cage = Cage(pristine_path=r'geom.xyz')
# ==============================================================================
test_Gen = Generator(gen_para={
    'group': 'OH',
    'geom_mode': 'pre_defined',},
    pst_cage=pristine_cage)

test_Gen.pre_build_unit(prev_xyz_path=r'C60Cl2.xyz')
new_cage = test_Gen.build_unit(new_addon_sites=[37, 44])
write(filename=r'./output/C60Cl2(OH)2.xyz', images=new_cage, format='xyz')

# ==============================================================================
test_Gen = Generator(gen_para={
    'group': 'CF3',
    'geom_mode': 'pre_defined',},
    pst_cage=pristine_cage)

test_Gen.pre_build_unit(prev_xyz_path=r'C60Cl2.xyz')
new_cage = test_Gen.build_unit(new_addon_sites=[37, 44])
write(filename=r'./output/C60Cl2(CF3)2.xyz', images=new_cage, format='xyz')

# ==============================================================================
test_Gen = Generator(gen_para={
    'group': 'CH3',
    'geom_mode': 'pre_defined',},
    pst_cage=pristine_cage)

test_Gen.pre_build_unit(prev_xyz_path=r'C60Cl2.xyz')
new_cage = test_Gen.build_unit(new_addon_sites=[37, 44])
write(filename=r'./output/C60Cl2(CH3)2.xyz', images=new_cage, format='xyz')

# ==============================================================================
test_Gen = Generator(gen_para={
    'group': 'F',
    'geom_mode': 'pre_defined',},
    pst_cage=pristine_cage)

test_Gen.pre_build_unit(prev_xyz_path=r'C60Cl2.xyz')
new_cage = test_Gen.build_unit(new_addon_sites=[37, 44])
write(filename=r'./output/C60Cl2F2.xyz', images=new_cage, format='xyz')

# ==============================================================================
test_Gen = Generator(gen_para={
    'group': 'Br',
    'geom_mode': 'pre_defined',},
    pst_cage=pristine_cage)

test_Gen.pre_build_unit(prev_xyz_path=r'C60Cl2.xyz')
new_cage = test_Gen.build_unit(new_addon_sites=[37, 44])
write(filename=r'./output/C60Cl2Br2.xyz', images=new_cage, format='xyz')

# ==============================================================================
test_Gen = Generator(gen_para={
    'group': 'Cl',
    'geom_mode': 'pre_defined',},
    pst_cage=pristine_cage)

test_Gen.pre_build_unit(prev_xyz_path=r'C60Cl2.xyz')
new_cage = test_Gen.build_unit(new_addon_sites=[37, 44])
write(filename=r'./output/C60Cl4.xyz', images=new_cage, format='xyz')
