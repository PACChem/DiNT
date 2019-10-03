"""
input file writing routines for parallel DiNT runs
"""

import os
import subprocess
import random
import moldr
import elstruct
import pydint.interface


def write_input(job_dir_path, nsamp,
                target_name='input.dat', bath_name='bath,xyz',
                smin=2, smax=6):
    """ write the input file
    """

    ranseed = random.randrange(1E8, 1E9)
    inp_str = py1dmin.interface.writer.onedmin_input(
        ranseed, nsamp, target_name, bath_name, smin, smax)

    job_file_path = os.path.join(job_dir_path, 'input.dat')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)


def write_inp_traj()

def write_inp_ags()

def write_inp_coll()

def combine_inp_files()

def write_xyz(job_dir_path, target_geo, bath_geo):
    """ write the target and bath xyz files
    """

    job_file_path = os.path.join(job_dir_path, 'target.xyz')
    with open(job_file_path, 'w') as xyz_file:
        xyz_file.write(target_geo)

    job_file_path = os.path.join(job_dir_path, 'bath.xyz')
    with open(job_file_path, 'w') as xyz_file:
        xyz_file.write(bath_geo)

