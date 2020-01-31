"""
Executes the automation part of DiNT
"""

import os
import sys
import moldr
import autofile
import pydint.dintroutines
import pydint.dintparser


# Set paths to various working directories
DRIVE_PATH = os.getcwd()
GEOM_PATH = os.path.join(DRIVE_PATH, 'geoms')

# Read the job submit input file into a string
with open(os.path.join(DRIVE_PATH, 'submit.dat'), 'r') as subfile:
    SUBMIT_FILE = subfile.read()

pydint.dintparser.check

subfile.close()

# Read the DiNT input file into a string
with open(os.path.join(DRIVE_PATH, 'input.dat'), 'r') as infile:
    INPUT_STRING = infile.read()

# Check to see if the parameters are defined
pydint.dintparser.check_defined_sections(INPUT_STRING)
pydint.dintparser.check_defined_dint_keywords(INPUT_STRING)

NJOBS = pydint.dintparser.read_njobs(INPUT_STRING)

# Read the trajectory parameters
POTENTIAL = pydint.dintparser.read_potential(INPUT_STRING)
POTFLAG = pydint.dintparser.read_pot_flag(INPUT_STRING)
NSURF0 = pydint.dintparser.read_nsurf(INPUT_STRING,0)
NSURFT = pydint.dintparser.read_nsurf(INPUT_STRING,t)
METHFLAG = pydint.dintparser.read_method(INPUT_STRING)
REPFLAG = pydint.dintparser.read_rep(INPUT_STRING)
INTEGRATOR = pydint.dintparser.read_integrator(INPUT_STRING)
HSTEP = pydint.dintparser.read_hstep(INPUT_STRING)
NPRINT = pydint.dintparser.read_nprint(INPUT_STRING)
#RANSEED = pydint.dintparser.read_ranseed(INPUT_STRING)
NTRAJ = pydint.dintparser.read_ntraj(INPUT_STRING)
for i in range(1,5):
    TFLAG[i] = pydint.dintparser.read_tflags(INPUT_STRING,i)

# Read the atom groups from the input
NMOL = pydint.dintparser.read_nmol(INPUT_STRING)
EZERO = pydint.dintparser.read_ezero(INPUT_STRING)
AG_DCTS = []
NATOM = []
INITX = []
INITP = []
INITJ = []
EZERO_I = []
for i in range(NMOL):
    AG_DCTS[i] = pydint.dintparser.read_ags(INPUT_STRING,i+1)
    NATOM[i] = pydint.dintparser.read_natom_ag(INPUT_STRING,i+1)
    INITX[i] = pydint.dintparser.read_initx_ag(INPUT_STRING,i+1)
    INITP[i] = pydint.dintparser.read_initp_ag(INPUT_STRING,i+1)
    INITJ[i] = pydint.dintparser.read_initj_ag(INPUT_STRING,i+1)
    EZERO_I[i] = pydint.dintparser.read_ezero_ag(INPUT_STRING,i+1)

# Read the collision parameters
TERMINATION = pydint.dintparser.read_term_flag(INPUT_STRING)
TNSTEP = pydint.dintparser.read_tnstep(INPUT_STRING)
#ORIENTATION = pydint.dubtparser.read_orient_conds(INPUT_STRING)
IOUTPUT = pydint.dintparser.read_ioutput(INPUT_STRING)
ILIST = pydint.dintparser.read_ilist(INPUT_STRING)

# Gather info for each AG
# Write each new input file
make dint_inp_1
print('Writing Trajectory input...')
pydint.dintroutines.write_input()

make dint_inp_2 for each AG
for i in range(NMOL):
    print('Writing AG' + NMOL + ' input...')
    pydint.dintroutines.write_input()

make dint_inp_3
print('Writing Collission input...')
pydint.dintroutines.write_input()

combine dint_inp


# Run an instance of DiNT for each input file
pydint.dintroutines.submit_job()

# Loop over the species and launch all of the desired jobs
for target, target_info in TARGET_DCTS.items():

    # Set the new information
    RUN_CHG = target_info[1] + BATH_LST[1]
    RUN_MLT = max([target_info[2], BATH_LST[2]])
    BATH_INFO = BATH_LST

    # Write the params to the run file system
    FS_THEORY_INFO = [THEORY_INFO[1],
                      THEORY_INFO[2],
                      moldr.util.orbital_restriction(
                          target_info, THEORY_INFO)]
    tgt_run_fs = autofile.fs.species(RUN_PREFIX)
    tgt_run_fs.leaf.create(target_info)
    tgt_run_path = tgt_run_fs.leaf.path(target_info)
    etrans_run_fs = autofile.fs.energy_transfer(tgt_run_path)
    etrans_run_path = etrans_run_fs.leaf.path(FS_THEORY_INFO)
    etrans_run_fs.leaf.create(FS_THEORY_INFO)

    # Run an instance of 1DMin for each processor
    for i in range(NJOBS):

        # Build run directory
        job_dir_path = os.path.join(
            etrans_run_path, 'run{0}'.format(str(i+1)))
        os.mkdir(job_dir_path)
        print('\n\nWriting files to'+job_dir_path)

        # Write the 1DMin input file
        print('  Writing input files...')
        py1dmin.ljroutines.write_input(
            job_dir_path, NSAMPS,
            target_name='target.xyz', bath_name='bath.xyz',
            smin=SMIN, smax=SMAX)

        # Write the electronic structure sumbission script
        print('  Writing electronic structure submission script...')
        py1dmin.ljroutines.write_elstruct_sub(
            job_dir_path, DRIVE_PATH, RUN_PROG)

    # Submit the job
    print('\n\nRunning each OneDMin job...')
    py1dmin.ljroutines.submit_job(DRIVE_PATH, etrans_run_path, NJOBS)

    # Read the lj from all the outputs
    print('\n\nAll OneDMin jobs finished.')
    print('\nReading the Lennard-Jones parameters...')
    SIGMA, EPSILONS = py1dmin.ljroutines.lj_parameters(etrans_run_path)
    if SIGMA is None and EPSILON is None:
        print('\nNo Lennard-Jones parameters found.')
        print('\n\nExiting OneDMin...')
        sys.exit()
    else: 
        # Grab the geometries and zero-energies
        print('\nReading the Lennard-Jones Potential Well Geometries...')
        GEOMS = py1dmin.ljroutines.lj_well_geometries(
            etrans_run_path)
        print('\nReading the Zero-Energies...')
        ENES = py1dmin.ljroutines.zero_energies(
            etrans_run_path)
