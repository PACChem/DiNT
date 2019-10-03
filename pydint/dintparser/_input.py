"""
parses the input file for keywords
"""

from autoparse.find import first_capture
from autoparse.find import all_captures
from autoparse.pattern import capturing
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_or_more
from autoparse.pattern import escape
from autoparse.pattern import NONSPACE
from autoparse.pattern import SPACE
from autoparse.pattern import WILDCARD
from autoparse.pattern import INTEGER
from autoparse.pattern import LINE_FILL
from autoparse.pattern import NONNEWLINE
from autoparse.pattern import NEWLINE


INPUT_SUPPORTED_SECTIONS = [
    'trajectory',
    'collision',
    'atom_groups'
]
INPUT_REQUIRED_SECTIONS = [
    'trajectory',
    'collision',
    'atom_groups'
]

DINT_SUPPORTED_KEYWORDS = [
    'pot_flag',					#Record 1
    'nsurf0', 'nsurft', 'methflag', 'repflag',	#Record 2
    'intflag',					#Record 3
    'hstep', 'eps', 'nprint',			#Record 3.0/1
    'ranseed',					#Record 4
    'ntraj', 'tflag1', 'tflag2', 'tflag3', 'tflag4',#Record 5
#    'ramptime', 'rampfact', 'nramp',		#Record 5.1
#    'andersen_temp', 'andersen_freq', 'scandth',#Record 5.2
#    'trajlist',					#Record 5.3
#    'ntarget', 'ephoton', 'wphoton',		#Record 5.4
    'nmol', 'ezero',				#Record 6
    'natom', 'initx', 'initp', 'initj', 'ezero_i',#Record 7
    'target_mass_xyz',				#Record 8.0
#    'rdum',					#Record 8.1.1
    'target_mass',				#Record 8.1.2
#    'lreadhess', 'nmtype', 'nmqn',		#Record 8.2.1
    'escatad', 'vvad', 'jjad', 'rrad', 'arrad',	#Record 8.3.2
    'temp0im', 'scale0im',			#Record 8.5.1
#    'samptot', 'lbinsamp', 'sampfilexx', 'sampfilepp',#Record 8.6.1
    'escale0im',				#Record 9.0
    'samptarg', 'letot', 'sampjmin', 'sampjmax', 
    'sampjtemp1', 'sampjtemp2', 'sampbrot1', 'sampbrot2',#Record 10.1
    'ejsc',					#Record 10.1.2
    'iorient', 'ldofrag',			#Record 11
#    'xx', 'pp',					#Record 11.0.1/2
    'rel0qc', 'ttt', 'bminqc', 'bmaxqc',	#Record 11.1
    'termflag', 'tnstep',			#Record 12
#    'tstime',					#Record 12.1
#    'tgradmag',					#Record 12.2
    'tnoutcome', 'tsymb1', 'tsymb2', 'distcut',	#Record 12.3
    'ioutput', 'ilist'				#Record 13
]
DINT_REQUIRED_KEYWORDS = [
    'pot_flag',					#Record 1
    'nsurf0', 'nsurft', 'methflag', 'repflag',	#Record 2
    'intflag',					#Record 3
    'hstep0', 'nprint',				#Record 3.0/1
#    'ranseed',					#Record 4
    'ntraj', 'tflag1', 'tflag2', 'tflag3', 'tflag4',#Record 5
    'nmol', 'ezero',				#Record 6
    'natom', 'initx', 'initp', 'initj', 'ezero_i',#Record 7
#    'iorient', 'ldofrag',			#Record 11
    'termflag', 'tnstep',			#Record 12
    'ioutput'					#Record 13
]


#----------------------------------------------------------------------------#
# Read the keywords from the trajectory setup section

def read_potential(input_string):
    """ obtain the potential to be used """
    pattern = ('potential' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(one_or_more(NONSPACE)))
    block = _get_potential_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword == 'tbplusexp6all'

    return keyword

def read_pot_flag(input_string):
    """ obtain the potential type to be used """
    pattern = ('pot_flag' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_potential_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_nsurf(input_string,input_ns): #CHECK
    """ obtain the electronic surfaces to be used """
    pattern = ('nsurf' + str(input_ns) +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_method(input_string):
    """ obtain the nonadiabatic method flag to be used """
    pattern = ('methflag' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_rep(input_string):
    """ obtain the adiabatic/diabatic representation flag to be used """
    pattern = ('repflag' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_integrator(input_string):
    """ obtain the integrator flag to be used """
    pattern = ('intflag' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_hstep(input_string):
    """ obtain the initial time step (in fs) to be used """
    pattern = ('hstep' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_nprint(input_string):
    """ obtain the frequency of printing every nprint steps to be used """
    pattern = ('nprint' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_ranseed(input_string):
    """ obtain the ranseed for SPRNG to be used """
    pattern = ('ranseed' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_ntraj(input_string):
    """ obtain the number of trajectories to be run per input file """
    pattern = ('ntraj' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_tflags(input_string,input_integer):
    """ obtain the trajectory flag options to be used """
    pattern = ('tflag' + str(nput_integer) +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_trajectory_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def _get_trajectory_options_section(input_string):
    """ grabs the section of text containing all of traj keywords """
    pattern = (escape('$trajectory') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


#----------------------------------------------------------------------------#
# Read the atom groups sections

def read_nmol(input_string):
    """ obtain the number of molecules (atom groups) in trajectory """
    pattern = ('nmol' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_ags_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_ezero(input_string):
    """ obtain the total minimum of energy (in eV) of all atom groups"""
    pattern = ('ezero' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_ags_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_ags(input_string,ag_number):
    """ builds a dictionary containing all needed info for the atom groups """
    ags_section = _get_ags_section(input_string)
    atom_group = _get_atom_group(ags_section,ag_number)

    ags_dct = {}
    for line in atom_group.splitlines():
        tmp = line.strip().split()
        assert len(tmp) >= 5
        name, mass, x, y, z = tmp[0], tmp[1], tmp[2], tmp[3], tmp[4]
        ags_dct[name] = [float(mass), float(x), float(y), float(z)]

    assert ags_dct

    return ags_dct


def read_natom_ag(input_string,ag_number):
    """ obtain the number of atoms in specified atom group """
    ags_section = _get_ags_section(input_string)
    atom_group = _get_atom_group(ags_section,ag_number)

    pattern = ('natom' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))

    keyword = first_capture(pattern, atom_group)

    assert keyword is not None
    keyword = int(keyword)

    return keyword

def read_initx_ag(input_string,ag_number):
    """ obtain the initx flag for specified atom group """
    ags_section = _get_ags_section(input_string)
    atom_group = _get_atom_group(ags_section,ag_number)

    pattern = ('initx' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))

    keyword = first_capture(pattern, atom_group)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_initp_ag(input_string,ag_number):
    """ obtain the initp flag for specified atom group """
    ags_section = _get_ags_section(input_string)
    atom_group = _get_atom_group(ags_section,ag_number)

    pattern = ('initp' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))

    keyword = first_capture(pattern, atom_group)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_initj_ag(input_string,ag_number):
    """ obtain the initj flag for specified atom group """
    ags_section = _get_ags_section(input_string)
    atom_group = _get_atom_group(ags_section,ag_number)

    pattern = ('initj' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))

    keyword = first_capture(pattern, atom_group)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_ezero_ag(input_string,ag_number):
    """ obtain the minimum of energy (in eV) of specified atom group """
    ags_section = _get_ags_section(input_string)
    atom_group = _get_atom_group(ags_section,ag_number)

    pattern = ('ezero_i' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))

    keyword = first_capture(pattern, atom_group)

    assert keyword is not None
    keyword = float(keyword)

    return keyword


def _get_ags_section(input_string):
    """ grabs the section of text containing all of the atom groups """
    pattern = (escape('$atom_groups') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


def _get_atom_group(input_string,input_integer):
    """ grabs the section of text containing all of the atom groups """
    pattern = (escape('$AG') + input_integer +
               + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$endAG'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


#----------------------------------------------------------------------------#
# Read the collision parameter section

def read_termflag(input_string):
    """ obtain the termination flag option """
    pattern = ('termflag' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_collision_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_tnstep(input_string):
    """ obtain the number of steps before termination """
    pattern = ('tnstep' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_collision_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_ioutput(input_string):
    """ obtain the ioutput flag option """
    pattern = ('ioutput' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_collision_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def read_ilist(input_string):
    """ obtain the list of output files option """
    pattern = ('ilist' +
               zero_or_more(SPACE) + '=' + zero_or_more(SPACE) +
               capturing(INTEGER))
    block = _get_collision_options_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    keyword = int(keyword)

    return keyword


def _get_collision_options_section(input_string):
    """ grabs the section of text containing all collision keywords """
    pattern = (escape('$collision') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


#----------------------------------------------------------------------------#
# Functions to check for errors in the input file

def check_defined_sections(input_string):
    """ verify all defined sections have been defined """
    pattern = (escape('$') + capturing(one_or_more(NONNEWLINE)))

    matches = all_captures(pattern, input_string)

    # See if each section has an paired end and is a supported keywords
    defined_sections = []
    for i, match in enumerate(matches):
        if (i+1) % 2 == 0:
            if match != 'end':
                raise ValueError
        else:
            defined_sections.append(match)

    # Check if sections are supported
    if not all(section in INPUT_SUPPORTED_SECTIONS
               for section in defined_sections):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(section in defined_sections
               for section in INPUT_REQUIRED_SECTIONS):
        raise NotImplementedError


def check_defined_dint_keywords(input_string):
    """ obtains the keywords defined in the input by the user """
    section_string = _get_dint_options_section(input_string)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in DINT_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in DINT_REQUIRED_KEYWORDS):
        raise NotImplementedError


def _get_defined_keywords(section_string):
    """ gets a list of all the keywords defined in a section """
    defined_keys = []
    for line in section_string.splitlines():
        if '=' in line:
            tmp = line.strip().split('=')[0]
            defined_keys.append(tmp.strip())

    return defined_keys
