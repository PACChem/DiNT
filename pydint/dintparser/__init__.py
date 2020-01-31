"""
Functions to parse the input
"""

# Functions for trajectory section
from pydint.dintparser._input import read_potential
from pydint.dintparser._input import read_pot_flag
from pydint.dintparser._input import read_nsurf
from pydint.dintparser._input import read_method
from pydint.dintparser._input import read_rep
from pydint.dintparser._input import read_integrator
from pydint.dintparser._input import read_hstep
from pydint.dintparser._input import read_nprint
from pydint.dintparser._input import read_ranseed
from pydint.dintparser._input import read_ntraj
from pydint.dintparser._input import read_tflags
# Functions for atom group sections
from pydint.dintparser._input import read_nmol
from pydint.dintparser._input import read_ezero
from pydint.dintparser._input import read_ags
from pydint.dintparser._input import read_natom_ag
from pydint.dintparser._input import read_initx_ag
from pydint.dintparser._input import read_initp_ag
from pydint.dintparser._input import read_initj_ag
from pydint.dintparser._input import read_ezero_ag
# Functions for collision section
from pydint.dintparser._input import read_termflag
from pydint.dintparser._input import read_tnstep
from pydint.dintparser._input import read_ioutput
from pydint.dintparser._input import read_ilist
# Functions to check for errors in the input file
from pydint.dintparser._input import check_defined_sections
from pydint.dintparser._input import check_defined_dint_keywords


__all__ = [
    'read_potential',
    'read_pot_flag',
    'read_nsurf',
    'read_method',
    'read_rep',
    'read_integrator',
    'read_hstep',
    'read_nprint',
    'read_ranseed',
    'read_ntraj',
    'read_tflags',
   
    'read_nmol',
    'read_ezero',
    'read_ags',
    'read_natom_ag',
    'read_initx_ag',
    'read_initp_ag',
    'read_initj_ag',
    'read_ezero_ag',
   
    'read_termflag',
    'read_tnstep',
    'read_ioutput',
    'read_ilist',
    'check_defined_sections',
    'check_defined_dint_keywords'
]
