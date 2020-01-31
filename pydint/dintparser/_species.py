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


DINT_SUPPORTED_SECTIONS = [
    'lennard_jones',
    'properties',
    'baths',
    'targets'
]
DINT_REQUIRED_SECTIONS = [
    'lennard_jones',
    'baths',
    'targets'
]

DINT_SUPPORTED_KEYWORDS = [
    'theory_level',
    'potential',
    'nsamps',
    'njobs',
    'smin',
    'smax',
    'run_prefix',
    'save_prefix',
    'confs'
]
DINT_REQUIRED_KEYWORDS = [
    'theory_level',
    'potential',
    'nsamps',
    'njobs',
    'run_prefix',
    'save_prefix',
    'confs'
]

DINT_SUPPORTED_KEYWORDS = [
    'inchi',
    'smiles',
    'geom',
    'charge',
    'spin',
    'method',
    'basis',
]

