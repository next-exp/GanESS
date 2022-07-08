"""
code: containers.py
description: namedtuples describing miscellaenous containers to
pass info around.

credits: see ic_authors_and_legal.rst in /doc

last revised: JJGC, 10-July-2017
"""

import sys

from collections import namedtuple

this_module = sys.modules[__name__]

def _add_namedtuple_in_this_module(name, attribute_names):
    new_nametuple = namedtuple(name, attribute_names)
    setattr(this_module, name, new_nametuple)

for name, attrs in (('SensorData'     , 'NPMT PMTWL '),):
    _add_namedtuple_in_this_module(name, attrs)

# Leave nothing but the namedtuple types in the namespace of this module
del name, namedtuple, sys, this_module, _add_namedtuple_in_this_module