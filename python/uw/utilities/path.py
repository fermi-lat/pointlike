""" Routines for file paths

    Author: Joshua Lande
"""
import os
import numpy as np

def expand(file):
    """ Simple routine to expand env. variables using gtlike convetion. """
    file = file.replace('(','{').replace(')','}')
    return os.path.expandvars(file)
