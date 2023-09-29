# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 10:23:38 2021

@author: Philipp
"""

import pyteomics.mass 
import pyteomics.electrochem


class ProteomicTools:
    def __init__(self):
        pass 
    
    def calculate_mass(self, sequence):
        return pyteomics.mass.calculate_mass(sequence)
    
    def calculate_gravy(self, sequence):
        if "U" in sequence: return 0
        try:
            return pyteomics.electrochem.gravy(sequence)
        except:
            return 0
    
    def calculate_pI(self, sequence):
        try:
            return pyteomics.electrochem.pI(sequence)
        except:
            return 0