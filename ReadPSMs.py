# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 19:49:41 2021

@author: Phili

ReadPSMs exported PSMs as txt from ProteomeDiscoverer 2.2
"""


import pandas as pd

from collections import Counter

import MyTools



class ReadPSMs:
    def __init__(self, file_name, dataset_name = ""):
        # variables
        self.spectrum_file = "Spectrum File"
        self.name = "name"
        self.sequence = "Annotated Sequence"
        self.mass = "MH+ [Da]"
        self.charge = "Charge"
        self.retention_time = "RT [min]"
        self.xcorr = "XCorr"
        self.mz = "m/z [Da]"
        self.delta_m_ppm = "DeltaM [ppm]"
        self.isolation_interference = "Isolation Interference [%]"
        self.ion_injection_time = "Ion Inject Time [ms]"
        self.missed_cleavages = "# Missed Cleavages"
        self.master_protein = "Master Protein Accessions"
        
        self.all_features = [self.spectrum_file, self.sequence, self.mass, 
                             self.charge, self.retention_time, self.xcorr,
                             self.mz, self.delta_m_ppm, self.isolation_interference,
                             self.ion_injection_time, self.missed_cleavages, 
                             self.master_protein]
        self.dataset_name = dataset_name
        
        self.MyTools = MyTools.MyTools()
        # load data
        self.psms_df = self._load_psm_file(file_name)
    
    def _load_psm_file(self, file_name):
        psms_df = pd.read_csv(file_name ,sep="\t")
        # TODO
        # Exclude PSMs ohne Proteine Accession
        psms_df["name"] = list(psms_df[self.spectrum_file])
        psms_df["replicate"] = [None for i in psms_df.index]
        psms_df["condition"] = [None for i in psms_df.index]
        return psms_df
    
    
    def give_all_columns(self, psms_df):
        columns = list(psms_df.columns)
        return columns
    
    
    def give_all_digit_columns(self, psms_df):
        all_columns = self.give_all_columns(psms_df)
        digit_columns = []
        for column in all_columns:
            if self.MyTools.string_is_number(psms_df[column][0]):
                digit_columns.append(column)
        return digit_columns
    
    
    def give_file_name_list(self, psms_df):
        file_names = sorted(list(set(psms_df["name"])))
        return file_names
    
    def give_all_peptides(self):
        peptides = sorted(list(set(self.psms_df[self.sequence])))
        return peptides
    
    def set_dataset_name(self, name):
        self.name = name
    
    
    def replace_file_name(self, old_file_name, new_file_name):
        self.psms_df = self.psms_df.replace([old_file_name], 
                                            [new_file_name], regex=True)

    def set_replicate(self, raw_file, replicate_name):
        self.psms_df.loc[self.psms_df[self.spectrum_file] == raw_file, 
                         'replicate'] = replicate_name
    
    
    def set_condition(self, raw_file, condition_name):
        self.psms_df.loc[self.psms_df[self.spectrum_file] == raw_file, 
                         'condition'] = condition_name
        

    def set_name(self, raw_file, name):
        self.psms_df.loc[self.psms_df[self.spectrum_file] == raw_file, 
                         'name'] = name

    
    
    def translate_feature(self, text):
        if text == "Mass":
            return self.mass
        elif text == "Retention time":
            return self.retention_time
        elif text == "XCorr":
            return self.xcorr
        elif text == "m/z":
            return self.mz
        elif text == "Delta mass":
            return self.delta_m_ppm
        elif text == "Isolation intereference":
            return self.isolation_interference
        elif text == "Ion injection time":
            return self.ion_injection_time
        elif text == "#Missed Cleavage":
            return self.missed_cleavages
        elif text == "peptide":
            return self.sequence
        elif text == "protein":
            return self.master_protein
        else:
            return None 
    
