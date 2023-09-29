# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 10:11:02 2021

@author: Philipp
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 21:32:12 2021

@author: Phili
"""

from ProteomicTools import ProteomicTools
from Bio import SeqIO



class ReadFastaFile:
    def __init__(self):
        self.ProteomicTools = ProteomicTools()
    
    
    def read_fasta_file(self, fasta_file: str) -> dict:
        """ """
        fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
        fasta_dict = {}
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            accession = name.split("|")[1]
            fasta_dict[accession] = sequence
        return fasta_dict

    def find_sequence(self, fasta_dict: dict, seq: str) -> str:
        acc_list = []
        for acc in fasta_dict:
            if seq in fasta_dict[acc]:
                acc_list.append(acc)
        return acc_list 
                
    
    def find_protein(self, fasta_dict: dict, accession: str) -> str:
        return fasta_dict[accession]
    
    def analyze_protein(self, fasta_dict: dict, accession: str):
        sequence = fasta_dict[accession]
        mass = self.ProteomicTools.calculate_mass(sequence)
        pI = self.ProteomicTools.calculate_pI(sequence)
        gravy = self.ProteomicTools.calculate_gravy(sequence)
        return sequence, mass, pI, gravy
        
    
    def analyze_all_proteins(self, fasta_dict: dict) -> dict:
        protein_properties = {} #Accession: sequence, mass, pI, GRAVY
        for accession in fasta_dict:
            sequence, mass, pI, gravy = self.analyze_protein(fasta_dict, accession)
            protein_properties[accession] = sequence, mass, pI, gravy
        return protein_properties
            




