# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 08:36:43 2022

@author: Phili
"""

import matplotlib.pyplot as plt


class Settings:
    ## General
    version = "1.0.2"
    date = "December, 2022"
    
    ## Default Paths
    save_path = "/Results"
    
    ### Text
    gravy = "GRAVY score"
    pi = "pI"
    fractions_feature_distribution = "Fractions Feature Distribution"
    all_feature_distribution = "All Feature Distribution"
    fractions_elution_peptide = "Fractions Elution Peptide"
    fractions_identifications = "Fractions Identifications"
    fractions_peptide_overlap = "Fractions Peptide Overlap"
    fractions_protein_overlap = "Fractions Protein Overlap"
    fractions_elution_peptides_top10 = "Fraction Elution Top 10 Peptides"
    peptides_in_number_fractions = "Peptide in #Fractions"
    peptide_elution_difference = "Peptide Elution Difference"
    all_feature_distribution_all_datasets = "All Feature Distribution all Datasets"
    all_identifications_overlap = "All Identifications Overlap"
    number_identifications = "Fractions number of Identifications"
    fractions_protein_feature_distribution = "Fractions Protein Feature Distribution"
    analyze_orthogonality = "Analyze Orthogonality"
    compare_all_datasets = "Compare all Datasets"

    
    selection_analysis = [fractions_feature_distribution, all_feature_distribution, 
                          fractions_elution_peptide, fractions_identifications, 
                          fractions_peptide_overlap, fractions_protein_overlap, 
                          fractions_elution_peptides_top10, peptides_in_number_fractions]
    
    ###
    fata_file_open_dialog = 'Select Fasta File'
    fasta_file_filter = 'Fasta file (*.fasta);; All Files (*)'
    fasta_file_initial_filter = 'Fasta file (*.fasta)'
    psm_file_open_dialog = 'Select PSM File'
    psm_file_filter = 'txt file (*.txt);; All Files (*)'
    psm_file_initial_filter = 'txt file (*.txt)'
    
    export_file_open_dialog = 'Save results as'
    export_file_filter = 'csv file (*.csv);; All Files (*)'
    export_file_initial_filter = 'csv file (*.csv)'


    ###
    pvalue_thresholds = [[1e-4, "***"], [1e-3, "**"], [1e-2, "*"], [1, "ns"]]
     
    ### Default Settings 
    heatmap_colormap = "Greys"
    heatmap_colormap_options = plt.colormaps()
    
    
    # 
    
            # possible_features = ["MH+ [Da]", "Charge", "RT [min]", 
            #                      "XCorr", "m/z [Da]", "DeltaM [ppm]", 
            #                      "Isolation Interference [%]",
            #                      "Ion Inject Time [ms]", "# Missed Cleavages",
            #                      "GRAVY score", "pI"]