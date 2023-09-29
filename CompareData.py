# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 22:47:36 2022

@author: Phili
"""


import pandas as pd
import seaborn as sns
from collections import Counter
import statistics

import matplotlib.pyplot as plt

import createGraphs
import MyTools


class CompareData:
    def __init__(self):

        
        self.createGraphs = createGraphs.createGraphs("Results")
        self.MyTools = MyTools.MyTools()
        
        
    def all_feature_distributions(self, dataset_1, dataset_2, feature):
        feature_list_1 = dataset_1.analyze_distribution_all_identifications( \
                                    dataset_1.ReadPSMs.psms_df, feature = feature, 
                                    show_plot = False)
        feature_list_2 = dataset_2.analyze_distribution_all_identifications( \
                                    dataset_2.ReadPSMs.psms_df, feature = feature, 
                                    show_plot = False)
        
        dataset_name_1 = dataset_1.ReadPSMs.dataset_name
        dataset_name_2 = dataset_2.ReadPSMs.dataset_name
        self.createGraphs.plot_comparison_distributions_histogram(feature_list_1, 
            feature_list_2, feature, labels = [dataset_name_1, dataset_name_2], 
            show_boxplot = True)
        self.createGraphs.plot_comparison_distributions_histogram(feature_list_1, 
            feature_list_2, feature, density = True, 
            labels = [dataset_name_1, dataset_name_2], show_boxplot = False)
        
        
        
    def fractions_identifications(self, dataset_1, dataset_2):
        identifications_df_1 = dataset_1.analyze_fractions_by_identifications( \
                                    dataset_1.ReadPSMs.psms_df, False)
        identifications_df_2 = dataset_2.analyze_fractions_by_identifications( \
                                    dataset_2.ReadPSMs.psms_df, False)
        if list(identifications_df_1.index) != list(identifications_df_2.index):
            print("Datasets not comparable. Please change the names of the " 
                  "individual files.")
            return False
        dataset_name_1 = dataset_1.ReadPSMs.dataset_name
        dataset_name_2 = dataset_2.ReadPSMs.dataset_name
        for ident in list(identifications_df_1.columns):
            self.createGraphs.plot_comparison_fractions_identifications( 
                identifications_df_1[ident],  identifications_df_2[ident],
                list(identifications_df_1.index), ident, 
                labels = [dataset_name_1, dataset_name_2])
            
    
    def _combine_datasets(self, dataset_1_df, dataset_2_df):
        df = pd.concat([dataset_1_df, dataset_2_df], ignore_index=True, 
                       sort=False)
        return df
    
    
            
    def fractions_feature_identifications(self, dataset_1, dataset_2, feature, 
                                          show_data = True, show_statistics = True): 
        dataset_1_df = dataset_1.ReadPSMs.psms_df
        dataset_1_df["Dataset"] = ["Dataset_1" for i in list(dataset_1_df.index)]
        dataset_2_df = dataset_2.ReadPSMs.psms_df
        dataset_2_df["Dataset"] = ["Dataset_2" for i in list(dataset_2_df.index)]
        if feature in ["GRAVY score", "pI"]:
            if "GRAVY score" not in list(dataset_1_df.index):
                dataset_1_df = dataset_1.calculate_properties(dataset_1_df)
            if "GRAVY score" not in list(dataset_2_df.index):
                dataset_2_df = dataset_1.calculate_properties(dataset_2_df)
        dataset_name_1 = dataset_1.ReadPSMs.dataset_name
        dataset_name_2 = dataset_2.ReadPSMs.dataset_name
        self.createGraphs.plot_comparison_fractions_feature(dataset_1_df, 
                dataset_2_df, feature, labels = [dataset_name_1, dataset_name_2],
                show_data = show_data, show_statistics = show_statistics)
        


    def all_datasets(self, datasets, feature):
        all_datasets_df = []
        all_datasets_name = []
        for index, data in enumerate(datasets):
            dataset = data.ReadPSMs.psms_df
            dataset_name = data.ReadPSMs.dataset_name
            dataset["Dataset"] = ["{}".format(dataset_name) for i in list(dataset.index)]
            if feature in ["GRAVY score", "pI"]:
                dataset = dataset.calculate_properties(dataset)
            all_datasets_df.append(dataset)
            all_datasets_name.append(dataset_name)
        self.createGraphs.plot_comparison_distributions_boxplot(all_datasets_df, feature)


        
            
    def all_identifications_overlap(self, dataset_1, dataset_2):
        """ """
        dataset_name_1 = dataset_1.ReadPSMs.dataset_name
        dataset_name_2 = dataset_2.ReadPSMs.dataset_name
        #peptides
        identifications_proteins_1 = dataset_1.give_all_proteins(\
                                            dataset_1.ReadPSMs.psms_df)
        identifications_proteins_2 = dataset_2.give_all_proteins( \
                                            dataset_2.ReadPSMs.psms_df)
        self.createGraphs.plot_comparison_overlap(identifications_proteins_1, 
                identifications_proteins_2, labels = [dataset_name_1, 
                                                      dataset_name_2])
        #Peptide
        identifications_peptides_1 = dataset_1.give_all_peptides( \
                                            dataset_1.ReadPSMs.psms_df)
        identifications_peptides_2 = dataset_2.give_all_peptides( \
                                            dataset_2.ReadPSMs.psms_df) 
        self.createGraphs.plot_comparison_overlap(identifications_peptides_1, 
                identifications_peptides_2, labels = [dataset_name_1, 
                                                      dataset_name_2])
                
        
    def elution_profile_peptide(self, dataset_1, dataset_2, peptide): 
        """ """
        dataset_name_1 = dataset_1.ReadPSMs.dataset_name
        dataset_name_2 = dataset_2.ReadPSMs.dataset_name
        #peptides
        elution_1 = dataset_1.analyze_peptides(dataset_1.ReadPSMs.psms_df, 
                                               peptide, False)
        elution_2 = dataset_2.analyze_peptides(dataset_2.ReadPSMs.psms_df, 
                                               peptide, False)
        self.createGraphs.plot_comparison_elution_distribution(elution_1, 
                        elution_2, peptide, [dataset_name_1, dataset_name_2])
        
        
        
    def elution_profile_top10_peptides(self, dataset_1, dataset_2):
        """ """
        most_abundant_peptides, _ = self.MyTools.extract_most_abundant_in_dataframe(
            dataset_1.ReadPSMs.psms_df, dataset_1.ReadPSMs.sequence, 10)
        for peptide in most_abundant_peptides:
            self.elution_profile_peptide(dataset_1, dataset_2, peptide)
            
    
            
    def elution_difference_peptides(self, dataset_1, dataset_2):
        dataset_1_df = dataset_1.ReadPSMs.psms_df
        dataset_2_df = dataset_2.ReadPSMs.psms_df
        overlap_peptides = list(set(dataset_1_df[dataset_1.ReadPSMs.sequence]) & set(dataset_2_df[dataset_2.ReadPSMs.sequence]))
        all_differences = []
        for index, peptide in enumerate(overlap_peptides):
            mean_1 = dataset_1.calculate_mean_elution_fraction(dataset_1_df, peptide)
            mean_2 = dataset_1.calculate_mean_elution_fraction(dataset_2_df, peptide)
            difference = mean_2 - mean_1
            all_differences.append(difference)
            if index%500 == 0: print(index)
        self.createGraphs.plot_elution_difference(all_differences)
        
        
    def calculate_orthogonality(self, dataset_1, dataset_2):
        dataset_1_df = dataset_1.ReadPSMs.psms_df
        dataset_2_df = dataset_2.ReadPSMs.psms_df
        overlap_peptides = list(set(dataset_1_df[dataset_1.ReadPSMs.sequence]) & set(dataset_2_df[dataset_2.ReadPSMs.sequence]))
        retention_times_1 = []
        retention_times_2 = []
        for index, peptide in enumerate(overlap_peptides):
            if index%500 == 0: print(index)
            rt1_list = dataset_1_df[dataset_1_df[dataset_1.ReadPSMs.sequence] == peptide][dataset_1.ReadPSMs.retention_time]
            rt2_list = dataset_2_df[dataset_2_df[dataset_2.ReadPSMs.sequence] == peptide][dataset_2.ReadPSMs.retention_time]
            rt1 = statistics.mean(rt1_list)
            rt2 = statistics.mean(rt2_list)
            retention_times_1.append(rt1)
            retention_times_2.append(rt2)
        self.createGraphs.plot_orthogonality(retention_times_1, retention_times_2)
        # plt.figure()
        # plt.scatter(retention_times_1, retention_times_2)
        # plt.show()