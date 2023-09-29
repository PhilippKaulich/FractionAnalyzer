# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:41:12 2020

@author: Philipp
"""


# to-do:
# - implement Replicate analysis as one condition (maybe new class?)


import pandas as pd
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import statistics
import itertools

from MyTools import MyTools
from ReadPSMs import ReadPSMs
from createGraphs import createGraphs
from ReadFastaFile import ReadFastaFile
from ProteomicTools import ProteomicTools
from Exporter import Exporter


f1 = "C:\\Users\\Philipp\\OneDrive\\Documents\\Python\\PDAnalyzer\\Dataset\\20201003_XL-Mmazei_DSSO_Chymo_SEC_Superdex_HCDstepped_1h_10-23_all_PSMs.txt"
f2 = "D:\\OneDrive\\Documents\\Python\\PDAnalyzer\\Dataset\\20201003_XL-Mmazei_DSSO_Chymo_SEC_Superdex_HCDstepped_1h_10-23_all_PSMs.txt"


f3 = "Dataset\\20211119_CaCo2_PLRP_HighPH_BU_Fr_01-34_FullTryptic_ALL_PSMs.txt"


f4 = "D:\\OneDrive\\Documents\\Python\\PDAnalyzer\\Dataset\\20201003_XL-Mmazei_DSBU_SEC_Superdex_HCDstepped_1h_7-25_multiConsensus_PSMs.txt"
f5 = "C:\\Users\\Philipp\\OneDrive\\Documents\\Python\\PDAnalyzer\\Dataset\\20201003_XL-Mmazei_DSBU_SEC_Superdex_HCDstepped_1h_7-25_multiConsensus_PSMs.txt"
f6 = "Dataset\\20201003_XL-Mmazei_DSSO_Chymo_SEC_Superdex_HCDstepped_1h_10-23_all_PSMs.txt"
f7 = "T:\\PTK\\Projects\\2DLC\\lowPH_lowPH\\BU\\20211119_CaCo2_PLRP_BU_Fr_01-34_FullTryptic_ALL_PSMs.txt"
f8 = "C:\\Users\\Philipp\\OneDrive\\Documents\\Python\\PDAnalyzer\\Dataset\\20211119_CaCo2_PLRP_HighPH_BU_Fr_01-34_FullTryptic_ALL_PSMs.txt"
f9 = "C:\\Users\\Philipp\\OneDrive\\Documents\\Python\\PDAnalyzer\\Dataset\\20211202_CaCo2_PLRP_HighPH_BU_Fr_01-35_PSMs.txt"

sp = "C:\\Users\\Phili\\Desktop\\temp"
sp2 = "D:\\TEMP\\"
sp3 = "Results"

fasta_file_name = "Dataset\\20210709_Uniprot_Human_reviewed.fasta"


class PDAnalyzer:
    def __init__(self, file_name = f3, fasta_file = fasta_file_name, 
                 save_path = sp3, dataset_name = ""):
        self.MyTools = MyTools()
        self.ReadPSMs = ReadPSMs(file_name, dataset_name)
        self.createGraphs = createGraphs(save_path)
        self.ReadFastaFile = ReadFastaFile()
        self.fasta_file_name = fasta_file
        self.ProteomicTools = ProteomicTools()



    def analyze_orthogonality(self, psms_df):
        print ("Fraction Nr. must be numeric!")
        rt_list = list(psms_df[self.ReadPSMs.retention_time])
        fraction_list = list(psms_df[self.ReadPSMs.name])
        fraction_list = [int(i) for i in fraction_list]
        print(fraction_list[0:10])
        print(rt_list[0:10])
        number_of_fractions = len(self.ReadPSMs.give_file_name_list(psms_df))
        print("number_of_fractions", number_of_fractions)
        self.createGraphs.plot_orthogonality(fraction_list, rt_list, 
                                             number_of_fractions)
        
        

    def _analyze_protein_properties(self, accession_list):
        fasta_dict = self.ReadFastaFile.read_fasta_file(self.fasta_file_name)
        protein_prop = {}
        print(fasta_dict)
        for accession in accession_list:
            if accession not in fasta_dict: print(accession); continue
            sequence, mass, pI, gravy = \
                self.ReadFastaFile.analyze_protein(fasta_dict, accession)
            protein_prop[accession] = [mass, pI, gravy]
        return protein_prop



    def analyze_fractions_by_protein_feature(self, psms_df, show_plot = True):
        feature = self.ReadPSMs.master_protein
        file_name_list, protein_list_all = self._extract_feature_from_filenames( \
                                                            psms_df, feature)
        accession_list = list(set([j for i in protein_list_all for j in i]))
        protein_prop = self._analyze_protein_properties(accession_list)
        properties = {"pi":[], "gravy":[], "mass":[]}
        for file_name, protein_list in zip(file_name_list, protein_list_all):
            properties["pi"].append([protein_prop[accession][1] for \
                                    accession in protein_list \
                                    if accession in protein_prop])
            properties["gravy"].append([protein_prop[accession][2] for \
                                    accession in protein_list \
                                    if accession in protein_prop])
            properties["mass"].append([protein_prop[accession][0] for \
                                    accession in protein_list \
                                    if accession in protein_prop])
        for feature in properties:
            feature_list = properties[feature]
            if show_plot:
                self.createGraphs.plot_feature_distribution(file_name_list, 
                                                        feature_list, feature)

            
    
    def _enhanced_annotation(self, sequence):
        return "." in sequence
    
    
    def calculate_properties(self, psms_df): 
        sequence_list = list(psms_df[self.ReadPSMs.sequence])
        pi_list = []
        gravy_list = []
        for sequence in sequence_list:
            sequence = sequence.upper()
            if self._enhanced_annotation(sequence): 
                sequence = sequence.split(".")[1]
            pi_list.append(self.ProteomicTools.calculate_pI(sequence))
            gravy_list.append(self.ProteomicTools.calculate_gravy(sequence))
        df = psms_df.copy()
        df["pI"] = pi_list
        df["GRAVY score"] = gravy_list
        return df
                
    

  
    def _extract_feature_from_filenames(self, psms_df, feature):
        """ """
        file_name_list = self.ReadPSMs.give_file_name_list(psms_df)
        feature_list_all = [list(set(psms_df[psms_df[self.ReadPSMs.name] == \
                        file_name][feature])) for file_name in file_name_list]
        return file_name_list, feature_list_all
  
  
    def analyze_overlap_coefficient(self, psms_df, feature, show_plot = True):
        """ t = peptide or protein """
        file_name_list, feature_list_all = self._extract_feature_from_filenames( \
                                                            psms_df, feature)
        overlap_coefficient_df = pd.DataFrame(0.00, 
                            columns=file_name_list, index=file_name_list)
        for file_name_1, feature_list_1 in zip(file_name_list, feature_list_all):
            for file_name_2, feature_list_2 in zip(file_name_list, feature_list_all):
                overlap_coefficient = \
                    self.MyTools.calculate_overlap_coefficient(feature_list_1, 
                                                               feature_list_2)
                overlap_coefficient_df[file_name_1][file_name_2] = overlap_coefficient
        if show_plot:
            self.createGraphs.plot_overlap_heatmap(overlap_coefficient_df)
        return overlap_coefficient_df
                

    def analyze_unique_peptide_fraction(self, psms_df, feature, show_plot = True):
        """ """
        _, feature_list_all = self._extract_feature_from_filenames( \
                                                            psms_df, feature)
        all_features = [j for i in feature_list_all for j in i] 
        frequency_features = Counter(all_features)
        frequency_number_of_peptides_in_fractions =  \
            Counter(list(frequency_features.values()))
        if show_plot:
            self.createGraphs.plot_number_peptides_in_fractions(frequency_number_of_peptides_in_fractions)
        return frequency_number_of_peptides_in_fractions



    def analyze_distribution_all_identifications(self, psms_df, feature, show_plot = True):
        if feature in ["GRAVY score", "pI"] and "gravy" not in psms_df.index: 
            psms_df = self.calculate_properties(psms_df)
        feature_list = [float(str(i).replace(",", ".")) for i in list(psms_df[feature])]
        if show_plot:
            self.createGraphs.plot_distributions_histogram(feature_list, feature)
        return feature_list
        
        

    def analyze_fractions_by_feature(self, psms_df, feature = "MH+ [Da]", 
                                     show_plot = True, show_data_points= True):
        """ Analysis mass average etc. """
        #if feature not in self.ReadPSMs.all_features: 
         #   return None 
        files_name_list = self.ReadPSMs.give_file_name_list(psms_df)
        feature_list_all = []
        if feature in ["GRAVY score", "pI"]: 
            psms_df = self.calculate_properties(psms_df)
        for file_name in files_name_list:
            print (file_name)
            feature_list = self.MyTools.extract_dataframe(psms_df, 
                        self.ReadPSMs.name, file_name)[feature]
            feature_list = [float(str(i).replace(",",".")) for i in feature_list]
            feature_list_all.append(list(feature_list))
        av, med, sem, sd = self.MyTools.calculate_average_and_sd_from_multiple_lists(feature_list_all)
        if show_plot:
            self.createGraphs.plot_feature_distribution(files_name_list, feature_list_all, feature, show_data_points)
        return feature_list_all

    
    def analyze_fractions_by_identifications(self, psms_df, show_plot = True):
        """ analysis #peptides, #psms, #proteins """
        files_name_list = self.ReadPSMs.give_file_name_list(psms_df)
        # {file_name: [seq1, seq2, seq3, ..]}
        numbers_file_name_dict = {}
        cumulative_unique_sequences = []
        cumulative_unique_proteins = []
        for file_name in files_name_list:
            extracted_df = self.MyTools.extract_dataframe(psms_df,
                self.ReadPSMs.name, file_name)
            sequence_list = extracted_df[self.ReadPSMs.sequence]
            unique_sequence_list = list(set(sequence_list))
            unique_protein_list = list(set(extracted_df[self.ReadPSMs.master_protein]))
            cumulative_unique_sequences = list(set(cumulative_unique_sequences 
                                                   + unique_sequence_list))
            cumulative_unique_proteins = list(set(cumulative_unique_proteins 
                                                  + unique_protein_list))
            numbers_file_name_dict[file_name] = [len(sequence_list), 
                                                 len(unique_sequence_list),
                                                 len(unique_protein_list),
                                                 len(cumulative_unique_sequences),
                                                 len(cumulative_unique_proteins)]
        df_result = pd.DataFrame.from_dict(numbers_file_name_dict, 
                    orient='index', columns= ["#PSMs", "#Peptides","#Proteins",
                    "#CumulativePeptides","#CumulativeProteins"])
        df_result.index = files_name_list
        if show_plot:
            self.createGraphs.plot_feature_identifications(df_result)
        return df_result


    def analyze_peptides(self, psms_df, peptide, show_plot = True):
        df = psms_df[psms_df[self.ReadPSMs.sequence] == peptide]
        count = Counter(list(df[self.ReadPSMs.name]))
        all_names = sorted(list(set(list(psms_df[self.ReadPSMs.name]))))
        for name in all_names:
            if name not in count: count[name] = 0 
        count = dict(sorted(count.items()))
        if show_plot:
            self.createGraphs.plot_elution_distribution(count, peptide)
        return count


    def analyze_topn_peptide_distribution(self, psms_df, topn=5, show_plot = True):
        #häufigsten Peptide (count PSMs)
        most_abundant_peptides, _ = self.MyTools.extract_most_abundant_in_dataframe(
            psms_df, self.ReadPSMs.sequence, topn)
        all_peptides = []
        all_counts = []
        for peptide in most_abundant_peptides:
            count = self.analyze_peptides(psms_df, peptide, show_plot)
            all_peptides.append(peptide)
            all_counts.append(count)
        return all_peptides, all_counts


        


    def calculate_mean_elution_fraction(self, psms_df, peptide):
        count = self.analyze_peptides(psms_df, peptide, False)
        h1 = []
        for index, fraction in enumerate(count):
            if count[fraction] != 0:
                h1.extend([(index+1)]*count[fraction])
        mean = statistics.mean(h1)
        return mean
        #mean = statistics.mean

    def give_all_proteins(self, psms_df):
        proteins = list(set(psms_df[self.ReadPSMs.master_protein]))
        return proteins
    
    def give_all_peptides(self, psms_df): 
        peptides = list(set(psms_df[self.ReadPSMs.sequence]))
        return peptides
        
    
    
    
    
    
    
    # def analyze_overlap_coefficient(self, psms_df, t="peptide"):
    #     """ t = peptide or protein """
    #     files_name_list, feature_list_dict, _ = \
    #         self._extract_feature_to_filename(psms_df, t)
    #     #return feature_list_dict
    #     overlap_coefficient_df = pd.DataFrame(0.00, 
    #                         columns=files_name_list, index=files_name_list)
    #     for file_name1 in files_name_list:
    #         for file_name2 in files_name_list:
    #             list_file_1 = feature_list_dict[file_name1]
    #             list_file_2 = feature_list_dict[file_name2]
    #             overlap_coefficient = \
    #                 self.MyTools.calculate_overlap_coefficient(list_file_1, 
    #                                                             list_file_2)
    #             overlap_coefficient_df[file_name1][file_name2] = overlap_coefficient
    #     # self.createGraphs.createHeatmap(overlap_coefficient_df, "Overlap coefficient", 
    #     #                                 "overlap_coefficient", "", "", "Greys")
    #     self.createGraphs.plot_overlap_heatmap(overlap_coefficient_df)
    #     print(overlap_coefficient_df)
    #     return(overlap_coefficient_df)

        
        
    # def analyze_fractions_by_features_gravy_pi(self, psms_df, feature = "pi"):
    #     """ """
    #     files_name_list, feature_list_dict, feature_list_dict_all = \
    #         self._extract_feature_to_filename(psms_df, t="peptide")
    #     peptide_mass_list, peptide_pI_list, peptide_gravy_list = [], [], []
    #     for file_name in files_name_list:
    #         sequence_list = feature_list_dict[file_name]
    #         mass_list, pI_list, gravy_list = [], [], []
    #         for sequence in sequence_list:
    #             sequence = sequence.upper()
    #             if self._enhanced_annotation(sequence): sequence = sequence.split(".")[1]
    #             mass_list.append(self.ProteomicTools.calculate_mass(sequence))
    #             pI_list.append(self.ProteomicTools.calculate_pI(sequence))
    #             gravy_list.append(self.ProteomicTools.calculate_gravy(sequence))
    #         peptide_mass_list.append(mass_list)
    #         peptide_pI_list.append(pI_list)
    #         peptide_gravy_list.append(gravy_list)
    #     if feature == "pI":
    #         self.createGraphs.create_plots_from_multiple_list(files_name_list, 
    #                 peptide_pI_list, "Peptide Isoelectric Point", "Fraction", "pI")
    #         res_list = [item for sublist in peptide_pI_list for item in sublist]
    #     elif feature == "GRAVY score":
    #         res_list = self.createGraphs.create_plots_from_multiple_list(files_name_list, 
    #                 peptide_gravy_list, "Peptide GRAVY", "Fraction", "GRAVY score")
    #         res_list = [item for sublist in peptide_gravy_list for item in sublist]
    #     return res_list 





    # def analyze_unique_peptide_fraction(self):
    #     """ x: #Fraction, y: #unique peptides """
    #     all_peptides = self.ReadPSMs.give_all_peptides()
    #     _, sequence_list_dict, _ = self._extract_feature_to_filename(self.ReadPSMs.psms_df)
    #     print("Number of peptides:", len(all_peptides))
    #     sequence_list_list = list(sequence_list_dict.values())
    #     number_fractions_unique_peptides: dict = {}
    #     for peptide in all_peptides:
    #         peptide_in_number_of_fractions = self.MyTools.calculate_number_value_in_multiple_lists(sequence_list_list, peptide)
    #        # print(peptide_in_number_of_fractions)
    #         number_fractions_unique_peptides = \
    #             self.MyTools.add_value_to_dict(number_fractions_unique_peptides, 
    #                                            peptide_in_number_of_fractions, 1)
    #     self.createGraphs.create_bar_plot_from_dict(number_fractions_unique_peptides, 
    #                 "Number of peptides identified in # fractions", 
    #                 "unique_peptides", "#Fractions", "#Peptides")
    #     return number_fractions_unique_peptides
    
    
    
    # def _extract_feature_to_filename(self, psms_df, t="peptide"):
    #     """ feature: peptide or protein """
    #     feature = self.ReadPSMs.translate_feature(t)
    #     files_name_list = self.ReadPSMs.give_file_name_list(psms_df)
    #     # {file_name: [seq1, seq2, seq3, ..]}
    #     feature_list_dict = {}
    #     feature_list_dict_all = {}
    #     for file_name in files_name_list:
    #         feature_list = self.MyTools.extract_dataframe(psms_df,
    #             self.ReadPSMs.name, file_name)[feature]
    #         feature_list_dict_all[file_name] = feature_list
    #         unique_feature_list = list(set(feature_list))
    #         feature_list_dict[file_name] = unique_feature_list
    #     #print(feature_list_dict)
    #     return files_name_list, feature_list_dict, feature_list_dict_all

        



    # def analyze_fractions_by_proteins(self):
    #     fasta_file_dict = self.ReadFastaFile.read_fasta_file(self.fasta_file_name)
    #     fasta_file = self.ReadFastaFile.analyze_all_proteins(fasta_file_dict)
    #     files_name_list, feature_list_dict, feature_list_dict_all = \
    #         self._extract_feature_to_filename(self.ReadPSMs.psms_df, t="protein")
    #     protein_mass_list, protein_pI_list, protein_gravy_list = [], [], []
    #     for file_name in files_name_list:
    #         accession_list = feature_list_dict[file_name]
    #         mass_list, pI_list, gravy_list = [], [], []
    #         for accession in accession_list:
    #             if accession not in fasta_file:  
    #                 print(accession, "nicht in Datanbank vorhanden")
    #                 continue
    #             sequence, mass, pI, gravy = fasta_file[accession]
    #             mass_list.append(mass)
    #             pI_list.append(pI)
    #             gravy_list.append(gravy)
    #         protein_mass_list.append(mass_list)
    #         protein_pI_list.append(pI_list)
    #         protein_gravy_list.append(gravy_list)
    #     self.createGraphs.create_plots_from_multiple_list(files_name_list, 
    #             protein_mass_list, "Protein Mass", "Fraction", "Mass / Da")
    #     self.createGraphs.create_plots_from_multiple_list(files_name_list, 
    #             protein_pI_list, "Protein Isoelectric Point", "Fraction", "pI")
    #     self.createGraphs.create_plots_from_multiple_list(files_name_list, 
    #             protein_gravy_list, "Protein GRAVY", "Fraction", "GRAVY score")
    #     return protein_mass_list, protein_pI_list, protein_gravy_list


    
    # def temp(self):
    #     #print(len(self.ReadPSMs.give_all_peptides()))
    #     df = self.ReadPSMs.psms_df
    #     feature = "MH+ [Da]"
    #     df[feature] = [float(str(i).replace(",",".")) for i in list(df[feature])]
    #     sns.boxplot(data = df, x ="name", y = feature)
    
    
    
    # def analyze_peptide_distribution(self, psms_df, topn=5):
    #     #häufigsten Peptide (count PSMs)
    #     most_abundant_peptides, _ = self.MyTools.extract_most_abundant_in_dataframe(
    #         psms_df, self.ReadPSMs.sequence, topn)
    #     peptide_distribution_df = pd.DataFrame(0, 
    #                             columns=self.ReadPSMs.give_file_name_list(psms_df), 
    #                             index=most_abundant_peptides)
    #     for peptide in most_abundant_peptides:
    #         spectrum_files_where_peptide_was_identified = \
    #             self.MyTools.extract_dataframe(psms_df, 
    #             self.ReadPSMs.sequence, peptide)[self.ReadPSMs.name]
    #         count_spectrum_files = Counter(spectrum_files_where_peptide_was_identified)
    #         # Create BarPlots
    #         self.createGraphs.create_bar_plot_from_dict(dict(count_spectrum_files), 
    #                 "Peptide Distribution", 
    #                 "{}_peptide_distribution".format(peptide), 
    #                 "#Fractions", "#PSMs")
    #         peptide_distribution_df.loc[peptide] = count_spectrum_files
    #         print(count_spectrum_files)
    #     peptide_distribution_df = peptide_distribution_df.fillna(0)
        
    #     return peptide_distribution_df
    
    
    
        # for peptide in most_abundant_peptides:
        #     spectrum_files_where_peptide_was_identified = \
        #         self.MyTools.extract_dataframe(psms_df, 
        #         self.ReadPSMs.sequence, peptide)[self.ReadPSMs.name]
        #     count_spectrum_files = Counter(spectrum_files_where_peptide_was_identified)
        #     # Create BarPlots
        #     self.createGraphs.create_bar_plot_from_dict(dict(count_spectrum_files), 
        #             "Peptide Distribution", 
        #             "{}_peptide_distribution".format(peptide), 
        #             "#Fractions", "#PSMs")
        #     peptide_distribution_df.loc[peptide] = count_spectrum_files
        #     print(count_spectrum_files)
        # peptide_distribution_df = peptide_distribution_df.fillna(0)