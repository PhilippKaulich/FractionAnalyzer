# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 21:46:25 2021

@author: Phili
"""

from collections import Counter
import statistics
import scipy.stats as stats
import math
import csv
import os 


class MyTools:
    def __init__(self):
        pass
    
    
    def export_multiple_lists(self, ll, file_path):
        with open(file_path, "w") as file:
            writer = csv.writer(file)
            writer.writerows(ll)
    
    
    
    def create_folder(self, path):
        if not os.path.exists(path):
            os.makedirs(path)
    
    
    
    def calculate_cohens_d(self, l1, l2):
        if len(l1) < 1 and len(l2) < 1: return 0
        average_l1 = statistics.mean(l1)
        std_l1 = statistics.stdev(l1)
        average_l2 = statistics.mean(l2)
        std_l2 = statistics.stdev(l2)
        n_l1 = len(l1)
        n_l2 = len(l2)
        std = math.sqrt( ( (n_l1 - 1) * std_l1**2 + (n_l2 - 1) * std_l2**2 ) \
                        /(n_l1 + n_l2 - 2 ))
        cohens_d = (average_l2 - average_l1) / std
        return cohens_d
                # av_dataset_1 = df_1_fr[feature].mean() 
                # av_dataset_2 = df_2_fr[feature].mean()
                # sd_dataset_1 = df_1_fr[feature].std()   
                # sd_dataset_2 = df_2_fr[feature].std()
                # n_1 = len(df_1_fr[feature])
                # n_2 = len(df_2_fr[feature])
                # sd_gepoolt = ((((len(df_1_fr)-1)*sd_dataset_1**2 + (len(df_2_fr)-1)*sd_dataset_2**2))/(len(df_1_fr) + len(df_2_fr) - 2 ))**(1/2)                 
                # cohens_d = (av_dataset_2 - av_dataset_1) / sd_gepoolt


    def string_is_number(self, string):
        return str(string).replace(",","").replace(".","").isdigit()


    def calculate_overlap_coefficient(self, l1, l2):
        overlap = list(set(l1) & set(l2))
        minimal_length = len(l1) if len(l1) < len(l2) else len(l2)
        if minimal_length == 0: return 0
        overlap_coefficient = len(overlap) / minimal_length
        return overlap_coefficient
    
    
    
    def calculate_overlap(self, l1, l2): 
        overlap = overlap = list(set(l1) & set(l2))
        unique_l1 = [e for e in list(set(l1)) if e not in overlap]
        unique_l2 = [e for e in list(set(l2)) if e not in overlap]
        return unique_l1, unique_l2, overlap
        
    
    
    def extract_dataframe(self, df, column, value):
        return df[df[column] == value]
    
    def extract_most_abundant_in_dataframe(self, df, column, topn):
        count_dict = Counter(list(df[column]))
        print(count_dict)
        most_abundant_elements = count_dict.most_common(topn)
        keys, values = zip(*most_abundant_elements)
        return list(keys), list(values)
    
    def _convert_str_list_in_float_list(self, l):
        if not l: return []
        elif type(l[0]) == float: return l
        elif type(l[0]) == int: return l
        else: return [float(element.replace(",",".")) for element in l]
        
    def calculate_average_and_sd_from_multiple_lists(self, l):
        """ l = [[...], [...], [...] """
        av, med, sd, sem = [], [], [], []
        for current_list in l:
            current_list = self._convert_str_list_in_float_list(current_list)
            if len(current_list):
                av.append(statistics.mean(current_list))
                med.append(statistics.median(current_list))
                sem.append(stats.sem(current_list))
                sd.append(statistics.stdev(current_list))
            else: 
                av.append(0)
                med.append(0)
                sd.append(0)
                sem.append(0)
        return av, med, sem, sd
    
    def sort_dictionary(self, d):
        return dict(sorted(d.items()))
    
    def calculate_number_value_in_multiple_lists(self, l, value):
        """ """
        count = 0
        for current_list in l:
            if value in current_list: count += 1
        return count
    
    
    def add_value_to_dict(self, d, key, value):
        if key not in d: d[key] = value
        else: d[key] += value
        return d