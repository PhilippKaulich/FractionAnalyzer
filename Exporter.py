# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 20:00:03 2022

@author: Phili
"""

import pandas as pd 


class Exporter:
    def __init__(self, save_path):
        self.save_path = save_path + "/"
    
    def export_dataframe(self, df, file_name):
        """ """
        df.to_csv(self.save_path + file_name)
        
    def export_list(self, l, file_name):
        """ """
        df = pd.DataFrame(data={" ": l})
        df.to_csv(self.save_path + file_name, sep=',',index=False)
        
    def export_list_of_lists(self, ll, file_name, column_names = []):
        """ """
        number_of_lists = len(ll)
        maximal_number_elements = max([len(l) for l in ll])
        if not column_names: column_names = [i for i in range(number_of_lists)]
        new_l = []
        for index, l in enumerate(ll):
            if len(l) < maximal_number_elements:
                l.extend([" "]*(maximal_number_elements-len(l)))
            new_l.append(l)

        data = {column_names[index]: new_l[index] for index in range(number_of_lists)}

        df = pd.DataFrame(data=data)
        save_name = "{}{}.csv".format(self.save_path, file_name)
        df.to_csv(save_name, sep=',',index=False)
        print("SAVED", self.save_path + file_name)
        
        
        
    def export_dict(self, d, file_name):
        df = pd.DataFrame(data = list(d.items()))
        df.to_csv(self.save_path + file_name, sep=',', index=False, header=False)
        