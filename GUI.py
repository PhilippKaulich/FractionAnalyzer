# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 20:41:57 2021

@author: Phili
"""

import pandas as pd
import sys
import pickle
from collections import Counter
import os

from PyQt5 import QtWidgets, QtCore, uic
from PyQt5.QtWidgets import QWidget, QAbstractItemView, QFileDialog, QTextEdit, QVBoxLayout, QScrollArea, QListWidgetItem, QCompleter, QTreeWidgetItem, QTableWidgetItem

from PDAnalyzer import PDAnalyzer
from CompareData import CompareData
from Settings import Settings
from MyTools import MyTools
#from ReadPSMs import ReadPSMs

# enable scaling high DPI mode 
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True) 
# if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
#     QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

# if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
#     QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)



# fasta = "Dataset\\PTK_Ecoli_FASTA_all.fasta"
# a = AnalyzeCleavageSites(fasta)
# r = ReadProteoforms.ReadProteoforms(pd_version="v4.0")
# df = r.read_file("Dataset\\test_3_Proteoforms.txt")
# aa_df = a.find_subsequence_termini(df["Sequence"], df["# PrSMs"])



class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('GUI_new.ui', self)
        title = "Fraction Analyzer"
        self.setWindowTitle(title)
        self.show()
        
        self.MyTools = MyTools()
        
        ##
        
        self.combo_analysis.setEnabled(False)
        self.combo_comparison.setEnabled(False)

        
        self.table_file_names.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        #self.table_file_names.resizeColumnsToContents()
        
        ###
        self.data_sets = {}
        self.data_sets_file_names = {}
        
        
        ## Menu
        self.menue_file_about.triggered.connect(self.show_about)


    def temp(self):
        print("TEMP")
        ### read table
        print(self._read_table())
        self._set_replicates_and_conditions()
        
        print("""""")
        
        
    def selection_show_plot(self, is_checked):
        if not is_checked:
            self.check_show_data.setChecked(False)


    def _return_number_of_all_file_names(self):
        n = 0
        for e in self.data_sets_file_names:
            n += len(self.data_sets_file_names[e])
        return n

        
        
    def _read_table(self):
        nb_row = self._return_number_of_all_file_names() 
        nb_col = 4
        table = [[]]
        for row in range (nb_row):
            table.append([])
            for col in range(nb_col):
                txt = self.table_file_names.item(row, col).text()
                table[row].append(txt)
        return table
        
    
    def _set_replicates_and_conditions(self):
        table = self._read_table()
        all_replicates = []
        all_conditions = []
        for element in table: 
            if not element: continue
            raw_file, name, replicate, condition = element
            for data_set in self.data_sets_file_names:
                if raw_file in self.data_sets_file_names[data_set]:
                    dataset = self.data_sets[data_set]
                    dataset.ReadPSMs.set_replicate(raw_file, replicate)
                    dataset.ReadPSMs.set_condition(raw_file, condition)
                    dataset.ReadPSMs.set_name(raw_file, name)
                    all_replicates.append(replicate)
                    all_conditions.append(condition)
        previous_conditions = [self.combo_selection_condition.itemText(i) for i in range(self.combo_selection_condition.count())]
        previous_replicates = [self.combo_selection_replicate.itemText(i) for i in range(self.combo_selection_replicate.count())]
        
        for r in list(set(all_replicates)):
            if r not in previous_replicates:
                self.combo_selection_replicate.addItem(r)
        for c in list(set(all_conditions)):
            if c not in previous_conditions:
                self.combo_selection_condition.addItem(c)

        
    def assign_replicate_conditions(self):
        self._set_replicates_and_conditions()

    
    def delete_all_replicate(self):
        self.list_file_names.clear()
        self.initialization()
        
    def delete_all_condition(self):
        self.list_file_names.clear()
        self.initialization()
        
        
    def _filter_psms(self, psms_df, condition, replicate):
        if condition not in ["All", "all", None, "NONE"]:
            psms_df = psms_df[psms_df["condition"] == condition]
        if replicate not in ["ALL","all", None, "NONE"]:
            psms_df = psms_df[psms_df["replicate"] == replicate]
        return psms_df
        
        
    def process_selection(self):
        data_set_file = self.combo_selection_dataset.currentText()
        analysis = self.combo_selection_analysis.currentText()
        condition = self.combo_selection_condition.currentText()
        replicate = self.combo_selection_replicate.currentText()
        data_set = self.data_sets[data_set_file]
        psms_df = data_set.ReadPSMs.psms_df
        psms_df = self._filter_psms(psms_df, condition, replicate)
        show_plot = self.check_show_plot.isChecked()
        show_data_points = self.check_show_data.isChecked()
        export_data = self.check_export.isChecked()
        if analysis == Settings.fractions_identifications: 
            result = data_set.analyze_fractions_by_identifications(psms_df, 
                                                                   show_plot)
        if analysis == Settings.fractions_peptide_overlap: 
            result = data_set.analyze_overlap_coefficient(psms_df, 
                                        data_set.ReadPSMs.sequence, show_plot)
        if analysis == Settings.fractions_peptide_overlap: 
            result = data_set.analyze_overlap_coefficient(psms_df, 
                                    data_set.ReadPSMs.master_protein, show_plot)
        if analysis == Settings.peptides_in_number_fractions: 
            result = data_set.analyze_unique_peptide_fraction(psms_df,
                                        data_set.ReadPSMs.sequence, show_plot)
        if analysis == Settings.fractions_elution_peptides_top10: 
            result = data_set.analyze_topn_peptide_distribution(psms_df, 
                                                                10, show_plot)
        if analysis == Settings.fractions_feature_distribution: 
            feature = self.combo_analysis.currentText() 
            result = data_set.analyze_fractions_by_feature(psms_df, feature,
                                                    show_plot,show_data_points)
            # m = MyTools()
            # m.export_multiple_lists(result, "test.csv")
        if analysis == Settings.all_feature_distribution: 
            feature = self.combo_analysis.currentText() 
            result = data_set.analyze_distribution_all_identifications(psms_df, 
                                                            feature, show_plot)
        if analysis == Settings.fractions_elution_peptide:
            peptide = self.combo_analysis.currentText() 
            result = data_set.analyze_peptides(psms_df, peptide, show_plot)
        if analysis == Settings.fractions_protein_feature_distribution:
            data_set.analyze_fractions_by_protein_feature(psms_df, show_plot)
            
        if analysis == Settings.analyze_orthogonality:
            data_set.analyze_orthogonality(psms_df)
        # if analysis=="Total Protein Mass Distribution": result = data_set.
        # if analysis=="Total Peptide Mass Distribution": result = data_set.
        # if analysis=="Fractions Elution Peptide": result = data_set.
        # if analysis=="Fractions Elution Protein": result = data_set.
        # if analysis=="Analyze All Distributions": result = data_set.
        # if analysis=="Analyze All": result = data_set.
        # if analysis == "Fractions Identified Peptides": 
        #     result = data_set.analyze_fractions_by_identifications()
        # if analysis == "Fractions Identified Proteins": 
        #     result = data_set.analyze_fractions_by_identifications()
        # if analysis == "Fractions Identified PSMs": 
        #     result = data_set.analyze_fractions_by_identifications()


# analyze_fractions_by_proteins()
# analyze_fractions_by_peptides()



# analyze_fractions_by_feature()
# analyze_fractions_by_identifications()


        
    def initialization(self):
        dataset_name = self.text_psm_file_name.toPlainText()
        save_path = self.text_save_path.toPlainText()
        # add data set to combobox
        if dataset_name not in self.data_sets_file_names:        
            pd_analyzer = PDAnalyzer(file_name = self.psm_file_name, 
                                fasta_file = self.fasta_file_name, 
                                save_path = save_path,
                                dataset_name = dataset_name)
            self.data_sets[dataset_name] = pd_analyzer
            
            self.combo_selection_dataset.addItem(dataset_name)
            self.combo_comparison_1.addItem(dataset_name)
            self.combo_comparison_2.addItem(dataset_name)
            # add file name lists to text box
            file_names = pd_analyzer.ReadPSMs.give_file_name_list(pd_analyzer.ReadPSMs.psms_df)
            self.data_sets_file_names[dataset_name]= file_names
            for index, file_name in enumerate(file_names):
                self.table_file_names.insertRow(index)
                raw_file  = QTableWidgetItem()
                raw_file.setText(file_name)
                na = QTableWidgetItem()
                na.setText("NONE")
                na2 = QTableWidgetItem()
                na2.setText("NONE")
                na3 = QTableWidgetItem()
                na3.setText(file_name)
                self.table_file_names.setItem(index, 0, raw_file)
                self.table_file_names.setItem(index, 1, na3)
                self.table_file_names.setItem(index, 2, na)
                self.table_file_names.setItem(index, 3, na2)
                # try:
                #     self.table_file_names.setItem(index, 0).setText(file_name)
                # except:
                #     print("passt wirgendwas nicht")
            # Multiple selection
            self.table_file_names.resizeColumnsToContents()
            self.table_file_names.setSelectionMode(2)
            self._set_replicates_and_conditions()
        
        
            
        
        
        

    def selection_analysis(self):
        data_set_file = self.combo_selection_dataset.currentText()
        data_set = self.data_sets[data_set_file]
        analysis = self.combo_selection_analysis.currentText()
        self.combo_analysis.setEnabled(True)
        # show combobox if necessary
        if analysis == Settings.fractions_feature_distribution or \
                analysis == Settings.all_feature_distribution: 
            possible_features = data_set.ReadPSMs.give_all_digit_columns(data_set.ReadPSMs.psms_df)
            possible_features.extend([Settings.gravy, Settings.pi])
            self.combo_analysis.clear()
            self.combo_analysis.addItems(possible_features)
        # show all possible peptides
        elif analysis == Settings.fractions_elution_peptide:
            possible_peptides = data_set.ReadPSMs.give_all_peptides()   
            self.combo_analysis.clear()
            self.combo_analysis.addItems(possible_peptides)
        else: 
            self.combo_analysis.clear()
            self.combo_analysis.setEnabled(False)
            
                        
            
            
    def comparison_analysis(self):
        data_set_file_1 = self.combo_comparison_1.currentText()
        data_set = self.data_sets[data_set_file_1]
        analysis = self.combo_comparison_analysis.currentText()
        self.combo_comparison.setEnabled(True)
        # show combobox if necessary
        if analysis == Settings.all_feature_distribution or \
                analysis == Settings.fractions_feature_distribution or \
                analysis == Settings.all_feature_distribution_all_datasets or \
                analysis == Settings.compare_all_datasets:
            possible_features = data_set.ReadPSMs.give_all_digit_columns(data_set.ReadPSMs.psms_df)
            possible_features.extend([Settings.gravy, Settings.pi])
            self.combo_comparison.clear()
            self.combo_comparison.addItems(possible_features)
        elif analysis == Settings.fractions_elution_peptide:
            possible_peptides = data_set.ReadPSMs.give_all_peptides()   
            self.combo_comparison.clear()
            self.combo_comparison.addItems(possible_peptides)
            
        else: 
            self.combo_comparison.clear()
            self.combo_comparison.setEnabled(False)
            


    def comparison_process(self):
        compare_data = CompareData()
        data_set_file_1 = self.combo_comparison_1.currentText()
        data_set_file_2 = self.combo_comparison_2.currentText()
        data_set_1 = self.data_sets[data_set_file_1]
        data_set_2 = self.data_sets[data_set_file_2]
        analysis = self.combo_comparison_analysis.currentText()
        show_data = self.check_comparison_data.isChecked()
        show_statistics = self.check_comparison_statistic.isChecked()
        if analysis == Settings.all_feature_distribution:
            feature = self.combo_comparison.currentText()
            compare_data.all_feature_distributions(data_set_1, data_set_2, feature)
        elif analysis ==  Settings.fractions_feature_distribution:
            feature = self.combo_comparison.currentText()
            compare_data.fractions_feature_identifications(data_set_1, 
                    data_set_2, feature, show_data, show_statistics)
        elif analysis == Settings.fractions_elution_peptide:
            peptide = self.combo_comparison.currentText()
            compare_data.elution_profile_peptide(data_set_1, data_set_2, peptide)
        elif analysis == Settings.fractions_elution_peptides_top10:
            compare_data.elution_profile_top10_peptides(data_set_1, data_set_2)
        elif analysis == Settings.peptide_elution_difference:
            compare_data.elution_difference_peptides(data_set_1, data_set_2)
        elif analysis == Settings.all_feature_distribution_all_datasets:
            feature = self.combo_comparison.currentText()
            all_datasets = list(self.data_sets.values())
            compare_data.all_datasets(all_datasets, feature)
        elif analysis == Settings.all_identifications_overlap:
            compare_data.all_identifications_overlap(data_set_1, data_set_2)
        elif analysis == Settings.number_identifications:
            compare_data.fractions_identifications(data_set_1, data_set_2)
        elif analysis == Settings.analyze_orthogonality:
            compare_data.calculate_orthogonality(data_set_1, data_set_2)
        elif analysis == Settings.compare_all_datasets:
            feature = self.combo_comparison.currentText()
            all_datasets = list(self.data_sets.values())
            compare_data.all_datasets(all_datasets, feature)
        else: 
            pass

        
        
        
        

        
    def load_fasta_file(self):
        self.fasta_file_name = QFileDialog.getOpenFileName(caption = Settings.fata_file_open_dialog,
            directory = '', filter = Settings.fasta_file_filter, 
            initialFilter = Settings.fasta_file_initial_filter)[0]
        self.text_fasta_file.setText(self.fasta_file_name.split("/")[-1])

        
    
    def load_psm_file(self):
        self.psm_file_name = QFileDialog.getOpenFileName(caption=Settings.psm_file_open_dialog,
            directory = '', filter = Settings.psm_file_filter, 
            initialFilter = Settings.psm_file_initial_filter)[0]
        path_name, file_name = os.path.split(self.psm_file_name) 
        self.text_psm_file.setText(file_name)
        self.text_psm_file_name.setText(file_name)
        save_path = path_name + Settings.save_path
        self.text_save_path.setText(save_path)
        self.MyTools.create_folder(save_path)
    
    
    def load_save_path(self):
        save_path = QFileDialog.getExistingDirectory(caption=Settings.psm_file_open_dialog,
            directory = '')
        self.MyTools.create_folder(save_path)
        self.text_save_path.setText(save_path)
        
    
    
    def change_level(self):
        level = self.combobox_level.currentText()
    
    



    def export_data(self):
        file_name = QFileDialog.getSaveFileName(caption = Settings.export_file_open_dialog,
            directory = '', filter = Settings.export_file_filter, 
            initialFilter = Settings.export_file_initial_filter)[0]
        self.potential_cleavage.to_csv(file_name)
        print(file_name, "erfolgerich gespeichert")
        
        

    def show_about(self):
        text = """
            Cleavage Analysis TDP
            {} ({})
            by Philipp T. Kaulich 
            
            Manual: -
            Contact: p.kaulich@iem.uni-kiel.de        
        """.format(Settings.version, Settings.date)
        self.w = ResultWindow()
        self.w.textedit.setText(text)
        self.w.show()


class ResultWindow(QWidget):
    """ New window used for peak detection results """
    def __init__(self):
        super().__init__()
        self.setFixedSize(640, 480)
        self.setWindowTitle("Peaklist - Results")
        layout = QVBoxLayout()
        self.textedit = QTextEdit()
        layout.addWidget(self.textedit)
        self.setLayout(layout)


app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()