# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 13:55:41 2021

@author: Phili
"""



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import numpy as np
import scipy.stats

import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn2_circles
from statannot import add_stat_annotation
from Settings import Settings
from Exporter import Exporter

import MyTools




class createGraphs:
    """ """
    def __init__(self,save_path):
        self.path = save_path
        self.MyTools = MyTools.MyTools()
        self.Exporter = Exporter(save_path)
    
    
    
    def plot_orthogonality(self, x, y, number_of_bins = 100):
        x = np.array(x)
        y = np.array(y)
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
        plt.figure(dpi=120)
        plt.plot(x, y, 'o', color='black')
        plt.plot(x, slope * x + intercept, color='red')
        # ticks  =[i for i in range(number_of_bins + 1)]
        # plt.xticks(ticks, rotation='horizontal')
        plt.show()
        print("R²=", r_value)
        
        res = np.histogram2d(x, y, bins = number_of_bins)
        hist2d = pd.DataFrame(res[0])
        x_edges, y_edges = res[1].round(0), res[2].round(0)
        print(hist2d)
        plt.figure(dpi=120)
        sns.heatmap(hist2d, cmap=Settings.heatmap_colormap)
        ticks  =[i for i in range(number_of_bins + 1)]
        # plt.xticks(ticks, x_edges, rotation='horizontal')
        # plt.yticks(ticks, y_edges, rotation='horizontal')
        plt.xlabel("1st Dimension")
        plt.ylabel("2nd Dimension")
        plt.show()
        number_of_occupied_bins = np.count_nonzero(hist2d)
        print(number_of_occupied_bins)
        orthogonality = (number_of_occupied_bins - number_of_bins) / (0.63 * 10000)
        print(orthogonality)
    
    
    
    def _combine_datasets(self, dataset_1_df, dataset_2_df):
        df = pd.concat([dataset_1_df, dataset_2_df], ignore_index=True, sort=False)
        return df
    
    
    def plot_number_peptides_in_fractions(self, count_dict, export = True):
        plt.figure(dpi=120)
        count_dict = dict(sorted(count_dict.items()))
        x, y = list(count_dict.keys()), list(count_dict.values())
        ax = sns.barplot(x=x, y=y, color="black")       
        patches = ax.patches
        percentage = [i/sum(y)*100 for i in y]
        for i in range(len(patches)):
            x1 = patches[i].get_x() + patches[i].get_width()/2
            y1 = patches[i].get_height() + 50
            ax.annotate('{:.1f}%'.format(percentage[i]), (x1, y1), ha='center')
        plt.xlabel("#Fractions")
        plt.ylabel("Count")
        plt.show()
        if export: 
            self.Exporter.export_dict(dict(count_dict), 
                                        Settings.peptides_in_number_fractions)
    
    
    
    
    def plot_elution_difference(self, difference_list):
        plt.figure(dpi=120)
        sns.histplot(difference_list)
        plt.show()
    
    
    
    
    def plot_comparison_distributions_histogram(self, data_1, data_2, feature, 
                                        density = False, labels = ["1", "2"],
                                        show_boxplot = False):
        plt.figure(dpi=120)
        sns.reset_orig()
        n, bins, patches = plt.hist(data_1, 50, density=density, 
                                color="black", alpha=0.55,label=labels[0])
        n, bins, patches = plt.hist(data_2, 50, density=density, 
                                color="orange", alpha=0.55, label=labels[1])
        plt.legend()
        plt.xlabel(feature)
        ylabel = "Density" if density else "Count"
        plt.ylabel(ylabel)
        plt.title('Histogram of {}'.format(feature))
        plt.grid(True)
        plt.show()
        
        if show_boxplot:
            data = [data_1, data_2]
            plt.figure(dpi=120)
            # palette = "colorblind"
            bp = sns.boxplot(data = data, showfliers = False)
            # sp = sns.stripplot(data = data, alpha=0.1, size=2, dodge=True)
            plt.show()

            
        
        
        
        
        
    def _combine_multiple_datasets(self, dataset_list):
        df = pd.concat(dataset_list, ignore_index=True, sort=False)
        return df
        
    
    def plot_comparison_distributions_boxplot(self, datasets, feature):
        """ """
        combined_df = self._combine_multiple_datasets(datasets)
        print(combined_df.head(10))
        plt.figure(dpi=120)
        x = "Dataset"
        y = feature
        # palette = "colorblind"
        bp = sns.boxplot(data = combined_df, x=x, y=y,
           showfliers = False)
        sp = sns.stripplot(data = combined_df, x=x, y=y,
            alpha=0.1, size=2, dodge=True)
        plt.xlabel("Dataset")
        plt.ylabel(feature)
        plt.show()
        
        
    def plot_comparison_fractions_identifications(self, list_1, list_2, 
                                        name_list, ident, labels = ["1", "2"]):
        # Peptides, #PSMs
        index = name_list
        x = np.arange(len(index))  # the label locations
        width = 0.35  # the width of the bars
        fig, ax = plt.subplots(dpi=120)
        rects1 = ax.bar(x - width/2, list_1, width, label=labels[0], color="black")
        rects2 = ax.bar(x + width/2, list_2, width, label=labels[1], color="red")
        # rects2 = ax.scatter(index, df["#Proteins"], label='#Proteins', color="blue")
        plt.ylabel(ident)
        plt.xlabel("Fraction")
        plt.xlim([-0.5, len(index) - 0.5])
        plt.legend()
        plt.show()
        
        
    
    
    def plot_comparison_fractions_feature(self, df_1, df_2, feature, 
                                          labels = ["1", "2"], 
                                          show_data = True, 
                                          show_statistics = True):
        combined_df = self._combine_datasets(df_1, df_2)
        plt.figure(dpi=120)
        order = sorted(list(set(list(df_1["name"]))))
        x = "name"
        y = feature
        hue = "Dataset"
        palette = "colorblind"
        bp = sns.boxplot(data = combined_df, x=x, y=y,
            hue=hue, palette = palette,  order=order, showfliers = False)
        if show_data:
            sp = sns.stripplot(data = combined_df, x=x, y=y,
                hue=hue, palette = palette, alpha=0.1, size=2, 
                dodge=True, order=order)
        plt.xlabel("Fraction")
        plt.ylabel(feature)
        handles, _ = bp.get_legend_handles_labels()
        plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(0,1.02,1,0.2), 
                   loc="lower left", borderaxespad=0, title = "", ncol=3)
        fractions = order
        datasets = ["Dataset_1", "Dataset_2"]
        box_pairs = []
        cohens_ds = []
        for fraction in fractions:
                box_pairs.append(((fraction, datasets[0]), (fraction, datasets[1])))
                feature_list_1 = list(df_1[df_1["name"] == fraction][feature])
                feature_list_2 = list(df_2[df_2["name"] == fraction][feature])
                cohens_d = self.MyTools.calculate_cohens_d(feature_list_1, feature_list_2)
                cohens_ds.append(cohens_d)
        if show_statistics:
            ax, test_result_list = add_stat_annotation(bp, data=combined_df, x=x,
                        y=y, hue=hue, order=order,
                        box_pairs=box_pairs, pvalue_thresholds = Settings.pvalue_thresholds, 
                        test='t-test_ind', text_format='star', loc='inside', verbose=2)
            for index, r in enumerate(test_result_list):
                print(box_pairs[index], r, cohens_ds[index])
            
        plt.show()
        
    
    
    def plot_comparison_overlap(self, l1, l2 , labels = ["1", "2"]):
        plt.figure(dpi=120)
        unique_1, unique_2, overlap = self.MyTools.calculate_overlap(l1, l2)
        subset = (len(unique_1), len(unique_2), len(overlap))
        subset_labels = (labels[0], labels[1])
        total = sum(subset)
        d = venn2(subsets = subset, set_labels = subset_labels, 
              set_colors=('orange', 'LightBlue'), alpha = 0.4,
              subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/total):1.0%}" + ")")
    #    fonzsize = 24
    #    font = 'serif'
    #    for p in ['01', '11', '10']:
    #        label = d.get_label_by_id(p)
    #        label.set_fontsize(fonzsize) 
    #        label.set_family(font)
        # Label size
    #    d.get_label_by_id('A').set_fontsize(22)
    #    d.get_label_by_id('B').set_fontsize(22)
        venn2_circles(subsets = subset, linewidth=1, color='k');
        plt.show()
    

        
    def plot_feature_distribution(self, name_list, feature_distribution, 
                                  feature, show_data_points = True, 
                                  export = True):
        plt.figure()
        df = pd.DataFrame(feature_distribution, index=name_list)
        df = df.replace(",",".")
        df = df.T
        sns.boxplot(data = df, showfliers = False)#, color="white")
        if show_data_points:
            sns.stripplot(data = df, alpha=0.1, size=2)
        plt.xlabel("Fraction")
        plt.ylabel(feature)
        plt.show()
        if export: 
            self.Exporter.export_list_of_lists(feature_distribution,
                                        Settings.fractions_feature_distribution,
                                        name_list)
        
    
    def plot_elution_distribution(self, count_dict, feature):
        plt.figure(dpi=120)
        x, y = list(count_dict.keys()), list(count_dict.values())
        sns.barplot(x=x, y=y, color="black")        
        plt.xlabel("Fraction")
        plt.ylabel(feature)
        plt.show()
        
        
    def plot_comparison_elution_distribution(self, count_dict_1, count_dict_2, feature, labels):
        plt.figure(dpi=120)
        x1, y1 = list(count_dict_1.keys()), list(count_dict_1.values())
        x2, y2 = list(count_dict_2.keys()), list(count_dict_2.values())
        if x1 != x2: return False
        X_axis = np.arange(len(x1))
        plt.bar(X_axis - 0.2, y1, 0.4, label = labels[0], color="black")
        plt.bar(X_axis + 0.2, y2, 0.4, label = labels[1], color="red")
        plt.legend()
        plt.xlabel("Fraction")
        plt.ylabel(feature)
        plt.show()
        
    
    
    def plot_overlap_heatmap(self, df):
        plt.figure(dpi=120)
        index = [i for i in range(1, len(df.index)+1)]
        columns = [i for i in range(1, len(df.index)+1)]
        annot = True if len(df.index)<11 else False
        plot = sns.heatmap(df*100, vmin=0, vmax=100, cmap="Greys", annot = annot,
                    cbar=True, xticklabels=True, yticklabels=True, fmt = '.0f',
                    cbar_kws={"shrink": 0.8, "ticks":[0,50,100], 
                              'format': '%.0f%%', "label": "Overlap Coefficient"})
        for t in plot.texts: t.set_text(t.get_text() + " %")
        plt.xticks(x = list(index), rotation=90)
        plt.yticks(y= list(index), rotation=0)
        plt.xlabel("Fraction", fontdict={'size': 12})
        plt.ylabel("Fraction", fontdict={'size': 12})
        
        plt.show()
    
    
    def plot_feature_identifications(self, df):
        # Peptides, #PSMs
        index = [str(i) for i in range(1, len(df.index)+1)]
        x = np.arange(len(index))  # the label locations
        width = 0.35  # the width of the bars
        fig, ax = plt.subplots(dpi=120)
        rects1 = ax.bar(x - width/2, df["#PSMs"], width, label='#PSMs', 
                        color="black")
        rects2 = ax.bar(x + width/2, df["#Peptides"], width, label='#Peptides', 
                        color="red")
        # rects2 = ax.scatter(index, df["#Proteins"], label='#Proteins', color="blue")
        plt.ylabel("#Identifications")
        plt.xlabel("Fraction")
        # ax.set_xticklabels(df.index)
        # plt.xticks(x = list(index), label=df.index, rotation=90)
        # plt.xlim([-0.5, len(index) - 0.5])
        plt.legend()
        plt.show()


        fig, ax = plt.subplots(dpi=120)
        rects1 = ax.scatter(index, df["#CumulativePeptides"], 
                            label='#Peptides (cumul.)', color="black")
        rects2 = ax.scatter(index, df["#CumulativeProteins"], 
                            label='#Proteins (cumul.)', color="red")
        plt.ylabel("#Identifications")
        # plt.xlabel("Fraction")
        # plt.xlim([-0.5, len(index) - 0.5])
        # # plt.ylim([-0.5, plt.ylim()[1]])
        # plt.xticks(x = list(index), label=df.index, rotation=90)
        # labels = [item.get_text() for item in ax.get_xticklabels()]
        # ax.set_xticklabels(df.index)
        plt.legend()
        plt.show()
        
        # ["#PSMs", "#Peptides","#Proteins",
        # "#CumulativePeptides","#CumulativeProteins"]
    
    
    def plot_distributions_histogram(self, data, feature):
        plt.figure(dpi=120)
        sns.reset_orig()
        n, bins, patches = plt.hist(data, 50, facecolor='g', alpha=0.75)
        plt.xlabel(feature)
        plt.ylabel("Count")
        # plt.title('Histogram of {}'.format(feature))
        plt.grid(True)
        plt.show()
        
        
        
        


        

        
        
        
        
        
        
        
        
        
        
        
            
    def _setaxeslabels(self,axes, title, xlabel, ylabel):
        """  """
        axes.set_title(title)
        axes.set_xlabel(xlabel)   
        axes.set_ylabel(ylabel)
        
        
    def simple_plot(self, x_values, y_values):
        plt.figure()
        plt.plot(x_values, y_values)
        plt.show()
        
    def _save_figure(self, axes, file_name):
        figure = axes.get_figure()
        figure.savefig("{}.png".format(file_name), dpi=400, format='png')
        
    def _save_df(self, df, file_name):
        outfile = open(file_name+".tsv",'w')
        df.to_csv(outfile,sep='\t')
        outfile.close() 
        
    def _save_dict(self, d, file_name):
        data = "".join(["{},{}\n".format(str(key), str(value)) \
                                                for key, value in d.items()])
        with open(file_name, 'w') as f:
            print(data, file=f)
            
    def _save_list_list(self, l, file_name):
        #l =  zip(*l)
        with open(file_name, "w", newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerows(l)
        
    def create_plots_from_multiple_list(self, x, l, title="", xlabel="", ylabel=""):
        av, med, sem, sd = self.MyTools.calculate_average_and_sd_from_multiple_lists(l)
        self.createBarplot(x, av, sem, title, title, xlabel, ylabel)
        self.createBoxplot(l, showfliers=False, title=title, 
                      name=title, xlabel=xlabel, ylabel=ylabel)
        
    
    def createHeatmap(self, df, title: str="", 
                      name: str="Unknown",  
                      xlabel: str="", ylabel: str="", cmap: str="Greens", 
                      annot: bool =False, annot_kws: dict={"size": 10}, 
                      fmt: str="d", mask: bool=False):
        """  """
        if df.sum().sum()==0:
            print (title, "does not contain data")
            return None
        plt.figure()
        axes = sns.heatmap(df, cmap=cmap, annot=annot, annot_kws=annot_kws, 
                           mask=mask, fmt=fmt) 
        self._setaxeslabels(axes, title, xlabel, ylabel)
        file_name = "{}\\{}_heatmap".format(self.path, name) 
        self._save_figure(axes, file_name)
        self._save_df(df, file_name)
        plt.show()

    def createBarplot(self, x, y, yerr, title: str="", 
                      name: str="Unknown",  xlabel: str="", ylabel: str=""):
        plt.figure()
        plt.bar(x, y, yerr=yerr, align='center')
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        file_name = "{}\\{}_barplot".format(self.path, name) 
        plt.savefig("{}.png".format(file_name), dpi=400, format='png')
        d = {i:[j,k] for i,j,k in zip(x, y, yerr)}
        self._save_dict(d, "{}.txt".format(file_name))
        plt.show()
        
        
    def createBoxplot(self, data, showfliers, title: str="", 
                      name: str="Unknown",  xlabel: str="", ylabel: str=""):
        plt.figure()
        plt.boxplot(data, showfliers=showfliers) 
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        file_name = "{}\\{}_boxplot".format(self.path, name) 
        plt.savefig("{}.png".format(file_name), dpi=400, format='png')
        self._save_list_list(data, file_name)
        
    
    def create_bar_plot_from_dict(self, d, title: str="", 
                      name: str="Unknown",  xlabel: str="", ylabel: str=""):
        plt.figure()
        d = self.MyTools.sort_dictionary(d)
        plt.bar(range(len(d)), list(d.values()), align='center')
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(range(len(d)), list(d.keys()))
        file_name = "{}\\{}_barplot".format(self.path, name) 
        plt.savefig("{}.png".format(file_name), dpi=400, format='png')
        self._save_dict(d, "{}.txt".format(file_name))
        plt.show()
                

    def createHistogram(self,y ,x=False, bins=None, title="", 
                        name='Unknwon',   xlabel="",
                        ylabel="Anzahl", save_data = False):
        """ """         
        
        if len(y)==0:
            print (title, "does not contain data")
            return None
        
        plt.figure()
        if not x: x = range(len(y))
        axes = sns.histplot(y, bins=bins)
        self._setaxeslabels(axes, title, xlabel, ylabel)
        plt.show()
        
        file_name = "{}\\{}".format(self.path, name) 
        df = pd.DataFrame([h.get_height() for h in axes.patches], 
                          [h.get_x() for h in axes.patches])
        self._save_df(df, file_name+"_histogram_result")
        self._save_figure(axes, file_name)
        if save_data: self._save_df(pd.DataFrame(y), file_name+"_data")
            

       
       
    def createRelative2Dhist(self, list1, list2):
        """ Erzeugt relatives 2D Histogram ausgehend von klassischem 2DHist (hist= output from plt.hist2D) """
        plt.figure()
        hist = plt.hist2d(list1, list2, bins=(1000,500), 
                          norm=mpl.colors.LogNorm(), cmap=plt.cm.jet)
        new_array = []
        for line in hist[0]:
            max_list = max(line)
            new_line = [i/max_list for i in line] if max_list!=0 else line
            new_array.append(new_line)
        new_index = []
        z = 0
        for i in hist[1][:-1]:
            new_index.append(round((i+hist[1][z+1])/2))
            z+=1
        z2 = 0
        new_columns = []
        for i in hist[2][:-1]:
            new_columns.append((round((i+hist[2][z2+1])/2)))
            z2+=1
        df_hist2D_rel = pd.DataFrame(new_array,columns=new_columns,index=new_index)
        sns.heatmap(df_hist2D_rel)
        plt.show()
        return (df_hist2D_rel)