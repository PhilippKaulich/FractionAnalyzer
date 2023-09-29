# PD Fraction Analyzer 

**Summary**: The tool analyzes fractionated samples (bottom-up and top-down proteomics) analyzed by ProteomeDiscoverer. 

**Input**: Exported PSMs from analysis by ProteomeDiscoverer.

**Output**: Various visualizations for the analysis and characterization of the fractionation.  



## Contents

[TOC]

## Introduction

PD Fraction Analyzer is a tool for the easy analysis of fractionation datasets analyzed by Proteome Discoverer for both bottom-up and top-down proteomics datasets. As input only the exported Peptide Spectrum Matches (PSMs) or Proteoform-Spectrum-Matches (PrSMs) and a fasta file containing the respecetive proteins is needed. The Graphical User Interface of PD Fraction Analyzer allows multiple ways to analyze and compare different fractionated datasets in respect to the identified P(r)SMs, peptides and proteins. 

The following analysis are possible:

- Distributions of different properties of the identified peptides/proteoforms (e.g., mass, GRAVY score, isoelectric point, Xcorr, retention time, *m/z*, isolation interference, etc.) and proteins (mass, GRAVY score, pI) in the individual fractions.  
- Number of Identifications (P(r)SMs, Peptides/Proteoforms, Protein Groups) for the individual fractions. 
- Overlap of identified Peptides/Proteoforms and Proteins for the individual fractions. 
- Number of unique peptides in number of fractions. 
- Number of P(r)SMs for each fraction of individual peptides/proteoforms. 
- Orthogonality Analysis of first and second dimension. 
- Comparison of two different data sets with respect to all of the above analyses. 



## Technical Description 



## Requirements

- Python 3.9.13 or higher

### Packages

- Bio==1.5.3
- matplotlib==3.5.1
- matplotlib_venn==0.11.7
- numpy==1.21.5
- pandas==1.4.2
- PyQt5==5.15.7
- pyteomics==4.5
- scipy==1.7.3
- seaborn==0.11.2
- statannot==0.2.3





## How to run the Tool

Installation and use of Anaconda Distribution and its build-in command line prompt is highly recommended. In case you don't use Anaconda, make sure all required packages are installed upfront.

````powershell
$ cd <PATH/TO/SOURCES>
$ python GUI.py
````

 

## Input Data

The requried input data are:

- Fasta-file that was used for Proteome Discoverer search 
- List of the identified P(r)SMs exported as txt-file from Proteoform Discoverer results



## Graphical User Interface
<img src="Various\GUI.jpg" style="zoom:90%" alt="GUI"/>



### Initialization

| Name            | Widget     | Description                                                  |
| --------------- | ---------- | ------------------------------------------------------------ |
| Load Fasta File | Button     | Opens the file dialog to select fasta file  (see [Input Data](##Input Data)) |
| Load PSM File   | Button     | Opens the file dialog to select P(r)SMs results from Proteome Discoverer exportet as txt file (see [Input Data](##Input Data)) |
| PSM File name   | Text field | Custom name that will be used for labeling the dataset in figures. |
| Save path       | Button     | Opens File Dialog to select path where all results will be saved. |
| Initialize      | Button     | All data is loaded and prepared for further analysis.        |



### Data Assigment

| Name                  | Widget      | Description                                                  |
| --------------------- | ----------- | ------------------------------------------------------------ |
| Table with file names | Text fields | Contains all raw files in which P(r)SMs have been identified. The fraction number must be specified manually. If several replicates were measured, the replicate can also be specified. |
| Save Assignment       | Button      | Saves the changes made in the table.                         |
| Delete Assignment     | Button      | Deletes all manually set assignments (fraction numbers and replicates) and resets them to the original settings. |



### Dataset Analysis

| Name                | Widget         | Description                                                  |
| ------------------- | -------------- | ------------------------------------------------------------ |
| Selection           | Drop down menu | Selection of the dataset, fraction, and replicate that should be analyzed. If fraction and replicate are not specified, all fractions and replicates are considered. |
| Analysis            | Drop down menu | Selection of the Analysis type that should be performed. If necessary, the upper drop down menu specifies the analysis. |
| Checkbox: Plot      | Checkbox       | If checked, plots are generated and displayed.               |
| Checkbox: Show Data | Checkbox       | If checked, the data points in the plots are displayed.      |
| Checkbox: Export    | Checkbox       | If checked, the result files are saved in the path specified in "Initialization" |



### Dataset Comparison

| Name                 | Widget         | Description                                                  |
| -------------------- | -------------- | ------------------------------------------------------------ |
| Comparison           | Drop down menu | Selection of the datasets that should be performed.          |
| Analysis             | Drop down menu | Selection of the Analysis type that should be performed. If necessary, the upper drop down menu specifies the analysis. |
| Checkbox: Statistics | Checkbox       | If checked, statistics is performed for the dataset comparison and shown as asterisks in the plots. |
| Checkbox: Show Data  | Checkbox       | If checked, the data points in the plots are displayed.      |
| Checkbox: Export     | Checkbox       | If checked, the result files are saved in the path specified in "Initialization" |






## Datasets for Testing

A demo dataset and demo fast file is provided in Dataset/.





## Output

### Dataset Analysis

Analysis of a data set involves analyzing the properties of the identified peptides/proteoforms or proteins for each fraction. This includes not only the number of identifications but also the overlap between fractions and all properties that Proteome Discoverer provides as output (such as molecular mass, delta M, Xcorr) and the GRAVY score and isoelectric point. Thereby, the distribution of a selected property for the different fractions is shown as a boxplot diagram. In addition, the number of unique identifications in a given number of fractions can be determined. In addition, the overall distribution of certain properties of all identifications can be displayed as a histogram. Furthermore, an orthogonality analysis modified for offline fractionation according to Gilar et al. (2005) is possible. 

<img src="Various\Results_Fraction_Identifications_AbsolutAndCummulative.png" style="zoom:100%" alt="GUI"/>

Figure: Number of identifications in the individual identifications. (A) Number of identfied P(r)SMs and Peptides/Proteoforms and (B) cumulative number of Peptides/Proteoforms and Protein groups identified.  



<img src="Various\Results_Fraction_Overlap_PeptideAndProtein.png" style="zoom:100%" alt="GUI"/>

Figure: Overlap of identified (A) peptides and (B) proteins in the individual fractions. 



<img src="Various\Results_Fraction_PeptideMass.png" style="zoom:100%" alt="GUI"/>

Figure: Mass distribution of all individual fractions. 

<img src="Various\Results_All_Mass.png" style="zoom:90%" alt="GUI"/>

Figure: Mass distribution of all identifications. 



<img src="Various\Results_Fraction_UniquePeptideInFractions.png" style="zoom:85%" alt="GUI"/>

Figure: Number of peptides identified in one to ten different fractions. 

<img src="Various\Results_Fraction_Elution_Peptides.png" style="zoom:85%" alt="GUI"/>

Figure: Number of PSMs identified for the peptide RLEGSTIMDKDQVAIPLDRK identified in the individual fractions. 



<img src="Various\Results_Fraction_Orthogonality.png" style="zoom:85%" alt="GUI"/>

Figure: Orthogonality analysis between first dimension (offline fractionation) and second dimension (online LC-MS/MS). Berechnung der Orthogonalität nach Gilar et al. (2005). 





### Comparison of Two Datasets



<img src="Various\Results_Comparison_Fraction_Identifications_Proteins.png" style="zoom:85%" alt="GUI"/>

Figure: Comparison of identified proteins in individual fraction of two different datasets. 





<img src="Various\Results_Comparison_Fraction_PeptideMass.png" style="zoom:85%" alt="GUI"/>

Figure: Comparison of the mass distribution in individual fraction of two different datasets. 



<img src="Various\Results_Comparison_All_Mass_AbsolutAndDensity.png" style="zoom:85%" alt="GUI"/>



Figure: Comparison of the mass distribution of all identifications of two different datasets. (A) Density plot, (B) absolut numbers



<img src="Various\Results_Comparison_Fraction_ElutionDifference.png" style="zoom:85%" alt="GUI"/>

Figure: Elution difference of peptides identified in two datasets. 



<img src="Various\Results_Comparison_Overlap_PeptideAndProtein_Venny.png" style="zoom:30%" alt="GUI"/>

Figure: Overlap of identified (A) peptides and (B) proteins in two different datasets. 





<img src="Various\Results_Comparison_2ndDimension_Orthogonality_both_Correlation.png" style="zoom:25%" alt="GUI"/>



Figure: Orthogonality of the LC-MS/MS dimension between two datasets. 



<img src="Various\Results_Comparison_MultipleDatasets_PeptideMass.png" style="zoom:85%" alt="GUI"/>

Figure: Comparison of the mass distribution in individual fraction of all loaded datasets. 







## Referencs

Gilar, M.; Olivova, P.; Daly, A. E.; Gebler, J. C. Orthogonality of Separation in Two-Dimensional Liquid Chromatography. Anal. Chem. 2005, 77 (19), 6426–6434.



## Changelog

Jan, 2022: First Release.

Jun, 2022: Comparison of two datasets possible. 

Dez, 2022: Graphical User Interface. 

Mar, 2023: Minor bugfixes. 



## How to cite 





## License 

 



