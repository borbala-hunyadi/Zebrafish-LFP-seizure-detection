# Zebrafish-LFP-seizure-detection

##A.Introduction 

1. Application
EEG signal processing

2. Purpose
For visualizing single-channel LFP signals, their time-frequency characteristics; and detect epileptic seizures in zebrafish larvae. 

3. Description
The algorithm (feature extraction and classifier) was optimized for detecting seizures in a specific genetic and chemical model of epilepsy as described in the reference paper below. For optimal performance on different zebrafish models and signals recorded under different conditions, the classifier needs to be retrained on data from a similar setting.

4. Citation
Borbála Hunyadi, Aleksandra Siekierska, Jo Sourbron, Daniëlle Copmans, Peter A.M. de Witte, Automated analysis of brain activity for seizure detection in zebrafish models of epilepsy, Journal of Neuroscience Methods, Volume 287, 1 August 2017, Pages 13-24, ISSN 0165-0270, https://doi.org/10.1016/j.jneumeth.2017.05.024.

5. Key words/tags
EEG, local field potential (LFP), seizure detection, epilepsy, animal model, machine learning

##B. Describing the Code 

Contents:
1. zebrafish.m: Matlab code
2. zebrafish.fig (the figure for the GUI, no need to open this one)
3. zebrafish.pdf: the documentation
4. example.txt: the software expects the EEG/LFP data in such a format.
5. example.mat: labels for the example.txt file. Labels are required for each 100ms segment of the data. 1 means that the segment contains seizure, 0 means it does not.
6. class_svm_chemical.mat and class_svm_genetic.mat: the classifiers
7. train_svm_zebrafish.m: if you want to train your own classifier, you need to use this function
8. train_svm_zebrafish.pdf: documentation for the training function
9. training_files.mat: and example training file for the train_svm_zebrafish function.

The code (zebrafish.m) implements a graphical user interface for visualizing single-channel EEG, detecting seizures and performing simple time-frequency analysis. 
Moreover, a separate function (train_svm_zebrafish.m) is provided to train a new classifier based on labelled data. For detailed description of both codes, please see the pdf documentation. 

Example use: 
[SVMModel]=train_svm_zebrafish(50,10,'chemical',cd,'training_files')


##C. Running the Code 

To run the GUI, simply type zebrafish in the Matlab command line and press Enter. The GUI will initialize.
To run the 
  	

##D. Program Output 
The outputs of the program are displayed in the GUI window or written in user-specified csv files.

##E. Others
The software expects that the sampling rate of the EEG/LFP is 1000Hz. If your recordings are different, you will first need to resample the data.
