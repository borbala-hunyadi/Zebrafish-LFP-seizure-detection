# Zebrafish-LFP-seizure-detection

## A.Introduction 

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

## B. Describing the Code 

The code (code/seizure-detection/zebrafish.m) implements a graphical user interface for visualizing single-channel EEG, detecting seizures and performing simple time-frequency analysis. 
Moreover, a separate function (code/classifier-training/train_svm_zebrafish.m) is provided to train a new classifier based on labelled data. For detailed description of both codes, please see the pdf documentation within each folder (zebrafish.pdf and train_svm_zebrafish.pdf). 

# code/seizure-detection:
1. zebrafish.m: Matlab code implementing the detector GUI
2. zebrafish.fig: figure called by the GUI
3. zebrafish.pdf: detailed documentation
4. class_svm_chemical.mat and class_svm_genetic.mat: the classifiers used by the GUI
# code/classifier-training:
1. train_svm_zebrafish.m: Matlab code to train a new classifier if needed, instead of using the included model (class_svm_chemical.mat and class_svm_genetic.mat)
2. train_svm_zebrafish.pdf: documentation for the training function
3. training_files.mat: and example training file for the train_svm_zebrafish function.
# data:
1. example.txt: example data, the seizure-deteciton GUI software expects the EEG/LFP data in such a format.
2. example.mat: labels for the example.txt file. Labels are created for each 100ms segment of the data. 1 means that the segment contains seizure, 0 means it does not.


## C. Running the Code 

To run the GUI, simply type [>> zebrafish] in the Matlab command line and press Enter. The GUI will initialize.
An example to run the classidier training: [>> [SVMModel]=train_svm_zebrafish(50,10,'chemical',cd,'training_files')]
  	
## D. Program Output 
The outputs of the program are displayed in the GUI window or written in user-specified csv files.

## E. Others
The software expects that the sampling rate of the EEG/LFP is 1000Hz. If your recordings are different, you will first need to resample the data.
