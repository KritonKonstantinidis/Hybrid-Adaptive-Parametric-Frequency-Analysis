%%%% This script generates the spectrogram of an EEG dataset.

%%%% If the user wishes to monitor the convergence of the EM routine, 
%%%% please uncomment lines 100-114 in the HybridEM routine.
%%%% If the user wishes to monitor the complete data log-likelihood,
%%%% please use HybridEMwithQ in line 21 of the routine ModelSelection and
%%%% uncomment lines 65-75 in ModelSelection
%%%% If the user wishes to monitor information criteria, please uncomment
%%%% lines 51-103 in the routine ModelSelection.

clear all
close all
clc

fs = 250; % Sampling frequency [Hz] (index/sec)

y=load(fullfile('Datasets','SED974_lite.mat'));
y=y.data;
y=y(1,:);
dataEEG=y(306*60*250:316*60*250)';


HybridAlgorithm(dataEEG,'pmin',11,'pmax',15,'maxIter',50,'minibatch_size',1250,'online',0);
