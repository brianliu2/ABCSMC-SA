%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Code for ABCSMC+SA algorithm, all algorithmic settings are defined in
%     paraSettingHS which specifies the number of sample and
%     research curve for eFAST, dimensional information about system of
%     interest, prior distribution for generating particles for ABC and
%     relevant tolerances for ABC
%   
%     Indicator for implying which parameter is stiff or sloppy is given in
%     stiffInd, and it will be passed to sub-function abcsmcRepressilator.
%     In this function, perturbation will be added to particles based on
%     the stiff/sloppy indicators.
%
%     Author: Xin Liu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 
%   Load the algorithmic settings and the model settings for individual
%   biological system. Three settings including the heat shock,
%   repressilator and time-delay oscillatory are provided
paraSettingHS;

%   Define the global variables which will be use 
global sumSi;
global stiffInd;
global paraDim;
global sensLab;

%   Run the eFAST to get the sensitivity for guiding the selective
%   computation allocation

run SAMain_HS


%   Take the sensitivity for one of states in system as the clue for grouping stiff and
%   sloppy

stiffInd = find(sumSi(:,:,1)>0.15);
sensLab = cell(paraDim,1);
for i = 1:length(stiffInd)
    sensLab{stiffInd(i)} = ['stiff'];
end

DimInd = 1:paraDim;
DimInd([stiffInd]) = [];

for i = 1:length(DimInd)
    sensLab{DimInd(i)} = ['sloppy'];
end

%   Pass the indices for ABC-SMC algorithm by defining the indices as
%   global variable
run abcsmcHS





