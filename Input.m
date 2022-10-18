%% Input file for main_Pcca_macrophages.m
% Specifies general run conditions and filepathes;
%% Implemented by Anna-Simone Frank (anna-simone.frank@uib.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cPCCA+ Input data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: maximum number of mlecules (same for each species)
N=100;   % State-space size
% Perform clustering by PCCA+
target=1;               %target eigenvalue: For Q, target=0; for P, target=1
kmin=3;                 % minimum number of clusters
kmax=3;                 % maximum number of clusters
flag=1;                 % flag=1: full optimization, flag=0: unscaled initial guess
solver='nelder-mead';   % solver for unconstrained optimization problem
%no display of eigenvalue solver
opts.disp=0;

% Choose a time lage for the transition probability matrix: 
tau=100;
X=sprintf('The choosen time lag for the transition probability matrix is tau = %d', tau);
disp(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transition-path-theory (TPT) analysis Input data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PercTH = 0.99999999;  % Threshold defining the percental support area for the phenotype sets; 
%% Note: (!) Start and end sets for TPT analysis have to be changed manually for the various cases in TPTCases.m file for bistability (!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify file of PCCA results to save results in
filePCCA='.\Results_Main_PCCA\PCCA_results_File.txt';

% Specify file of TPT results to save results in
fileTPTBi='.\Results_Main_PCCA\TPTbistabilityCases_File.txt';

fileTPTTriSc1='.\Results_Main_PCCA\TPTtristability_SubCase1_File.txt';
fileTPTTriSc2='.\Results_Main_PCCA\TPTtristability_SubCase2_File.txt';
fileTPTTriSc3='.\Results_Main_PCCA\TPTtristability_SubCase3_File.txt';
fileTPTTriSc4='.\Results_Main_PCCA\TPTtristability_SubCase4_File.txt';
fileTPTTriSc5='.\Results_Main_PCCA\TPTtristability_SubCase5_File.txt';
fileTPTTriSc6='.\Results_Main_PCCA\TPTtristability_SubCase6_File.txt';

% Specify filepath for saving of figures
filepath='.\Results_Main_PCCA';
