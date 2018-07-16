function HMMtoText(HMM, name)
% Martijn Zeestraten, August 2015
% Function writes the data from HMM to the following files:
% - HSMM_<name>_priors.txt
% - HSMM_<name>_trans.txt
% - HSMM_<name>_mu.txt
% - HSMM_<name>_sigma.txt
% - HSMM_<name>_varnames.txt
% The supplied HMM should have the following fields:
% - Sigma         [nbVar x nbVar x nbStates] covariance matrices
% - Mu            [nbVar x nbStates] means
% - StatesPriors  [nbStates x 1]  State priors
% - Trans         [nbStates x nbStates] Transition matrix
% - varnames      {1 x nbStates} Cell of variable names

if isdir('./textModels/')==0
	mkdir('./textModels/');
end

dlmwrite(['./textModels/HMM_', name, '_sigma.txt'],HMM.Sigma, 'delimiter', ' ');
dlmwrite(['./textModels/HMM_', name, '_mu.txt'],HMM.Mu, 'delimiter', ' ');	
dlmwrite(['./textModels/HMM_', name, '_priors.txt'],HMM.StatesPriors', 'delimiter', ' ');
dlmwrite(['./textModels/HMM_', name, '_trans.txt'],HMM.Trans, 'delimiter', ' ');


% Write Varnames to file:
varnames = [];
for i = 1:length(HMM.varnames)-1;
	varnames = [varnames, HMM.varnames{i}, ' '];
end
fileID = fopen(['./textModels/HMM_', name, '_varnames.txt'],'w');
varnames = [varnames, HMM.varnames{end}, ' '];
fprintf(fileID,varnames);
fclose(fileID);



end