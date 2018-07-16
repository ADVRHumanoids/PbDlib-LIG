function HSMMtoText(HSMM, name)
% Martijn Zeestraten, August 2015
% Function writes the data from HSMM to the following files:
% - HSMM_<name>_priors.txt
% - HSMM_<name>_trans.txt
% - HSMM_<name>_durMu.txt
% - HSMM_<name>_durSigma.txt
% - HSMM_<name>_mu.txt
% - HSMM_<name>_sigma.txt
% - HSMM_<name>_varnames.txt
% The supplied HSMM should have the following fields:
% - Sigma         [nbVar x nbVar x nbStates] covariance matrices
% - Mu            [nbVar x nbStates] means
% - StatesPriors  [nbStates x 1]  State priors
% - Trans         [nbStates x nbStates] Transition matrix
% - Mu_Pd         [1 x nbStates]        Matrix of duration means
% - Sigma_Pd      [1 x 1 x nbStates]    Matrix of duration covariances
% - varnames      {1 x nbStates} Cell of variable names

if isdir('./textModels/')==0
	mkdir('./textModels/');
end

dlmwrite(['./textModels/HSMM_', name, '_sigma.txt'],HSMM.Sigma, 'delimiter', ' ');
dlmwrite(['./textModels/HSMM_', name, '_mu.txt'],HSMM.Mu, 'delimiter', ' ');	
dlmwrite(['./textModels/HSMM_', name, '_priors.txt'],HSMM.StatesPriors', 'delimiter', ' ');
dlmwrite(['./textModels/HSMM_', name, '_trans.txt'],HSMM.Trans, 'delimiter', ' ');
dlmwrite(['./textModels/HSMM_', name, '_durMu.txt'],HSMM.Mu_Pd, 'delimiter', ' ');
dlmwrite(['./textModels/HSMM_', name, '_durSigma.txt'],HSMM.Sigma_Pd, 'delimiter', ' ');

% Write Varnames to file:
varnames = [];
for i = 1:length(HSMM.varnames)-1;
	varnames = [varnames, HSMM.varnames{i}, ' '];
end
fileID = fopen(['./textModels/HSMM_', name, '_varnames.txt'],'w');
varnames = [varnames, HSMM.varnames{end}, ' '];
fprintf(fileID,varnames);
fclose(fileID);



end