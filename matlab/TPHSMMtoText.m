function TPHSMMtoText(HSMM, name)
% Martijn Zeestraten, August 2015
% Function writes the data from TP-HSMM to the following files:
% - HSMM_<name>_priors.txt
% - HSMM_<name>_trans.txt
% - HSMM_<name>_durMu.txt
% - HSMM_<name>_durSigma.txt
% - HSMM_<name>_mu_P*.txt
% - HSMM_<name>_sigma_P*.txt
% - HSMM_<name>_varnames.txt
% The supplied TP-HSMM should have the following fields:
% - Sigma         [nbVar x nbVar x nbFrames x nbStates] covariance matrices
% - Mu            [nbVar x nbFrames x nbStates] means
% - StatesPriors  [nbStates x 1]  State priors
% - Trans         [nbStates x nbStates] Transition matrix
% - Mu_Pd         [1 x nbStates]        Matrix of duration means
% - Sigma_Pd      [1 x 1 x nbStates]    Matrix of duration covariances
% - varnames      {1 x nbStates} Cell of variable names


if isdir('./textModels/')==0
	mkdir('./textModels/');
end


% c++ loadGMM from matlab requires seperarate files for each Frame. The
% filename needs to have the format *Pi.txt where i is the index of the
% frame.
for i = 1 : HSMM.nbFrames
	filename =  ['./textModels/HSMM_', name, '_sigma_P', num2str(i),'.txt'];
	dlmwrite(filename,squeeze(HSMM.Sigma(:,:,i,:)), 'delimiter', ' ');
	
	% Saving Mu
	filename =  ['./textModels/HSMM_', name, '_mu_P', num2str(i),'.txt'];
	dlmwrite(filename,squeeze(HSMM.Mu(:,i,:)), 'delimiter', ' ');	
	
end

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