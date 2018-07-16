function TrajGMMtoText(TrajGMM, name)
% 
% Leonel Rozo, August 2015
% 
% Function writes the data from HSMM to the following files:
% - TrajGMM_<name>_priors.txt
% - TrajGMM_<name>_mu.txt
% - TrajGMM_<name>_sigma.txt
% - TrajGMM_<name>_varnames.txt
% - TrajGMM_<name>_info.txt
% The supplied TrajGMM should have the following fields:
% - Sigma         [nbVar x nbVar x nbStates] covariance matrices
% - Mu            [nbVar x nbStates] means
% - Priors				[nbStates x 1]  State priors
% - varnames      {1 x nbStates} Cell of variable names
% - nbVarPos			[1 x 1] Dimension of position data
% - nbDeriv				[1 x 1] Number of static and dynamic features 
% - dt						[1 x 1] Time step

if isdir('./textModels/')==0
	mkdir('./textModels/');
end

dlmwrite(['./textModels/TrajGMM_', name, '_sigma.txt'],TrajGMM.Sigma, ...
	'delimiter', ' ');
dlmwrite(['./textModels/TrajGMM_', name, '_mu.txt'],TrajGMM.Mu, ...
	'delimiter', ' ');	
dlmwrite(['./textModels/TrajGMM_', name, '_priors.txt'],TrajGMM.Priors,...
	'delimiter', ' ');
info = [TrajGMM.nbVarPos ; TrajGMM.nbDeriv ; TrajGMM.dt];
dlmwrite(['./textModels/TrajGMM_', name, '_info.txt'], info,...
	'delimiter', ' ');

% Write Varnames to file:
varnames = [];
for i = 1:length(TrajGMM.varnames)-1;
	varnames = [varnames, TrajGMM.varnames{i}, ' '];
end
fileID = fopen(['./textModels/TrajGMM_', name, '_varnames.txt'],'w');
varnames = [varnames, TrajGMM.varnames{end}, ' '];
fprintf(fileID,varnames);
fclose(fileID);

end