function model = loadHSMMfromText(prefix)
% Martijn Zeestraten, August 2015
% Function loads the data from HSMM defined by the following files:
% - HSMM_<name>_priors.txt
% - HSMM_<name>_trans.txt
% - HSMM_<name>_durMu.txt
% - HSMM_<name>_durSigma.txt
% - HSMM_<name>_mu.txt
% - HSMM_<name>_sigma.txt
% - HSMM_<name>_varnames.txt
% The outputted HSMM will have the following fields:
% - Sigma         [nbVar x nbVar x nbStates] covariance matrices
% - Mu            [nbVar x nbStates] means
% - StatesPriors  [nbStates x 1]  State priors
% - Trans         [nbStates x nbStates] Transition matrix
% - Mu_Pd         [1 x nbStates]        Matrix of duration means
% - Sigma_Pd      [1 x 1 x nbStates]    Matrix of duration covariances
% - varnames      {1 x nbStates} Cell of variable names

model.Sigma_Pd     = dlmread([prefix, '_durSigma.txt']);
model.Mu_Pd        = dlmread([prefix, '_durMu.txt']);
model.Trans        = dlmread([prefix, '_trans.txt']);
model.StatesPriors = dlmread([prefix, '_priors.txt']);
model.Mu           = dlmread([prefix, '_mu.txt']);
model.nbStates     = size(model.Mu,2);
model.nbVar        = size(model.Mu,1);
model.Sigma        = reshape(dlmread([prefix, '_sigma.txt']),[model.nbVar,model.nbVar,model.nbStates]);

end