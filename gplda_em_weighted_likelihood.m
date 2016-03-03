function plda = gplda_em_weighted_likelihood(data1, data2, spk_labs1, spk_labs2, nphi, niter)
% trains a PLDA model given training data using EM algorithm. Speaker
% labels are given in spk_labs. nphi determines the dimensionality of the
% Eigenvoice subspace, and niter specifies the number of EM iteration.
%
% Technically, assuming a factor analysis (FA) model of the form:
%
%           x = m + Phi . y + e
%
% for i-vectors, the code computes the maximum likelihood estimate (MLE)of
% the factor loading matrix Phi (aka the Eigenvoice subspace). Here, x is
% the i-vector, m is the mean of training i-vectors, and y~N(0,I) is a
% vector of latent factors. The residual noise term e explains the
% variability not captured through the latent variables.
%
% Inputs:
%   data            : input data matrix, one observation per column
%   spk_labs        : class labels for observations in data matrix 
%   nphi            : dimensionality of the Eigenvoice subspace 
%   niter           : number of EM iterations for learning PLDA model
%
% Outputs:
%    plda           : a structure containing the PLDA hyperparameters
%					  as well as the mean of development data and a whitening 
%					  transform:(plda.Phi: Eigenvoice matrix, plda.Sigma: covariance
%					  matrix of the residual noise, plda.M: mean, plda.W: whitening transform)
%
% References:
%   [1] S.J.D. Prince and J.H. Elder, "Probabilistic linear discriminant analysis
%       for inferences about identity," in Proc. IEEE ICCV, Rio de Janeiro, Brazil,
%       Oct. 2007.
%   [2] D. Garcia-Romero and C.Y. Espy-Wilson, "Analysis of i-vector length 
%       normalization in speaker recognition systems," in Proc. INTERSPEECH,
%       Florence, Italy, Aug. 2011, pp. 249-252.
%   [3] P. Kenny, "Bayesian speaker verification with heavy-tailed priors," 
%       in Proc. Odyssey, The Speaker and Language Recognition Workshop, Brno, 
%       Czech Republic, Jun. 2010.
%
%
% Omid Sadjadi <s.omid.sadjadi@gmail.com>
% Microsoft Research, Conversational Systems Research Center

[ndim, nobs1] = size(data1);
[~, nobs2] = size(data2);

if ( nobs1 ~= length(spk_labs1) || nobs2 ~= length(spk_labs2) ),
	error('oh dear! number of data samples should match the number of labels!');
end

% make sure the labels are sorted
[spk_labs1, I] = sort(spk_labs1);
data1 = data1(:, I);
[~, ia, ic] = unique(spk_labs1, 'stable');
spk_counts1 = histc(ic, 1 : numel(ia)); % # sessions per speaker

[spk_labs2, I] = sort(spk_labs2);
data2 = data2(:, I);
[~, ia, ic] = unique(spk_labs2, 'stable');
spk_counts2 = histc(ic, 1 : numel(ia)); % # sessions per speaker

data = [data1,data2];
M = mean(data, 2);
data = bsxfun(@minus, data, M); % centering the data
data = length_norm(data); % normalizing the length
W1   = calc_white_mat(cov(data'));
data1 = W1' * data1; % whitening the data
data2 = W1' * data2;

fprintf('\n\nRandomly initializing the PLDA hyperparameters ...\n\n');
% Initialize the parameters randomly
[s1, s2] = RandStream.create('mrg32k3a', 'NumStreams', 2);
Sigma    = 100 * randn(s1, ndim); % covariance matrix of the residual term
Phi = randn(s2, ndim, nphi); % factor loading matrix (Eignevoice matrix)
    
Phi = bsxfun(@minus, Phi, mean(Phi, 2));
W2   = calc_white_mat(Phi' * Phi);
Phi = Phi * W2; % orthogonalize Eigenvoices (columns)

fprintf('Re-estimating the Eigenvoice subspace with %d factors ...\n', nphi);
for iter = 1 : niter
    fprintf('EM iter#: %d \t', iter);
    tim = tic;
    % expectation
    [Ey1, Eyy1] = expectation_plda(data1, Phi, Sigma, spk_counts1);
    [Ey2, Eyy2] = expectation_plda(data2, Phi, Sigma, spk_counts2);
    % maximization
    [Phi, Sigma] = maximization_plda(data1, data2, Ey1, Eyy1, Ey2, Eyy2);
    tim = toc(tim);
    fprintf('[elaps = %.2f s]\n', tim);
end

plda.Phi   = Phi;
plda.Sigma = Sigma;
plda.W     = W1;
plda.M     = M;

function [Ey, Eyy] = expectation_plda(data, Phi, Sigma, spk_counts)
% computes the posterior mean and covariance of the factors
nphi     = size(Phi, 2);
nsamples = size(data, 2);
nspks    = size(spk_counts, 1);

Ey  = zeros(nphi, nsamples);
Eyy = zeros(nphi);

% initialize common terms to save computations
uniqFreqs  	  = unique(spk_counts);
nuniq 		  = size(uniqFreqs, 1);
invTerms      = cell(nuniq, 1);
invTerms(:)   = {zeros(nphi)};
PhiT_invS_Phi = ( Phi'/Sigma ) * Phi;
I = eye(nphi);
for ix = 1 : nuniq
    nPhiT_invS_Phi = uniqFreqs(ix) * PhiT_invS_Phi;
    Cyy =  pinv(I + nPhiT_invS_Phi);
    invTerms{ix} = Cyy;
end

data = Sigma\data;
cnt  = 1;
for spk = 1 : nspks
    nsessions = spk_counts(spk);
    % Speaker indices
    idx = cnt : ( cnt - 1 ) + spk_counts(spk);
    cnt  = cnt + spk_counts(spk);
    Data = data(:, idx);
    PhiT_invS_y = sum(Phi' * Data, 2);
    Cyy = invTerms{ uniqFreqs == nsessions };
    Ey_spk  = Cyy * PhiT_invS_y;
    Eyy_spk = Cyy + Ey_spk * Ey_spk';
    Eyy     = Eyy + nsessions * Eyy_spk;
    Ey(:, idx) = repmat(Ey_spk, 1, nsessions);
end

function [Phi, Sigma] = maximization_plda(data1, data2, Ey1, Eyy1, Ey2, Eyy2)
% ML re-estimation of the Eignevoice subspace and the covariance of the
% residual noise (full).
alpha = 0.5;
nsamples1 = size(data1, 2);
nsamples2 = size(data2, 2);
nsamples = nsamples1 + nsamples2;

alpha1 = alpha/nsamples1;
alpha2 = (1-alpha)/nsamples2;

Data1_sqr = data1 * data1';
Data2_sqr = data2 * data2';
Phi      = (data1*Ey1' + data2*Ey2') / (nsamples*(alpha1*Eyy1 + alpha2*Eyy1));
Sigma    = (alpha1*Data1_sqr + alpha2*Data2_sqr) - Phi*(alpha1*Ey1*data1' + alpha2*Ey2*data2');

function W = calc_white_mat(X)
% calculates the whitening transformation for cov matrix X
[~, D, V] = svd(X);
W = V * diag(sparse(1./( sqrt(diag(D)) + 1e-10 )));
    
