clc
clear
addpath(genpath('/home/nxs113020/Downloads/MSR_Identity_dir/MSR_Identity_Toolkit_v1.0/code/'));

%% Step0: Opening MATLAB pool
nworkers = 12;
nworkers = min(nworkers, feature('NumCores'));
isopen = matlabpool('size')>0;
if ~isopen, matlabpool(nworkers); end

%% Step1: Training the UBM
ubm_location = 'Aging/mf_UBM';
tmp = load(ubm_location);
ubm = tmp.ubm;
nmix = length(ubm.w);
clear tmp;

%% Step2: Learning the total variability subspace from background data
tv_location = 'Aging/mf_T';
tmp = load(tv_location);
T = tmp.T;
tv_dim = size(T,1); 
clear tmp

%% Step3: Training the Gaussian PLDA model with development i-vectors
lda_dim = 200;
nphi    = 200;
niter   = 10;
dataList = 'lists/dev_list';
fid = fopen(dataList, 'rt');
C = textscan(fid, '%s %s');
fclose(fid);
feaFiles = C{1};
#dev_ivs = zeros(tv_dim, length(feaFiles));

#stats = cell(length(feaFiles), 1);
#parfor file = 1 : length(feaFiles),
#    [N, F] = compute_bw_stats(feaFiles{file}, ubm);
#    stats{file} = [N; F];
#end 

#parfor file = 1 : length(feaFiles),
#    dev_ivs(:, file) = extract_ivector(stats{file}, ubm, T);
#end

#save('dev_ivectors','dev_ivs')
#save('stats','stats')
#pause
% reduce the dimensionality with LDA
spk_labs = C{2};
V = lda(dev_ivs, spk_labs);
dev_ivs = V(:, 1 : lda_dim)' * dev_ivs;
%------------------------------------
plda = gplda_em(dev_ivs, spk_labs, nphi, niter);


%% Step4: Scoring the verification trials
fea_dir = '../features/';
fea_ext = '.htk';
fid = fopen('speaker_model_maps.lst', 'rt');
C = textscan(fid, '%s %s');
fclose(fid);
model_ids = unique(C{1}, 'stable');
model_files = C{2};
nspks = length(model_ids);
model_ivs1 = zeros(tv_dim, nspks);
model_ivs2 = model_ivs1;
parfor spk = 1 : nspks,
    ids = find(ismember(C{1}, model_ids{spk}));
    spk_files = model_files(ids);
    spk_files = cellfun(@(x) fullfile(fea_dir, [x, fea_ext]),...  %# Prepend path to files
                       spk_files, 'UniformOutput', false);
    N = 0; F = 0; 
    for ix = 1 : length(spk_files),
        [n, f] = compute_bw_stats(spk_files{ix}, ubm);
        N = N + n; F = f + F; 
        model_ivs1(:, spk) = model_ivs1(:, spk) + extract_ivector([n; f], ubm, T);
    end
    model_ivs2(:, spk) = extract_ivector([N; F]/length(spk_files), ubm, T); % stats averaging!
    model_ivs1(:, spk) = model_ivs1(:, spk)/length(spk_files); % i-vector averaging!
end
trial_list = 'trials.lst';
fid = fopen(trial_list, 'rt');
C = textscan(fid, '%s %s %s');
fclose(fid);
[model_ids, ~, Kmodel] = unique(C{1}, 'stable'); % check if the order is the same as above!
[test_files, ~, Ktest] = unique(C{2}, 'stable');
test_files = cellfun(@(x) fullfile(fea_dir, [x, fea_ext]),...  %# Prepend path to files
                       test_files, 'UniformOutput', false);
test_ivs = zeros(tv_dim, length(test_files));
parfor tst = 1 : length(test_files),
    [N, F] = compute_bw_stats(test_files{tst}, ubm);
    test_ivs(:, tst) = extract_ivector([N; F], ubm, T);
end
% reduce the dimensionality with LDA
model_ivs1 = V(:, 1 : lda_dim)' * model_ivs1;
model_ivs2 = V(:, 1 : lda_dim)' * model_ivs2;
test_ivs = V(:, 1 : lda_dim)' * test_ivs;
%------------------------------------
scores1 = score_gplda_trials(plda, model_ivs1, test_ivs);
linearInd =sub2ind([nspks, length(test_files)], Kmodel, Ktest);
scores1 = scores1(linearInd); % select the valid trials

scores2 = score_gplda_trials(plda, model_ivs2, test_ivs);
scores2 = scores2(linearInd); % select the valid trials

%% Step5: Computing the EER and plotting the DET curve
labels = C{3};
eer1 = compute_eer(scores1(linearInd), labels, true); % IV averaging
hold on
eer2 = compute_eer(scores2(linearInd), labels, true); % stats averaging
