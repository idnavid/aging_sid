addpath(genpath('/home/nxs113020/Downloads/MSR_Identity_dir/MSR_Identity_Toolkit_v1.0/code/'));
% 
% %% Step0: Opening MATLAB pool
% parpool;
% 
% nworkers = 12;
% 
% %% Step1: Training the UBM
% dataList = 'lists/ubm_male_nist0506.lst';
% nmix        = 1024;
% final_niter = 10;
% ds_factor   = 10;
% %ubm = gmm_em(dataList, nmix, final_niter, ds_factor, nworkers);
% 
% %save('models/ubm_male','ubm');
% 
% %% Step1: Training the UBM
% ubm_location = 'models/ubm_male';
% tmp = load(ubm_location);
% ubm = tmp.ubm;
% nmix = length(ubm.w);
% clear tmp;
% 
%% Step2: Learning the total variability subspace from background data
% tv_dim = 400; 
% niter  = 5;
% dataList = 'lists/ubm_male_nist0506.lst';
% fid = fopen(dataList, 'rt');
% C = textscan(fid, '%s');
% fclose(fid);
% feaFiles = C{1};
% stats = cell(length(feaFiles), 1);
% parfor file = 1 : length(feaFiles),
%     [N, F] = compute_bw_stats(feaFiles{file}, ubm);
%     stats{file} = [N; F];
% end
% T = train_tv_space(stats, ubm, tv_dim, niter, nworkers);
% 
%save('models/T_male','T');
%
% tv_location = 'models/T_male';
% tmp = load(tv_location);
% T = tmp.T;
% tv_dim = size(T,1); 
% clear tmp
% 
%% Step3: Training the Gaussian PLDA model with development i-vectors
lda_dim = 90;
nphi    = 90;
niter   = 10;
pldaLists = {'lists/plda_male_20anddown_htk.lst';'lists/plda_male_23andup_htk.lst'};
spk_labs_cell = cell(2,1);
dev_ivs_cell = cell(2,1);
for i = [1,2]
    dataList = pldaLists{i};
    fid = fopen(dataList, 'rt');
    C = textscan(fid, '%s %s');
    fclose(fid);
    feaFiles = C{1};
    dev_ivs = zeros(tv_dim, length(feaFiles));
    stats = cell(length(feaFiles), 1);
    parfor file = 1 : length(feaFiles),
        [N, F] = compute_bw_stats(feaFiles{file}, ubm);
        stats{file} = [N; F];
    end
    
    parfor file = 1 : length(feaFiles),
        dev_ivs(:, file) = extract_ivector(stats{file}, ubm, T);
    end
    
    % reduce the dimensionality with LDA
    spk_labs_cell{i} = C{2};
    V = lda(dev_ivs, spk_labs_cell{i});
    dev_ivs_cell{i} = V(:, 1 : lda_dim)' * dev_ivs;
end
%------------------------------------
%plda = gplda_em(dev_ivs, spk_labs, nphi, niter);

plda = gplda_em_weighted_likelihood(dev_ivs_cell{1}, dev_ivs_cell{2}, spk_labs_cell{1}, spk_labs_cell{2}, nphi, niter);

%% Step4: Scoring the verification trials
fea_dir = '/home/nxs113020/features/aging_features/';
fea_ext = '.htk';
fid = fopen('lists/speaker_models', 'rt');
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

trial_list = 'lists/trials_A';
fid = fopen(trial_list, 'rt');
C = textscan(fid, '%s %s %s');
fclose(fid);
[model_ids, ~, Kmodel] = unique(C{1}, 'stable'); % check if the order is the same as above!
[test_files, ~, Ktest] = unique(C{2}, 'stable');
test_files = cellfun(@(x) fullfile(fea_dir, [x, fea_ext]),...  %# Prepend path to files
                       test_files, 'UniformOutput', false);
test_ivs = zeros(tv_dim, length(test_files));
parfor tst = 1 : length(test_files)
    [N, F] = compute_bw_stats(test_files{tst}, ubm);
    test_ivs(:, tst) = extract_ivector([N; F], ubm, T);
end
% reduce the dimensionality with LDA
model_ivs1_lda = V(:, 1 : lda_dim)' * model_ivs1;
model_ivs2_lda = V(:, 1 : lda_dim)' * model_ivs2;
test_ivs_lda = V(:, 1 : lda_dim)' * test_ivs;
%------------------------------------

scores1 = score_gplda_trials(plda, model_ivs1_lda, test_ivs_lda);
linearInd = sub2ind([nspks, length(test_files)], Kmodel, Ktest);
scores1 = scores1(linearInd); % select the valid trials

scores2 = score_gplda_trials(plda, model_ivs2_lda, test_ivs_lda);
scores2 = scores2(linearInd); % select the valid trials

%% Step5: Computing the EER and plotting the DET curve
labels = C{3};
[eer1, dcf08, dcf10] = compute_eer(scores1(:), labels, true) % IV averaging
hold on
[eer2, dcf08, dcf10] = compute_eer(scores2(:), labels, true); % stats averaging
