addpath(genpath('/home/nxs113020/Downloads/MSR_Identity_dir/MSR_Identity_Toolkit_v1.0/code/'));
%%
ubm_location = 'models/ubm_male';
tmp = load(ubm_location);
ubm = tmp.ubm;
nmix = length(ubm.w);
clear tmp;
%%
tv_location = 'models/T_male';
tmp = load(tv_location);
T = tmp.T;
tv_dim = size(T,1); 
clear tmp
%%
deviv_location = 'models/dev_ivectors';
tmp = load(deviv_location);
dev_ivs = tmp.dev_ivs;
tv_dim = size(T,1); 
clear tmp
%%
lda_dim = 100;
nphi    = 100;
niter   = 10;
dataList = 'lists/plda_swb_20lower_23higher_htk.lst';
fid = fopen(dataList, 'rt');
C = textscan(fid, '%s %s');
fclose(fid);
feaFiles = C{1};
spk_labs = C{2};
V = lda(dev_ivs, spk_labs);
dev_ivs = V(:, 1 : lda_dim)' * dev_ivs;
plda = gplda_em(dev_ivs, spk_labs, nphi, niter);
%%
fea_dir = '/home/nxs113020/features/marp_male_60sec_chunks_postVAD/';
fea_ext = '.htk';
fid = fopen('lists/MARP_lists/speaker_models_1_19', 'rt');
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
%%
fea_dir = '/home/nxs113020/features/marp_male_60sec_chunks_postVAD/';
trial_list = 'lists/MARP_lists/m_trials_1_19.lst';
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
hold on 
c = 'r';
[eer2, dcf08, dcf10] = compute_eer(scores2, labels, true,c); % stats averaging
