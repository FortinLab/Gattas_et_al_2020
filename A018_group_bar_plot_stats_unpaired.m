cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
load(['group_animal_individual_trials_for_stats_'  band '_nonmatched_' task])

reps = 100000;
for elec = 1:4


clear InSeqCorr_OutSeqCorr InSeqCorr_InSeqInCorr InSeqCorr_OutSeqInCorr...
    OutSeqCorr_InSeqInCorr OutSeqCorr_OutSeqInCorr InSeqInCorr_OutSeqInCorr

% InSeq+ vs. OutSeq-
adata = [cond1_trials{elec,1:6}];
bdata = [cond2_trials{elec,1:6}];
InSeqCorr_OutSeqCorr = permutation_unpaired(adata, bdata, reps)

% InSeq+ vs. InSeq-
adata = [cond1_trials{elec,1:6}];
bdata = [cond3_trials{elec,1:6}];
InSeqCorr_InSeqInCorr = permutation_unpaired(adata, bdata, reps);

% InSeq+ vs. OutSeq-
adata =[cond1_trials{elec,1:6}];
bdata = [cond4_trials{elec,1:6}];
InSeqCorr_OutSeqInCorr = permutation_unpaired(adata, bdata, reps);

% OutSeq+ vs. InSeq-
adata =[cond2_trials{elec,1:6}];
bdata = [cond3_trials{elec,1:6}];
OutSeqCorr_InSeqInCorr = permutation_unpaired(adata, bdata, reps);

% OutSeq+ vs. OutSeq-
adata =[cond2_trials{elec,1:6}];
bdata = [cond4_trials{elec,1:6}];
OutSeqCorr_OutSeqInCorr = permutation_unpaired(adata, bdata, reps);

% InSeq - vs. Outseq -
adata =[cond3_trials{elec,1:6}];
bdata = [cond4_trials{elec,1:6}];
InSeqInCorr_OutSeqInCorr = permutation_unpaired(adata, bdata, reps);

p_values(elec,:) = [InSeqCorr_OutSeqCorr InSeqCorr_InSeqInCorr InSeqCorr_OutSeqInCorr OutSeqCorr_InSeqInCorr OutSeqCorr_OutSeqInCorr InSeqInCorr_OutSeqInCorr]


[p_fdr, p_masked] = fdr([InSeqCorr_OutSeqCorr InSeqCorr_InSeqInCorr InSeqCorr_OutSeqInCorr OutSeqCorr_InSeqInCorr OutSeqCorr_OutSeqInCorr InSeqInCorr_OutSeqInCorr], 0.05)
signif (elec,:) = p_masked
p_fdr_all (elec)= p_fdr
end