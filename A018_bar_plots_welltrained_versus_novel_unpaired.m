clear all;close all;clc
reps = 100000;
cd ('D:\Gattas\ephys_data_final\welltrained\group_plots')

%%
elec  = 4
% comparison 1
% welltrained
load('group_animal_individual_trials_for_stats_beta_nonmatched.mat')
adata = [cond1_trials{elec,1:4}];
clear cond1_trials 
% novel 2
load('group_animal_individual_trials_for_stats_beta_nonmatched_novel2.mat')
bdata = [cond1_trials{elec,1:4}];
p1 = permutation_unpaired(adata, bdata, reps)
clear cond1_trials
clear adata bdata

% comparison 2
% welltrained
load('group_animal_individual_trials_for_stats_beta_nonmatched.mat')
adata = [cond1_trials{elec,1:4}];
clear cond1_trials
% novel1
load('group_animal_individual_trials_for_stats_beta_nonmatched_novel1.mat')
bdata = [cond1_trials{elec,1:4}];
p2 = permutation_unpaired(adata, bdata, reps)
clear cond1_trials
clear adata bdata

% comparison 3
% novel 2 vs. 1
load('group_animal_individual_trials_for_stats_beta_nonmatched_novel2.mat')
adata = [cond1_trials{elec,1:4}];
clear cond1_trials
load('group_animal_individual_trials_for_stats_beta_nonmatched_novel1.mat')
bdata = [cond1_trials{elec,1:4}];
p3 = permutation_unpaired(adata, bdata, reps)


%%
[p_fdr, p_masked] = fdr([p1 p2 p3], 0.05)
                            % all poss p-values