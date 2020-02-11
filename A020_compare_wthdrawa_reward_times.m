% compare withdrawal vs. reward time
%clear all; close all
anim = 2
lock = 'FrontReward' %'p_o' % p_i: poke-in % p_o: poke-out
trial_selection =  'inseq_outseq_exp'; %'inseq_outseq_exp'  'otseqcorr_otseqincorr' %'corr_incorr_excld_A'  %'inseq_outseq_exp' %%'inseq_outseq_con' %inseqcorr_inseqincorr
match = 'yes';
task = 'welltrained' %novel1, novel2, welltrained
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

% load behavioral and neural data
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load([chan_name '5.mat']) % load exmp chan
load('BehaviorMatrix.mat')

% define eventWindow/lock
eventWindow = [-0.5 0.5];
pokeInAlignedBehavMatrix_1 = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'FrontReward');
pokeInAlignedBehavMatrix_2 = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut');
trial_times_reward = [pokeInAlignedBehavMatrix_1(:).TrialLogVect];
trial_times_p_o    = [pokeInAlignedBehavMatrix_2(:).TrialLogVect];
inSeqLog     = [pokeInAlignedBehavMatrix_1.TranspositionDistance]==0;
if anim == 6
    otSeqLog = [pokeInAlignedBehavMatrix_1.ItemItemDistance]~=1;
else
    otSeqLog = [pokeInAlignedBehavMatrix_1.TranspositionDistance]~=0;
end
corrTrlLog   = [pokeInAlignedBehavMatrix_1.Performance]==1;
inSeqCorrLog = inSeqLog&corrTrlLog==1;
otSeqCorrLog = otSeqLog&corrTrlLog==1;
fs = 1000;
% inseq corr
trial_times_reward_inseq = trial_times_reward(:,inSeqCorrLog);
trial_times_p_o_inseq    = trial_times_p_o(:,inSeqCorrLog);
for i = 1:size(trial_times_reward_inseq,2)
    rew_idx = find(trial_times_reward_inseq(:,i));
    p_o_idx = find(trial_times_p_o_inseq(:,i));
    rew_times(i) = rew_idx(1);
    poke_times(i) = p_o_idx(1);
    delay_inseq(i) = (rew_idx(1) - p_o_idx(1))/1000;
end
mx_inseq = max(delay_inseq)
mn_inseq = min(delay_inseq)
stdv_inseq = std(delay_inseq)
avg_inseq  = mean(delay_inseq)
% figure;
% plot(rew_times/fs, 'm*')
% hold on
% plot(poke_times/fs, 'b*')


% outseq corr
trial_times_reward_otseq = trial_times_reward(:,otSeqCorrLog);
trial_times_p_o_otseq    = trial_times_p_o(:,otSeqCorrLog);
for i = 1:size(trial_times_reward_otseq,2)
    rew_idx = find(trial_times_reward_otseq(:,i));
    p_o_idx = find(trial_times_p_o_otseq(:,i));
    rew_times_ot(i) = rew_idx(1);
    poke_times_ot(i) = p_o_idx(1);
    delay_otseq(i) = (rew_times_ot(1) - poke_times_ot(1))/1000;
end
mx_otseq = max(delay_otseq)
mn_otseq = min(delay_otseq)
stdv_otseq = std(delay_otseq)
avg_otseq  = mean(delay_otseq)

anim_6_avg_inseq = avg_inseq
anim_6_std_inseq = stdv_inseq
anim_6_avg_otseq = avg_otseq
anim_6_std_otseq = stdv_otseq

%
figure;
hold on
bar_vector1 = [anim_6_avg_inseq anim_6_avg_otseq]
bar_vector2 = [anim_6_std_inseq anim_6_std_otseq]
bar(bar_vector1 )
errorbar(1:2,bar_vector1,bar_vector2, 'rx')
title(['reward time - poke out time for anim ' num2str(anim)] )
ylabel('time in sec')
set(gca, 'XTick', 1:2, 'XTickLabel', {'inseq +' 'otseq +'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
