clear all; close all
anim  = 1
match = ''
lock  = 'p_o' %'p_o' % p_i: poke-in % p_o: poke-out
task  = 'welltrained' %novel1, novel2, welltrained
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

% load behavioral and neural data
[chan_name, chan_length ] = get_anim_info(anim, task);
chan_counter = length(chan_length);
load('chan_artifact_thresh.mat')
load('BehaviorMatrix.mat')
load([chan_name '4.mat']) % load exmp chan


% variables
fs         = 1000; % sampling rate

% define eventWindow/lock
if strcmp('p_i', lock)
    eventWindow = [-0.5 1];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn'); 
elseif strcmp('p_o', lock)
    eventWindow = [-0.5 0.5];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut');
end

% get inseq/outseq trials and correct trials
inSeqLog     = [pokeInAlignedBehavMatrix.TranspositionDistance]==0;

if anim == 6
    otSeqLog = [pokeInAlignedBehavMatrix.ItemItemDistance]~=1;
else
    otSeqLog = [pokeInAlignedBehavMatrix.TranspositionDistance]~=0;
end
odor_excld_A = [pokeInAlignedBehavMatrix.Odor]~=1;
corrTrlLog   = [pokeInAlignedBehavMatrix.Performance]==1;
times        = [pokeInAlignedBehavMatrix.PokeDuration];

% create notch filter
flt = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);

% get correctly identified inseq and outseq trials exluding A
times = [pokeInAlignedBehavMatrix.PokeDuration];
trial_num = 1:length(inSeqLog);
x = [1.1 1.2 1 1.1 1.1 .8];
response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);

% redefine response threshold
inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A&response_time==1;
otSeqCorrLog = otSeqLog==1 & corrTrlLog==1 & odor_excld_A==1 & response_time==0;

% get incorrectly identified inseq and outseq trials exclud A
inSeqWrongLog = inSeqLog==1 & corrTrlLog==0 & odor_excld_A==1 & response_time==0;
otSeqWrongLog = otSeqLog==1 & corrTrlLog==0 & odor_excld_A==1 & response_time==1;

% 2 conditions of interest set as trial_1 and trial_2
trial_log_1 = inSeqCorrLog;
trial_log_2 = otSeqCorrLog;


% match for trial number
trial_min = min([sum(trial_log_1) sum(trial_log_2)]);
if strcmp('yes',match)
if sum(trial_log_2)< sum(trial_log_1)
    trial_1_idx = find(trial_log_1);
    indices = randperm(length(trial_1_idx), trial_min);
    trial_log_1a = zeros(1,length(trial_log_2));
    trial_log_1a(trial_1_idx(indices)) = 1;
    trial_log_1 = trial_log_1a;
    trial_log_1 = logical (trial_log_1);
elseif  sum(trial_log_1)< sum(trial_log_2)
    trial_2_idx = find(trial_log_2);
    indices = randperm(length(trial_2_idx), trial_min);
    trial_log_2a = zeros(1,length(trial_log_1));
    trial_log_2a(trial_2_idx(indices)) = 1;
    trial_log_2 = trial_log_2a;
    trial_log_2 = logical(trial_log_2);
end
end

% wavelet params
dt = 1/fs;
NumVoices = 32;
a0 = 2^(1/NumVoices);
wavCenterFreq = 6/(2*pi);
minfreq = 3;
maxfreq = 250;
minscale = wavCenterFreq/(maxfreq*dt);
maxscale = wavCenterFreq/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;
freq = wavCenterFreq./(fs*scales.*dt);

% calculating continuous wavelet transform of exemplar to get param t
LFP   = ExtractTrialData_SM(pokeInAlignedBehavMatrix, statMatrix(:,2));
LFP_1 = cell2mat(LFP(trial_log_1))';

cwt = cwtft({LFP_1(1, :),dt},...
    'scales',scales,'wavelet','morl');
cwt_power_exemp = abs(cwt.cfs).^2;

% init var
abs_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_1), chan_counter);
abs_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_2), chan_counter);

trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);

%
tic
for chan = 1:chan_counter % loop thru chans
    
    % LFP data
    data = load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
    
    % extract time-locked LFP data
    LFP = ExtractTrialData_SM(pokeInAlignedBehavMatrix, filtfilt(flt,data.statMatrix(:,2)) );
    
    % extract trials' time-locked LFP data
    LFP_1 = cell2mat(LFP(trial_log_1))';
    LFP_2 = cell2mat(LFP(trial_log_2))';

    trial_log_1_artifact_temp= zeros(sum(trial_log_1),1);
    trial_log_2_artifact_temp= zeros(sum(trial_log_2),1);
    
    
    %%  inseq
    parfor trial = 1:size(LFP_1, 1)
        % if trail has artifact = nan
        idx_above_thresh = (LFP_1(trial, :)>chan_artifact_thresh(chan,2));
        idx_below_thresh = (LFP_1(trial, :)<chan_artifact_thresh(chan,1));
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_1_artifact_temp(trial)= sum(total_idx);
        
        % wavelet spectrogram
        cwt = cwtft({LFP_1(trial, :),dt},...
            'scales',scales,'wavelet','morl');
        cwt_power = abs(cwt.cfs).^2;
        if sum(total_idx)<1 % if trial is artifact free
            norm_freq = cwt_power;
        else
            norm_freq = nan(size(cwt_power));
        end
        abs_freq_acrs_chan_cond_1(:, :, trial, chan) =  norm_freq;
        
    end
    
    trial_log_1_artifact{chan}=  trial_log_1_artifact_temp';
    
    
    
    %%  otseq
    parfor trial = 1:size(LFP_2, 1)
        % if trail has artifact = nan
        idx_above_thresh = (LFP_2(trial, :)>chan_artifact_thresh(chan,2));
        idx_below_thresh = (LFP_2(trial, :)<chan_artifact_thresh(chan,1));
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_2_artifact_temp(trial)= sum(total_idx);
        
        % wavelet spectrogram
        cwt = cwtft({LFP_2(trial, :),dt},...
            'scales',scales,'wavelet','morl');
        cwt_power = abs(cwt.cfs).^2;
        if sum(total_idx)<1 % if trial is artifact free
            norm_freq = cwt_power; 
        else
            norm_freq = nan(size(cwt_power));
        end
        abs_freq_acrs_chan_cond_2(:, :, trial, chan) =  norm_freq;
    end
    
    trial_log_2_artifact{chan}=  trial_log_2_artifact_temp';
    
   
    
end
toc

%% num of clean trials
for a = 1:size(trial_log_1_artifact,1)
    chan_trials_log_1 (a)= sum(trial_log_1_artifact{a}==0);
end
total_tiral_log_1 = sum(chan_trials_log_1);

for a = 1:size(trial_log_2_artifact,1)
    chan_trials_log_2 (a)= sum(trial_log_2_artifact{a}==0);
end
total_trial_log_2 = sum(chan_trials_log_2);

% get chans along axis
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;
freq_range= freq>19 & freq<36;
band_value = '19-36 hz';

%% get all trials that are not nans
% condition 1, all 4 chans
abs_freq_acrs_chan_cond_1_matched_chan_1 = cell(1,size(abs_freq_acrs_chan_cond_1,3));
counter1 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_1,3)
        if ~isnan(abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(1)))
            counter1 = counter1+1;
            abs_freq_acrs_chan_cond_1_matched_chan_1{counter1} = abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(1));
        end
    end
abs_freq_acrs_chan_cond_1_matched_chan_2 = cell(1,size(abs_freq_acrs_chan_cond_1,3));
counter2 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_1,3)
        if ~isnan(abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(2)))
            counter2 = counter2+1;
            abs_freq_acrs_chan_cond_1_matched_chan_2{counter2} = abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(2));
        end
    end    
    
    
abs_freq_acrs_chan_cond_1_matched_chan_3 = cell(1,size(abs_freq_acrs_chan_cond_1,3));
counter3 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_1,3)
        if ~isnan(abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(3)))
            counter3 = counter3+1;
            abs_freq_acrs_chan_cond_1_matched_chan_3{counter3} = abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(3));
        end
    end    
    
abs_freq_acrs_chan_cond_1_matched_chan_4 = cell(1,size(abs_freq_acrs_chan_cond_1,3));
counter4 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_1,3)
        if ~isnan(abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(4)))
            counter4 = counter4+1;
            abs_freq_acrs_chan_cond_1_matched_chan_4{counter4} = abs_freq_acrs_chan_cond_1(:, :, a, direc_chan(4));
        end
    end    
    
% condition 2, all 4 chans 
abs_freq_acrs_chan_cond_2_matched_chan_1 = cell(1,size(abs_freq_acrs_chan_cond_2,3));
counter5 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_2,3)
        if ~isnan(abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(1)))
            counter5 = counter5+1;
            abs_freq_acrs_chan_cond_2_matched_chan_1{counter5} = abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(1));
        end
    end

abs_freq_acrs_chan_cond_2_matched_chan_2 = cell(1,size(abs_freq_acrs_chan_cond_2,3));
counter6 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_2,3)
        if ~isnan(abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(2)))
            counter6 = counter6+1;
            abs_freq_acrs_chan_cond_2_matched_chan_2{counter6} = abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(2));
        end
    end    
    
abs_freq_acrs_chan_cond_2_matched_chan_3 = cell(1,size(abs_freq_acrs_chan_cond_2,3));
counter7 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_2,3)
        if ~isnan(abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(3)))
            counter7 = counter7+1;
            abs_freq_acrs_chan_cond_2_matched_chan_3{counter7} = abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(3));
        end
    end    
    
abs_freq_acrs_chan_cond_2_matched_chan_4 = cell(1,size(abs_freq_acrs_chan_cond_2,3));
counter8 = 0;
    for a = 1:size(abs_freq_acrs_chan_cond_2,3) % loop thru trials
        if ~isnan(abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(4))) % save trial if no artifact (not nan)
            counter8 = counter8+1;
            abs_freq_acrs_chan_cond_2_matched_chan_4{counter8} = abs_freq_acrs_chan_cond_2(:, :, a, direc_chan(4));
        end
    end  
min_trial_after_artifact_rejection =  min([counter1 counter2 counter3 counter4 counter5 counter6 counter7 counter8]);

% avg power, an across 6 trials
abs_beta_freq_acrs_cond_1_matched = nan(min_trial_after_artifact_rejection, 1001, length(direc_chan));
abs_beta_freq_acrs_cond_2_matched = nan(min_trial_after_artifact_rejection, 1001, length(direc_chan));

for a = 1:min_trial_after_artifact_rejection
    % cond1
abs_beta_freq_acrs_cond_1_matched(a,:,1)  = nanmean(abs_freq_acrs_chan_cond_1_matched_chan_1{a}(freq_range,:),1);
abs_beta_freq_acrs_cond_1_matched(a,:,2)  = nanmean(abs_freq_acrs_chan_cond_1_matched_chan_2{a}(freq_range,:),1);
abs_beta_freq_acrs_cond_1_matched(a,:,3)  = nanmean(abs_freq_acrs_chan_cond_1_matched_chan_3{a}(freq_range,:),1);
abs_beta_freq_acrs_cond_1_matched(a,:,4)  = nanmean(abs_freq_acrs_chan_cond_1_matched_chan_4{a}(freq_range,:),1);
    % cond2
abs_beta_freq_acrs_cond_2_matched(a,:,1)  = nanmean(abs_freq_acrs_chan_cond_2_matched_chan_1{a}(freq_range,:),1);
abs_beta_freq_acrs_cond_2_matched(a,:,2)  = nanmean(abs_freq_acrs_chan_cond_2_matched_chan_2{a}(freq_range,:),1);
abs_beta_freq_acrs_cond_2_matched(a,:,3)  = nanmean(abs_freq_acrs_chan_cond_2_matched_chan_3{a}(freq_range,:),1);
abs_beta_freq_acrs_cond_2_matched(a,:,4)  = nanmean(abs_freq_acrs_chan_cond_2_matched_chan_4{a}(freq_range,:),1);
end

% start and end time idx
if sum(trial_log_2)==0 || min(times(trial_log_2))>.5
    start_time = 1;
else
    start_time = (.5- min(times(trial_log_2)))*fs;
end
end_time = 500;

%% keep trial values and error bars

Cond_1_trials_power = squeeze(mean(abs_beta_freq_acrs_cond_1_matched(:,start_time:end_time,:),2));
Cond_2_trials_power = squeeze(mean(abs_beta_freq_acrs_cond_2_matched(:,start_time:end_time,:),2));
cond_diff = Cond_1_trials_power - Cond_2_trials_power
cond_sum = Cond_1_trials_power + Cond_2_trials_power

weighted_diff = cond_diff./cond_sum
figure
bar(mean(weighted_diff,1))
hold on
errorbar(1:4, mean(weighted_diff,1), (std(weighted_diff, 0,1))/sqrt(size(weighted_diff,1)), 'rx')
ylim([-.2 .4])
p = permutation_paired(weighted_diff(:,1), weighted_diff(:,3), 10000)
 
 %%
% average across a fixed num of trials across condition and elect
abs_beta_freq_acrs_cond_1_matched_avg_trials = (squeeze(mean(abs_beta_freq_acrs_cond_1_matched,1)))';
abs_beta_freq_acrs_cond_2_matched_avg_trials = (squeeze(mean(abs_beta_freq_acrs_cond_2_matched,1)))';



abs_beta_pr_cond1 = mean(abs_beta_freq_acrs_cond_1_matched_avg_trials(:,start_time:end_time),2);
abs_beta_pr_cond2 = mean(abs_beta_freq_acrs_cond_2_matched_avg_trials(:,start_time:end_time),2);

diff_cond = [abs_beta_pr_cond1 - abs_beta_pr_cond2]
sum_cond  = [abs_beta_pr_cond1 + abs_beta_pr_cond2]
figure
bar(diff_cond./sum_cond)
ylim([-.2 .4])


% save data per anim
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
%cond_1_2_diff = cell(1,6)
load ('cond_1_2_diff_absolute_freq')
cond_1_2_diff{anim} = [abs_beta_pr_cond1 - abs_beta_pr_cond2];
save('cond_1_2_diff_absolute_freq', 'cond_1_2_diff')

 