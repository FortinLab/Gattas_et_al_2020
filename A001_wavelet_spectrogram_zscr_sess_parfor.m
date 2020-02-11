clear all; %close all
anim = % 1:5
lock = 'p_i' %'p_o' % p_i: poke-in % p_o: poke-out
trial_selection = 'inseq_outseq_exp'; %'inseq_outseq_exp'  'otseqcorr_otseqincorr' %'corr_incorr_excld_A'  %'inseq_outseq_exp' %%'inseq_outseq_con' %inseqcorr_inseqincorr
match = '';
plot_type = '';
task = 'welltrained' %novel1, novel2, welltrained
per = 0.5; % caxis variable for plotting
scale = 'db'
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

% load behavioral and neural data
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load([chan_name '5.mat']) % load exmp chan
load('baseline_info_wavelet_32num_10logdb_3hz_250hz_notched_artifact_reject.mat')
load('chan_artifact_thresh.mat')
load('BehaviorMatrix.mat')

% create notch filter
fs  = 1000;
flt = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);


% define eventWindow/lock
if strcmp('p_i', lock)
    eventWindow = [-0.5 1];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn');
    
elseif strcmp('p_o', lock)
    eventWindow = [-0.5 0.5];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut');
elseif strcmp('FrontReward', lock)
    eventWindow = [-0.5 0.5];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'FrontReward');
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

% get correctly identified inseq and outseq trials exluding A
times = [pokeInAlignedBehavMatrix.PokeDuration];
trial_num = 1:length(inSeqLog);

% anim 3 has a 47 min outlier trials = nan
if anim==3 && strcmp('novel2', task) || anim==4 && strcmp('novel2', task)  
    times(find(times==max(times))) = nan;
end
[ cutoff ] = get_response_time_cutoff( anim, task, times,inSeqLog ,otSeqLog, trial_num);
 response_time= [pokeInAlignedBehavMatrix.PokeDuration]>cutoff;

%% redefine response threshold

inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A&response_time==1;
otSeqCorrLog = otSeqLog==1 &corrTrlLog==1 & odor_excld_A==1 & response_time==0;

% get incorrectly identified inseq and outseq trials exclud A
inSeqWrongLog = inSeqLog==1 & corrTrlLog==0 & odor_excld_A==1 & response_time==0;
otSeqWrongLog = otSeqLog==1 & corrTrlLog==0 & odor_excld_A==1 & response_time==1;

% corr vs. incorr exld A
corrTrlLog_exld_A   = inSeqCorrLog+otSeqCorrLog;
incorrTrlLog_exld_A = inSeqWrongLog+otSeqWrongLog;


% 2 conditions of interest set as trial_1 and trial_2
if strcmp('corr_incorr_excld_A', trial_selection)
    trial_log_1 = logical(corrTrlLog_exld_A);
    trial_log_2 = logical(incorrTrlLog_exld_A);
elseif strcmp('corr_incorr_incld_A', trial_selection)
    trial_log_1 = corrTrlLog;
    trial_log_2 = ~corrTrlLog;
    fignam = '  corr vs.  incorr '
elseif strcmp('inseq_outseq_exp', trial_selection)
    trial_log_1 = inSeqCorrLog;
    trial_log_2 = otSeqCorrLog;
    fignam = ' inseq corr vs. outseq corr '
elseif strcmp('inseq_outseq_con', trial_selection)
    trial_log_1 = inSeqWrongLog;
    trial_log_2 = otSeqWrongLog;
    fignam = ' inseq incorr vs. outseq incorr '
elseif strcmp('inseqcorr_inseqincorr', trial_selection)
    trial_log_1 = inSeqCorrLog;
    trial_log_2 = inSeqWrongLog;
    fignam = ' inseq corr vs. inseq incorr '
elseif strcmp('otseqcorr_otseqincorr', trial_selection)
    trial_log_1 = otSeqCorrLog;
    trial_log_2 = otSeqWrongLog;
    fignam = ' otseq corr vs. otseq incorr '
elseif strcmp('inseqcorr_otseqincorr', trial_selection)
    trial_log_1 = inSeqCorrLog;
    trial_log_2 = otSeqWrongLog;
    fignam = ' inseq corr vs. otseq incorr '
end


%% control for trial num
trial_min = min(sum(trial_log_1), sum(trial_log_2));
if strcmp('yes', match)
    if trial_min>4
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
    else
        trls = 'non-matched'
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
norm_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_1), chan_counter);
norm_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_2), chan_counter);
absolute_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_1), chan_counter);
absolute_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_2), chan_counter);
trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);

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
        cwt_power = 20*log10(abs(cwt.cfs).^2);
    
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(cwt_power));
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/chan_powr_std(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(cwt_power));
        end
        norm_freq_acrs_chan_cond_1(:, :, trial, chan) =  norm_freq;
        absolute_freq_acrs_chan_cond_1(:, :, trial, chan) =  cwt_power;
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

        cwt_power = 20*log10(abs(cwt.cfs).^2);
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(cwt_power));
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/chan_powr_std(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(cwt_power));
        end
        norm_freq_acrs_chan_cond_2(:, :, trial, chan) =  norm_freq;
        absolute_freq_acrs_chan_cond_2(:, :, trial, chan) =  cwt_power;

    end
    
    trial_log_2_artifact{chan}=  trial_log_2_artifact_temp';
    
end
toc

% num of clean trials
for a = 1:size(trial_log_1_artifact,1)
    chan_trials_log_1 (a)= sum(trial_log_1_artifact{a}==0);
end
total_tiral_log_1 = sum(chan_trials_log_1);
for a = 1:size(trial_log_2_artifact,1)
    chan_trials_log_2 (a)= sum(trial_log_2_artifact{a}==0);
    
end
total_tiral_log_2 = sum(chan_trials_log_2);


% initialize matrix for all frequencies
mn_acrs_trials_inseq = zeros(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2), chan_counter);
mn_acrs_trials_outseq = zeros(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2), chan_counter);
absolute_mn_acrs_trials_inseq = zeros(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2), chan_counter);
absolute_mn_acrs_trials_outseq = zeros(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2), chan_counter);

% average across trials within chans
for elec = 1:size(norm_freq_acrs_chan_cond_1,4) %loop thru chans
    mn_acrs_trials_inseq (:, :, elec)  = nanmean(norm_freq_acrs_chan_cond_1(:, :, :, elec), 3);
    mn_acrs_trials_outseq (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_2(:, :, :, elec), 3);
    absolute_mn_acrs_trials_inseq(:, :, elec)  =  nanmean(absolute_freq_acrs_chan_cond_1(:, :, :, elec), 3);
    absolute_mn_acrs_trials_outseq(:, :, elec) =  nanmean(absolute_freq_acrs_chan_cond_2(:, :, :, elec), 3);
end

% mean acrs chans
indiv_freq_inseq_trls_chans = nanmean(mn_acrs_trials_inseq,3);
indiv_freq_otseq_trls_chans = nanmean(mn_acrs_trials_outseq,3);
absolute_indiv_freq_inseq_trls_chans = nanmean(absolute_mn_acrs_trials_inseq,3);
absolute_indiv_freq_otseq_trls_chans = nanmean(absolute_mn_acrs_trials_outseq,3);

% remove edge effects
edge = 100;
indiv_freq_inseq_trls_chans = indiv_freq_inseq_trls_chans(:, edge+1:end-edge);
indiv_freq_otseq_trls_chans = indiv_freq_otseq_trls_chans(:, edge+1:end-edge);
absolute_indiv_freq_inseq_trls_chans = absolute_indiv_freq_inseq_trls_chans(:, edge+1:end-edge);
absolute_indiv_freq_otseq_trls_chans = absolute_indiv_freq_otseq_trls_chans(:, edge+1:end-edge);

