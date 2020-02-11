clear all; close all;clc
anim = 5
lock = 'p_o' %'p_o' % p_i: poke-in % p_o: poke-out
match = '';
band = 'theta'
task = 'welltrained' %novel1, novel2, welltrained
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

% load behavioral and neural data
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load('baseline_info_wavelet_32num_20logdb_3hz_250hz_notched_artifact_reject.mat')
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
fs  = 1000;
flt = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);


% get correctly identified inseq and outseq trials exluding A
times = [pokeInAlignedBehavMatrix.PokeDuration];
trial_num = 1:length(inSeqLog);
figure; plot(trial_num(inSeqLog), times(inSeqLog), 'r*') % look at response distribution
hold on
plot(trial_num(otSeqLog), times(otSeqLog), 'b*')

if strcmp('novel1', task)
    x = [1.1 .8 1 1.1];
 response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);

elseif strcmp('novel2', task)
    x = [1 1 1 1.18];
 response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);

elseif strcmp('welltrained', task)
    x = [1.1 1.2 1 1.1 1.1 .8];
 response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);
end
line([0 length(trial_num)], [x(anim), x(anim)], 'color', 'k', 'LineWidth', 2)
line([0 length(trial_num)], [1.2, 1.2], 'color', 'm', 'LineWidth', 2)

ylim([0 3])
legend('inseq', 'outseq')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
title(['anim: '  num2str(anim) ' ' task])

% redefine response threshold
inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A&response_time==1;
otSeqCorrLog = otSeqLog==1 & corrTrlLog==1 & response_time==0;

% get incorrectly identified inseq and outseq trials exclud A
inSeqWrongLog = inSeqLog==1 & corrTrlLog==0 & response_time==0;
otSeqWrongLog = otSeqLog==1 & corrTrlLog==0 & response_time==1;

% 2 conditions of interest set as trial_1 and trial_2
trial_log_1 = inSeqCorrLog;
trial_log_2 = otSeqCorrLog;
trial_log_3 = inSeqWrongLog;
trial_log_4 = otSeqWrongLog;

% control for trial num
trial_val = [sum(trial_log_1) sum(trial_log_2) sum(trial_log_3) sum(trial_log_4)];
cond_low_trial = find(trial_val==min(trial_val))
trial_min = min([sum(trial_log_1) sum(trial_log_2) sum(trial_log_3) sum(trial_log_4)]);
order = sort(trial_val);
if strcmp('yes', match)
    if anim == 3
        %match cond 1
        trial_1_idx = find(trial_log_1);
        indices = randperm(length(trial_1_idx), trial_val(4));
        trial_log_1a = zeros(1,length(trial_log_2));
        trial_log_1a(trial_1_idx(indices)) = 1;
        trial_log_1 = trial_log_1a;
        trial_log_1 = logical (trial_log_1);
        
        %match cond 2
        trial_2_idx = find(trial_log_2);
        indices = randperm(length(trial_2_idx), trial_val(4));
        trial_log_2a = zeros(1,length(trial_log_2));
        trial_log_2a(trial_2_idx(indices)) = 1;
        trial_log_2 = trial_log_2a;
        trial_log_2 = logical(trial_log_2);
        
    elseif anim ==1
                %match cond 1
        trial_1_idx = find(trial_log_1);
        indices = randperm(length(trial_1_idx), trial_val(4));
        trial_log_1a = zeros(1,length(trial_log_2));
        trial_log_1a(trial_1_idx(indices)) = 1;
        trial_log_1 = trial_log_1a;
        trial_log_1 = logical (trial_log_1);
        
        %match cond 2
        trial_2_idx = find(trial_log_2);
        indices = randperm(length(trial_2_idx), trial_val(4));
        trial_log_2a = zeros(1,length(trial_log_2));
        trial_log_2a(trial_2_idx(indices)) = 1;
        trial_log_2 = trial_log_2a;
        trial_log_2 = logical(trial_log_2);
        
        
    elseif anim ==2
                       %match cond 1
        trial_1_idx = find(trial_log_1);
        indices = randperm(length(trial_1_idx), trial_val(4));
        trial_log_1a = zeros(1,length(trial_log_2));
        trial_log_1a(trial_1_idx(indices)) = 1;
        trial_log_1 = trial_log_1a;
        trial_log_1 = logical (trial_log_1);
        
        %match cond 2
        trial_2_idx = find(trial_log_2);
        indices = randperm(length(trial_2_idx), trial_val(4));
        trial_log_2a = zeros(1,length(trial_log_2));
        trial_log_2a(trial_2_idx(indices)) = 1;
        trial_log_2 = trial_log_2a;
        trial_log_2 = logical(trial_log_2);
    else
        if trial_min>4
            if cond_low_trial==4
                
                %match cond 1
                trial_1_idx = find(trial_log_1);
                indices = randperm(length(trial_1_idx), trial_min);
                trial_log_1a = zeros(1,length(trial_log_4));
                trial_log_1a(trial_1_idx(indices)) = 1;
                trial_log_1 = trial_log_1a;
                trial_log_1 = logical (trial_log_1);
                
                %match cond 2
                trial_2_idx = find(trial_log_2);
                indices = randperm(length(trial_2_idx), trial_min);
                trial_log_2a = zeros(1,length(trial_log_4));
                trial_log_2a(trial_2_idx(indices)) = 1;
                trial_log_2 = trial_log_2a;
                trial_log_2 = logical(trial_log_2);
                
                %match cond 3
                trial_3_idx = find(trial_log_3);
                indices = randperm(length(trial_3_idx), trial_min);
                trial_log_3a = zeros(1,length(trial_log_4));
                trial_log_3a(trial_3_idx(indices)) = 1;
                trial_log_3 = trial_log_3a;
                trial_log_3 = logical(trial_log_3);
            elseif  cond_low_trial==3
                %match cond 1
                trial_1_idx = find(trial_log_1);
                indices = randperm(length(trial_1_idx), trial_min);
                trial_log_1a = zeros(1,length(trial_log_3));
                trial_log_1a(trial_1_idx(indices)) = 1;
                trial_log_1 = trial_log_1a;
                trial_log_1 = logical (trial_log_1);
                
                %match cond 2
                trial_2_idx = find(trial_log_2);
                indices = randperm(length(trial_2_idx), trial_min);
                trial_log_2a = zeros(1,length(trial_log_3));
                trial_log_2a(trial_2_idx(indices)) = 1;
                trial_log_2 = trial_log_2a;
                trial_log_2 = logical(trial_log_2);
                
                %match cond 4
                trial_4_idx = find(trial_log_4);
                indices = randperm(length(trial_4_idx), trial_min);
                trial_log_4a = zeros(1,length(trial_log_3));
                trial_log_4a(trial_4_idx(indices)) = 1;
                trial_log_4 = trial_log_4a;
                trial_log_4 = logical(trial_log_4);
            else
                trls= 're-code--didnt match'
            end
        else
            trls = 'non-matched'
        end
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
norm_freq_acrs_chan_cond_3 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_3), chan_counter);
norm_freq_acrs_chan_cond_4 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_4), chan_counter);
trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);
trial_log_3_artifact = cell(chan_counter,1);
trial_log_4_artifact = cell(chan_counter,1);

OutSeq_mean_times = mean(times(otSeqCorrLog))
OutSeq_std_times = std(times(otSeqCorrLog))
InSeq_mean_times = mean(times(inSeqCorrLog))
InSeq_std_times = std(times(inSeqCorrLog))
[OutSeq_mean_times OutSeq_std_times InSeq_mean_times InSeq_std_times]

tic
for chan = 1:chan_counter % loop thru chans
    
    % LFP data
    data = load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
    
    % extract time-locked LFP data
    LFP = ExtractTrialData_SM(pokeInAlignedBehavMatrix, filtfilt(flt,data.statMatrix(:,2)) );
    
    % extract trials' time-locked LFP data
    LFP_1 = cell2mat(LFP(trial_log_1))';
    LFP_2 = cell2mat(LFP(trial_log_2))';
    LFP_3 = cell2mat(LFP(trial_log_3))';
    LFP_4 = cell2mat(LFP(trial_log_4))';
    trial_log_1_artifact_temp= zeros(sum(trial_log_1),1);
    trial_log_2_artifact_temp= zeros(sum(trial_log_2),1);
    trial_log_3_artifact_temp= zeros(sum(trial_log_3),1);
    trial_log_4_artifact_temp= zeros(sum(trial_log_4),1);
    
    
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
                    norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/ chan_powr_std(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(cwt_power));
        end
        norm_freq_acrs_chan_cond_1(:, :, trial, chan) =  norm_freq;
        
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
                    norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/ chan_powr_std(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(cwt_power));
        end
        norm_freq_acrs_chan_cond_2(:, :, trial, chan) =  norm_freq;
        
    end
    
    trial_log_2_artifact{chan}=  trial_log_2_artifact_temp';
    
    parfor trial = 1:size(LFP_3, 1)
        % if trail has artifact = nan
        idx_above_thresh = (LFP_3(trial, :)>chan_artifact_thresh(chan,2));
        idx_below_thresh = (LFP_3(trial, :)<chan_artifact_thresh(chan,1));
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_3_artifact_temp(trial)= sum(total_idx);
        
        % wavelet spectrogram
        cwt = cwtft({LFP_3(trial, :),dt},...
            'scales',scales,'wavelet','morl');
        cwt_power = 20*log10(abs(cwt.cfs).^2);
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(cwt_power));
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/ chan_powr_std(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(cwt_power));
        end
        norm_freq_acrs_chan_cond_3(:, :, trial, chan) =  norm_freq;
        
    end
    
    trial_log_3_artifact{chan}=  trial_log_3_artifact_temp';
    
    
    parfor trial = 1:size(LFP_4, 1)
        % if trail has artifact = nan
        idx_above_thresh = (LFP_4(trial, :)>chan_artifact_thresh(chan,2));
        idx_below_thresh = (LFP_4(trial, :)<chan_artifact_thresh(chan,1));
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_4_artifact_temp(trial)= sum(total_idx);
        
        % wavelet spectrogram
        cwt = cwtft({LFP_4(trial, :),dt},...
            'scales',scales,'wavelet','morl');
        cwt_power = 20*log10(abs(cwt.cfs).^2);
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(cwt_power));
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/ chan_powr_std(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(cwt_power));
        end
        norm_freq_acrs_chan_cond_4(:, :, trial, chan) =  norm_freq;
        
    end
    
    trial_log_4_artifact{chan}=  trial_log_4_artifact_temp';
    
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

for a = 1:size(trial_log_3_artifact,1)
    chan_trials_log_3 (a)= sum(trial_log_3_artifact{a}==0);
    
end
total_tiral_log_3 = sum(chan_trials_log_3);

for a = 1:size(trial_log_4_artifact,1)
    chan_trials_log_4 (a)= sum(trial_log_4_artifact{a}==0);
    
end
total_tiral_log_4 = sum(chan_trials_log_4);



% initialize matrix for all frequencies
mn_acrs_trials_cond1 = zeros(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2), chan_counter);
mn_acrs_trials_cond2 = zeros(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2), chan_counter);
mn_acrs_trials_cond3= zeros(size(norm_freq_acrs_chan_cond_3,1), size(norm_freq_acrs_chan_cond_3,2), chan_counter);
mn_acrs_trials_cond4 = zeros(size(norm_freq_acrs_chan_cond_4,1), size(norm_freq_acrs_chan_cond_4,2), chan_counter);

% average across trials within chans
for elec = 1:size(norm_freq_acrs_chan_cond_1,4) %loop thru chans
    mn_acrs_trials_cond1 (:, :, elec)  = nanmean(norm_freq_acrs_chan_cond_1(:, :, :, elec), 3);
    mn_acrs_trials_cond2 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_2(:, :, :, elec), 3);
    mn_acrs_trials_cond3 (:, :, elec)  = nanmean(norm_freq_acrs_chan_cond_3(:, :, :, elec), 3);
    mn_acrs_trials_cond4 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_4(:, :, :, elec), 3);
end

% mean acrs chans
indiv_freq_trls_chans_cond1 = nanmean(mn_acrs_trials_cond1,3);
indiv_freq_trls_chans_cond2 = nanmean(mn_acrs_trials_cond2,3);
indiv_freq_trls_chans_cond3 = nanmean(mn_acrs_trials_cond3,3);
indiv_freq_trls_chans_cond4 = nanmean(mn_acrs_trials_cond4,3);


% spatial map indications for each anim
[medial_elecs, lateral_elecs, i, i2,...
 subplot_row_lateral, subplot_row_medial , ...
 subplot_column_lateral, subplot_column_medial] = electrode_map( anim );


%% electrode map stim locked plots plots

% freq ranges
if strcmp('delta hz', band )
    freq_range= freq>3 & freq<5;
    band_value = '3-5 hz';
elseif strcmp('theta', band )  
%         freq_range= freq>6 & freq<9;
%         band_value = '6-9 hz';
%              freq_range= freq>9 & freq<12;
%              band_value = '9-12 hz'

        freq_range= freq>11 & freq<18;
            band_value = '11-18 hz'
%     
%         freq_range= freq>6 & freq<12;
%             band_value = '6-12 hz';
elseif strcmp('lo beta', band )
    freq_range= freq>12 & freq<25;
            band_value = '12-25 hz';

elseif strcmp('beta', band )
    freq_range= freq>19 & freq<36;  % experimental cond
            band_value = '19-36 hz';
elseif strcmp('beta_gamma', band )
    freq_range= freq>19 & freq<55;
            band_value = '19-55 hz';
elseif strcmp('lo_gamma', band )
    freq_range= freq>36 & freq<55;
            band_value = '36-55 hz';

elseif strcmp('hi_gamma', band )    
    freq_range= freq>60 & freq<100;
            band_value = '60-100 hz';

elseif strcmp('ripple', band )    
    freq_range= freq>126 & freq<250;
            band_value = '126-250 hz';
end

% start and end time idx
start_time = 300;
end_time   = 500;
    
if strcmp('p_i', lock)
    end_time = round((OutSeq_mean_times-OutSeq_std_times)*fs); 
    start_time = end_time-200;
end

% electrodes along CA1 axis for each animal for group average
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;

% value for each trial/animal
% cond 1
cond1_trials = cell(4,6)
cond2_trials = cell(4,6)
cond3_trials = cell(4,6)
cond4_trials = cell(4,6)
if strcmp('novel1',task)||strcmp('novel2',task)
     load(['group_animal_individual_trials_for_stats_'  band '_nonmatched_' task])
else
    if strcmp('yes',match)
    else
        load(['group_animal_individual_trials_for_stats_'  band '_nonmatched_' task])
    end
end

for chan = 1:length(direc_chan)
    counter = 0
    for trial = 1:size(norm_freq_acrs_chan_cond_1,3)
        temp  =  nanmean(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,trial,direc_chan(chan))))
        if ~isnan(temp)
            counter = counter+1;
            cond1_trials{chan,anim}(counter) = temp;
        end
    end
end

for chan = 1:length(direc_chan)
    counter = 0
    for trial = 1:size(norm_freq_acrs_chan_cond_2,3)
        temp = nanmean(nanmean(norm_freq_acrs_chan_cond_2(freq_range, start_time:end_time,trial,direc_chan(chan))));
        if ~isnan(temp)
            counter = counter+1;
            cond2_trials{chan,anim}(counter) = temp;
        end
    end
end

for chan = 1:length(direc_chan)
        counter = 0
    for trial = 1:size(norm_freq_acrs_chan_cond_3,3)
        temp = nanmean(nanmean(norm_freq_acrs_chan_cond_3(freq_range, start_time:end_time,trial,direc_chan(chan))));
        if ~isnan(temp)
            counter = counter+1;
            cond3_trials{chan,anim}(counter) = temp;
        end
    end
end

for chan = 1:length(direc_chan)
        counter = 0

    for trial = 1:size(norm_freq_acrs_chan_cond_4,3)
        temp = nanmean(nanmean(norm_freq_acrs_chan_cond_4(freq_range, start_time:end_time,trial,direc_chan(chan))));
        if ~isnan(temp)
            counter = counter+1;
            cond4_trials{chan,anim}(counter) = temp;
        end
    end
end

if strcmp('novel1',task)||strcmp('novel2',task)
    save(['group_animal_individual_trials_for_stats_'  band '_nonmatched_' task], 'cond1_trials', 'cond2_trials', 'cond3_trials', 'cond4_trials')
else
    if strcmp('yes',match)
        save(['group_animal_individual_trials_for_stats_' band '_matched_' task], 'cond1_trials', 'cond2_trials', 'cond3_trials', 'cond4_trials')
    else
        save(['group_animal_individual_trials_for_stats_'  band '_nonmatched_' task], 'cond1_trials', 'cond2_trials', 'cond3_trials', 'cond4_trials')
    end
end
% average across trials
%bar_plot_4cond_values_per_anim = cell(1,6);

if strcmp('novel1',task)||strcmp('novel2',task)
    load(['bar_plot_4cond_values_per_anim_' band '_nonmatched_' task])
elseif strcmp('welltrained',task)
    if strcmp('yes',match)
        load(['bar_plot_4cond_values_per_anim_' band '_' band_value '_nonmatched_' task])
    else
        load(['bar_plot_4cond_values_per_anim_' band '_' band_value '_nonmatched_' task])
    end
end

    counter =0;
for chan = 1:length(direc_chan)
    counter = counter+1
    bar_plot_4cond_values_per_anim{anim}(counter, 1) =  nanmean(nanmean(mn_acrs_trials_cond1(freq_range, start_time:end_time,direc_chan(chan)),2))
    bar_plot_4cond_values_per_anim{anim}(counter, 2) =  nanmean(nanmean(mn_acrs_trials_cond2(freq_range, start_time:end_time,direc_chan(chan)),2))
    bar_plot_4cond_values_per_anim{anim}(counter, 3) =  nanmean(nanmean(mn_acrs_trials_cond3(freq_range, start_time:end_time,direc_chan(chan)),2))
    bar_plot_4cond_values_per_anim{anim}(counter, 4) =  nanmean(nanmean(mn_acrs_trials_cond4(freq_range, start_time:end_time,direc_chan(chan)),2))
end

if strcmp('yes',match)
    save(['bar_plot_4cond_values_per_anim_' band '_' band_value '_matched_' task], 'bar_plot_4cond_values_per_anim')
else
    save(['bar_plot_4cond_values_per_anim_' band '_' band_value '_nonmatched_' task], 'bar_plot_4cond_values_per_anim')
end
%% preview indiv animal bar plots
figure
for elec = 1:4
subplot(1, 4, elec)
bar_vectora1 = [bar_plot_4cond_values_per_anim{anim}(elec, 1) bar_plot_4cond_values_per_anim{anim}(elec, 4) ...
                 bar_plot_4cond_values_per_anim{anim}(elec, 2) bar_plot_4cond_values_per_anim{anim}(elec, 3)];
bar(bar_vectora1)
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([-1 1])
end

%% calculate pairwise ratios
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
load('all_proportion_all_anims')
out_seq_corr_proportion = bar_plot_4cond_values_per_anim{anim}(:,2)./bar_plot_4cond_values_per_anim{anim}(:,1);
out_seq_incorr_proportion = bar_plot_4cond_values_per_anim{anim}(:,4)./bar_plot_4cond_values_per_anim{anim}(:,1);
all_out_seq_corr_proportion(:,anim) = out_seq_corr_proportion
all_out_seq_incorr_proportion(:,anim) = out_seq_incorr_proportion
save('all_proportion_all_anims', 'all_out_seq_corr_proportion', 'all_out_seq_incorr_proportion')


%% after running paiwise ratios for all animals, measure and plot percentage of conditions w/respect to InSeq
animal_length = 6
load('all_proportion_all_anims')
mean_out_seq_corr = -1*(1-nanmean(all_out_seq_corr_proportion,2))
mean_out_seq_incorr = -1*(1-nanmean(all_out_seq_incorr_proportion,2))
std_out_seq_corr = (nanstd(all_out_seq_corr_proportion,0,2))
std_out_seq_incorr = (nanstd(all_out_seq_incorr_proportion,0,2))
mn = -2
mx =  1

figure
for elec = 1:4
    subplot(1, 4, elec)
    bar_vectora1 = [mean_out_seq_corr(elec) mean_out_seq_incorr(elec)];
    bar(bar_vectora1)
    hold on
    errorbar(1:2,bar_vectora1,[std_out_seq_corr(elec)/sqrt(animal_length) std_out_seq_incorr(elec)/sqrt(animal_length)], 'rx')
    ylim([mn mx])
     
  %  set(gca, 'XTick', 1:2, 'XTickLabel', {'OutSeq +' 'OutSeq -' }, 'YTickLabel', [200 150 100 50 0 50], 'XTickLabelRotation',45)
    ylabel('% change from InSeq +')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
end

%% pairwise differences
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
clear out_seq_corr_diff out_seq_incorr_diff all_out_seq_corr_diff all_out_seq_incorr_diff
load('all_differences_all_anims')
out_seq_corr_diff = bar_plot_4cond_values_per_anim{anim}(:,2)-bar_plot_4cond_values_per_anim{anim}(:,1);
out_seq_incorr_diff = bar_plot_4cond_values_per_anim{anim}(:,4)-bar_plot_4cond_values_per_anim{anim}(:,1);
all_out_seq_corr_diff(:,anim) = out_seq_corr_diff
all_out_seq_incorr_diff(:,anim) = out_seq_incorr_diff
save('all_differences_all_anims', 'all_out_seq_corr_diff', 'all_out_seq_incorr_diff')

%% after running diff for all animals
clear out_seq_corr_diff out_seq_incorr_diff all_out_seq_corr_diff all_out_seq_incorr_diff
load('all_differences_all_anims')
mean_out_seq_corr_all_anim = nanmean(all_out_seq_corr_diff,2)
mean_out_seq_incorr_all_anim = nanmean(all_out_seq_incorr_diff,2)
std_out_seq_corr_all_anim = nanmean(all_out_seq_corr_diff,2)
std_out_seq_incorr_all_anim = nanmean(all_out_seq_incorr_diff,2)
mn = -.3
mx = .2
figure
for elec = 1:4
    subplot(1, 4, elec)
    bar_vectora1 = [mean_out_seq_corr_all_anim(elec) mean_out_seq_incorr_all_anim(elec)];
    bar(bar_vectora1)
    hold on
    errorbar(1:2,bar_vectora1,[std_out_seq_corr_all_anim(elec)/sqrt(animal_length) std_out_seq_incorr_all_anim(elec)/sqrt(animal_length)], 'rx')
    ylim([mn mx])   
    set(gca, 'XTick', 1:2, 'XTickLabel', {'OutSeq +' 'OutSeq -' }, 'XTickLabelRotation',45)
    ylabel('Difference from InSeq +')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
end

%% extract traces from each animal
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
clear inseq_corr_trace otseq_corr_trace inseq_incorr_trace otseq_incorr_trace
load('group_traces')
animal_length = 6
% inseq_corr_trace = nan(animal_length, size(indiv_freq_trls_chans_cond1,2), length(direc_chan));
% otseq_corr_trace = nan(animal_length, size(indiv_freq_trls_chans_cond1,2), length(direc_chan));
% inseq_incorr_trace = nan(animal_length, size(indiv_freq_trls_chans_cond1,2), length(direc_chan));
% otseq_incorr_trace = nan(animal_length, size(indiv_freq_trls_chans_cond1,2), length(direc_chan));

counter =0;
for chan = 1:length(direc_chan)
    counter = counter+1
    inseq_corr_trace   (anim,:, chan) = nanmean(mn_acrs_trials_cond1(freq_range, :,direc_chan(chan)),1);
    otseq_corr_trace   (anim,:, chan) = nanmean(mn_acrs_trials_cond2(freq_range, :,direc_chan(chan)),1);
    inseq_incorr_trace   (anim,:, chan) = nanmean(mn_acrs_trials_cond3(freq_range, :,direc_chan(chan)),1);
    otseq_incorr_trace   (anim,:, chan) = nanmean(mn_acrs_trials_cond4(freq_range, :,direc_chan(chan)),1);
    
end

save('group_traces', 'inseq_corr_trace', 'otseq_corr_trace', 'inseq_incorr_trace',  'otseq_incorr_trace')



%% after getting traces for all animals
clear inseq_corr_trace_all_anim otseq_corr_trace_all_anim inseq_incorr_trace_all_anim otseq_incorr_trace_all_anim
load('group_traces')
inseq_corr_trace_all_anim = squeeze(nanmean(inseq_corr_trace,1))';
otseq_corr_trace_all_anim = squeeze(nanmean(otseq_corr_trace,1))';
inseq_incorr_trace_all_anim = squeeze(nanmean(inseq_incorr_trace,1))';
otseq_incorr_trace_all_anim = squeeze(nanmean(otseq_incorr_trace,1))';

%% plot group time traces
win = 200
mn  = -.5
mx  = .8
figure
for elec = 1:4
subplot(1,4,elec)
plot(linspace(eventWindow(1), eventWindow(2),size(inseq_corr_trace,2)),conv(inseq_corr_trace_all_anim(elec,:), ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
hold on
plot(linspace(eventWindow(1), eventWindow(2),size(inseq_corr_trace,2)),conv(otseq_incorr_trace_all_anim(elec,:), ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(inseq_corr_trace,2)),conv(otseq_corr_trace_all_anim(elec,:), ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(inseq_corr_trace,2)),conv(inseq_incorr_trace_all_anim(elec,:), ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'InSeq +', 'OutSeq -' , 'OutSeq +' , 'InSeq -'  }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',.8)
line([-.25 -.25], [mn mx], 'color', 'k', 'LineWidth', .25)
ylim ([mn mx])
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
end

%% inseq incorr
load('all_in_seq_incorr_proportion')
in_seq_incorr_proportion = bar_plot_4cond_values_per_anim{anim}(:,3)./bar_plot_4cond_values_per_anim{anim}(:,1);
all_in_seq_incorr_proportion(:,anim) = in_seq_incorr_proportion
save('all_in_seq_incorr_proportion','all_in_seq_incorr_proportion')

InSeq_Incorr_proportion = nanmean(all_in_seq_incorr_proportion,2)
err = nanstd(all_in_seq_incorr_proportion,0,2)
%%
mn = -2.2
mx = .5
figure
for elec = 1:4
    subplot(1, 4, elec)
    bar(InSeq_Incorr_proportion(elec));
    hold on
    errorbar(InSeq_Incorr_proportion(elec),err(elec), 'rx')
    ylim([mn mx])
    set(gca, 'XTick', 1, 'XTickLabel', {'InSeq -' },'YTickLabel', [200 150 100 50 0 50], 'XTickLabelRotation',45) %, 
    ylabel('% change from InSeq +')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
end