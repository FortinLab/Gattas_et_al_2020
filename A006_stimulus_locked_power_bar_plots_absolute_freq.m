clear all; close all
anim = 1
lock = 'p_o' %'p_o' % p_i: poke-in % p_o: poke-out
match = '';
band = 'beta'
task = 'welltrained' %novel1, novel2, welltrained
per = 0.5; % caxis variable for plotting
scale = 'db'
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
ylim([0 3])
legend('inseq', 'outseq')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
title(['anim: '  num2str(anim) ' ' task])

%figure; plot( times, 'r*') % look at response distribution

% redefine response threshold

inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A&response_time==1;
otSeqCorrLog = otSeqLog&corrTrlLog&response_time==0;

% corr vs. incorr exld A
corrTrlLog_exld_A   = response_time&odor_excld_A==1;
incorrTrlLog_exld_A = response_time==0&odor_excld_A==1;

% get incorrectly identified inseq and outseq trials exclud A
inSeqWrongLog = inSeqLog==1 & corrTrlLog==0 & odor_excld_A==1 & response_time==0;
otSeqWrongLog = otSeqLog==1 & corrTrlLog==0 & odor_excld_A==1 & response_time==1;

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
abs_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_1), chan_counter);
abs_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_2), chan_counter);
abs_freq_acrs_chan_cond_3 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_3), chan_counter);
abs_freq_acrs_chan_cond_4 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_4), chan_counter);
trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);
trial_log_3_artifact = cell(chan_counter,1);
trial_log_4_artifact = cell(chan_counter,1);
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
        cwt_power = 20*log10(abs(cwt.cfs).^2);
        if sum(total_idx)<1 % if trial is artifact free
            norm_freq = cwt_power; 
        else
            norm_freq = nan(size(cwt_power));
        end
        abs_freq_acrs_chan_cond_2(:, :, trial, chan) =  norm_freq;
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
            norm_freq =cwt_power;
        else
            norm_freq = nan(size(cwt_power));
        end
        abs_freq_acrs_chan_cond_3(:, :, trial, chan) =  norm_freq;
        
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
            norm_freq =cwt_power;
        else
            norm_freq = nan(size(cwt_power));
        end
        abs_freq_acrs_chan_cond_4(:, :, trial, chan) =  norm_freq;
        
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
mn_acrs_trials_cond1 = zeros(size(abs_freq_acrs_chan_cond_1,1), size(abs_freq_acrs_chan_cond_1,2), chan_counter);
mn_acrs_trials_cond2 = zeros(size(abs_freq_acrs_chan_cond_2,1), size(abs_freq_acrs_chan_cond_2,2), chan_counter);
mn_acrs_trials_cond3= zeros(size(abs_freq_acrs_chan_cond_3,1), size(abs_freq_acrs_chan_cond_3,2), chan_counter);
mn_acrs_trials_cond4 = zeros(size(abs_freq_acrs_chan_cond_4,1), size(abs_freq_acrs_chan_cond_4,2), chan_counter);

% average across trials within chans
for elec = 1:size(abs_freq_acrs_chan_cond_1,4) %loop thru chans
    mn_acrs_trials_cond1 (:, :, elec)  = nanmean(abs_freq_acrs_chan_cond_1(:, :, :, elec), 3);
    mn_acrs_trials_cond2 (:, :, elec) = nanmean(abs_freq_acrs_chan_cond_2(:, :, :, elec), 3);
    mn_acrs_trials_cond3 (:, :, elec)  = nanmean(abs_freq_acrs_chan_cond_3(:, :, :, elec), 3);
    mn_acrs_trials_cond4 (:, :, elec) = nanmean(abs_freq_acrs_chan_cond_4(:, :, :, elec), 3);
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


% electrode map stim locked plots plots

% freq ranges
if strcmp('delta hz', band )
    freq_range= freq>3 & freq<5;
    band_value = '3-5 hz';
elseif strcmp('theta', band )  
%         freq_range= freq>6 & freq<9;
%         band_value = '6-9 hz';
%     
%         freq_range= freq>11 & freq<19;
%             band_value = '11-19 hz';
%     
%         freq_range= freq>6 & freq<12;
%             band_value = '6-12 hz';

             freq_range= freq>9 & freq<12;
             band_value = '9-12 hz';
elseif strcmp('lo beta', band )
    freq_range= freq>12 & freq<25;
            band_value = '12-25 hz';

elseif strcmp('beta', band )
    freq_range= freq>19 & freq<36;  % experimental cond
            band_value = '19-36 hz';

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
if sum(trial_log_2)==0 || min(times(trial_log_2))>.5
    start_time = 1;
else
    start_time = (.5- min(times(trial_log_2)))*fs;
end
 end_time = 500;
 
 % colorbar axis
minim = [-.5 -.5 -.2 -.5 -.5 -.5];
max = [1 .5 .9 1 .5 .5];
 %
%  figure
%     i_counter = 0;
%     for chan = lateral_elecs
%             i_counter = i_counter+1;
%             
%         low_gamma_cond1(i_counter) = mean(mean(mn_acrs_trials_cond1(freq_range, start_time:end_time,chan==chan_length),2));
%         low_gamma_cond2(i_counter) = mean(mean(mn_acrs_trials_cond2(freq_range,start_time:end_time,chan==chan_length),2));
%         low_gamma_cond3(i_counter) = mean(mean(mn_acrs_trials_cond3(freq_range, start_time:end_time,chan==chan_length),2));
%         low_gamma_cond4(i_counter) = mean(mean(mn_acrs_trials_cond4(freq_range,start_time:end_time,chan==chan_length),2));
%         subplot(subplot_row_lateral(anim), subplot_column_lateral(anim), i2(i_counter))
%         bar([low_gamma_cond1(i_counter)  low_gamma_cond4(i_counter) low_gamma_cond2(i_counter) low_gamma_cond3(i_counter)])
%         set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
%         ylim([minim(anim) max(anim)])
%         title (['chan' num2str(chan)])
%         set(gca, 'FontSize',14,'FontWeight','Bold')
%     end
% 
%         figure
%         i_counter = 0;
%     for chan = medial_elecs
%             i_counter = i_counter+1;
%         low_gamma_cond1(i_counter) = mean(mean(mn_acrs_trials_cond1(freq_range, start_time:end_time,chan==chan_length),2));
%         low_gamma_cond2(i_counter) = mean(mean(mn_acrs_trials_cond2(freq_range,start_time:end_time,chan==chan_length),2));
%         low_gamma_cond3(i_counter) = mean(mean(mn_acrs_trials_cond3(freq_range, start_time:end_time,chan==chan_length),2));
%         low_gamma_cond4(i_counter) = mean(mean(mn_acrs_trials_cond4(freq_range,start_time:end_time,chan==chan_length),2));
%         subplot(subplot_row_medial(anim), subplot_column_medial(anim), i(i_counter))
%         bar([low_gamma_cond1(i_counter)  low_gamma_cond4(i_counter) low_gamma_cond2(i_counter) low_gamma_cond3(i_counter)])
%         set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
%         ylim([minim(anim) max(anim)])
%         title (['chan' num2str(chan)])
%         set(gca, 'FontSize',14,'FontWeight','Bold')
%     end
    
% electrodes along CA1 axis for each animal for group average
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;

% value for each trial/animal
% cond 1
% cond1_trials = cell(4,6)
% cond2_trials = cell(4,6)
% cond3_trials = cell(4,6)
% cond4_trials = cell(4,6)

 %% 
% if strcmp('novel1',task)||strcmp('novel2',task)
%      load(['group_animal_individual_trials_for_stats_'  band '_nonmatched_' task])
% else
%     if strcmp('yes',match)
%     else
%         load(['group_animal_individual_trials_for_stats_'  band '_nonmatched_' task])
%     end
% end
for chan = 1:length(direc_chan)
    counter = 0
    for trial = 1:size(abs_freq_acrs_chan_cond_1,3)
        temp  =  nanmean(nanmean(abs_freq_acrs_chan_cond_1(freq_range, start_time:end_time,trial,direc_chan(chan))))
        if ~isnan(temp)
            counter = counter+1;
            cond1_trials{chan,anim}(counter) = temp;
        end
    end
end

for chan = 1:length(direc_chan)
    counter = 0
    for trial = 1:size(abs_freq_acrs_chan_cond_2,3)
        temp = nanmean(nanmean(abs_freq_acrs_chan_cond_2(freq_range, start_time:end_time,trial,direc_chan(chan))));
        if ~isnan(temp)
            counter = counter+1;
            cond2_trials{chan,anim}(counter) = temp;
        end
    end
end

for chan = 1:length(direc_chan)
        counter = 0
    for trial = 1:size(abs_freq_acrs_chan_cond_3,3)
        temp = nanmean(nanmean(abs_freq_acrs_chan_cond_3(freq_range, start_time:end_time,trial,direc_chan(chan))));
        if ~isnan(temp)
            counter = counter+1;
            cond3_trials{chan,anim}(counter) = temp;
        end
    end
end

for chan = 1:length(direc_chan)
        counter = 0

    for trial = 1:size(abs_freq_acrs_chan_cond_4,3)
        temp = nanmean(nanmean(abs_freq_acrs_chan_cond_4(freq_range, start_time:end_time,trial,direc_chan(chan))));
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
bar_plot_4cond_values_per_anim = cell(1,6);

if strcmp('novel1',task)||strcmp('novel2',task)
    load(['bar_plot_4cond_values_per_anim_' band '_nonmatched_' task '_absolute_freq'])
elseif strcmp('welltrained',task)
    if strcmp('yes',match)
        load(['bar_plot_4cond_values_per_anim_' band '_' band_value '_matched_' task '_absolute_freq'])
    else
        load(['bar_plot_4cond_values_per_anim_' band '_' band_value '_nonmatched_' task '_absolute_freq'])
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
    save(['bar_plot_4cond_values_per_anim_' band '_' band_value '_matched_' task '_absolute_freq'], 'bar_plot_4cond_values_per_anim')
else
    save(['bar_plot_4cond_values_per_anim_' band '_' band_value '_nonmatched_' task '_absolute_freq'], 'bar_plot_4cond_values_per_anim')
end