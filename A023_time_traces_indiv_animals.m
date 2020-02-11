clear all; close all
anim = 1
lock = 'p_o' %'p_o' % p_i: poke-in % p_o: poke-out
match = '';
band = 'lo_gamma'
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

% define eventWindow/lock
if strcmp('p_i', lock)
    eventWindow = [-0.5 1.5];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn');
    
elseif strcmp('p_o', lock)
    eventWindow = [-.5 0.5];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut');
end
x = [1.1 1.2 1 1.1 1.1 .8];
response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);

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

% redefine response threshold
inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A&response_time==1;
otSeqCorrLog = otSeqLog==1 & corrTrlLog==1 & response_time==0;

% exclude trials based on time threshold
inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A& response_time==1;
otSeqCorrLog = otSeqLog==1 & corrTrlLog==1 & response_time==0;

% get incorrectly identified inseq and outseq trials exclud A
inSeqWrongLog = inSeqLog==1 & corrTrlLog==0 & response_time==0;
otSeqWrongLog = otSeqLog==1 & corrTrlLog==0 & response_time==1;

% 2 conditions of interest set as trial_1 and trial_2
trial_log_1 = inSeqCorrLog;
trial_log_2 = otSeqCorrLog;
trial_log_3 = inSeqWrongLog;
trial_log_4 = otSeqWrongLog;


% wavelet params
fs = 1000;
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

% create notch filter
flt = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);

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
[sum(trial_log_1) sum(trial_log_2) sum(trial_log_3) sum(trial_log_4)]
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

%% num of clean trials
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


%% freq ranges
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

% electrodes along CA1 axis for each animal for group average
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;

%%
cond_chans = zeros(4,length(100:size(norm_freq_acrs_chan_cond_4,2)-100), 4);
counter = 0;
for chan = 1:length(direc_chan)
    counter = counter+1;
    cond_chans(1,:,counter) = nanmean(nanmean(norm_freq_acrs_chan_cond_1(freq_range, 100:end-100,:,direc_chan(chan)),3),1);
    cond_chans(2,:,counter) = nanmean(nanmean(norm_freq_acrs_chan_cond_2(freq_range, 100:end-100,:,direc_chan(chan)),3),1);
    cond_chans(3,:,counter) = nanmean(nanmean(norm_freq_acrs_chan_cond_3(freq_range, 100:end-100,:,direc_chan(chan)),3),1);
    cond_chans(4,:,counter) = nanmean(nanmean(norm_freq_acrs_chan_cond_4(freq_range, 100:end-100,:,direc_chan(chan)),3),1);
end
%% indiv animal bar plots
start= 300
end_time = 500
mn = -1
mx = 2
elec = 1
figure
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_chans(1,start:end_time,elec)) ...
                 nanmean(cond_chans(4,start:end_time,elec)) ...
                 nanmean(cond_chans(2,start:end_time,elec)) ...
                 nanmean(cond_chans(3,start:end_time,elec)) ];
bar(bar_vectora1)
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])

elec = 2
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_chans(1,start:end_time,elec)) ...
                 nanmean(cond_chans(4,start:end_time,elec)) ...
                 nanmean(cond_chans(2,start:end_time,elec)) ...
                 nanmean(cond_chans(3,start:end_time,elec)) ];
bar(bar_vectora1)
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])

elec = 3
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_chans(1,start:end_time,elec)) ...
                 nanmean(cond_chans(4,start:end_time,elec)) ...
                 nanmean(cond_chans(2,start:end_time,elec)) ...
                 nanmean(cond_chans(3,start:end_time,elec)) ];
bar(bar_vectora1)
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])

elec = 4
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_chans(1,start:end_time,elec)) ...
                 nanmean(cond_chans(4,start:end_time,elec)) ...
                 nanmean(cond_chans(2,start:end_time,elec)) ...
                 nanmean(cond_chans(3,start:end_time,elec)) ];
bar(bar_vectora1)
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])

%%
figure
mn = -1
mx = 1
begin_idx = eventWindow(1)+.1
end_idx = eventWindow(2)-.1
win = .1*fs %100 ms smoothing windowfigure
elec = 1
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(4,:,elec),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(3,:,elec),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 2
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(4,:,elec),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(3,:,elec),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 3
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(4,:,elec),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(3,:,elec),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 4
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(4,:,elec),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(3,:,elec),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])
linkaxes
%% smoothed inseq+ outseq-
figure
mn = -1
mx = 1
begin_idx = -.4
end_idx = eventWindow(2)-.1
win = .1*fs %100 ms smoothing windowfigure
elec = 1
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 2
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 3
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 4
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(1,:,elec),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), conv(nanmean(cond_chans(2,:,elec),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

%% non-smoothed inseq+ outseq-
figure
mn = -1
mx = 1
begin_idx = -.4
end_idx = eventWindow(2)-.1
elec = 1
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)),nanmean(cond_chans(1,:,elec),1), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), nanmean(cond_chans(2,:,elec),1), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 2
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)),nanmean(cond_chans(1,:,elec),1), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), nanmean(cond_chans(2,:,elec),1), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 3
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)),nanmean(cond_chans(1,:,elec),1), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), nanmean(cond_chans(2,:,elec),1), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])

elec = 4
subplot(1,4,elec)
hold on
plot(linspace(begin_idx, end_idx,size(cond_chans,2)),nanmean(cond_chans(1,:,elec),1), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_chans,2)), nanmean(cond_chans(2,:,elec),1), 'b', 'LineWidth', 1.8)
legend({'InSeq +' 'OutSeq +' }, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(-3:.05:3)), -3:.05:3,'Color','k','LineWidth',1)
ylim ([mn mx])



