clear all; close all; clc
anim            = 6
trial_selection = 'run'
task            = 'welltrained' %novel1, novel2, welltrained
thresh          = 1.5;
freq_states     =''% 'SeqMem'
NumVoices       = 32;
fs              = 1000;
lock            = 'p_i'

% load behavioral and neural data
[chan_name, chan_length ]= get_anim_info(anim, task);
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])
chan_counter = length(chan_length);
load([chan_name '4.mat']) % load exmp chan
if strcmp('SeqMem', freq_states)
    load('baseline_wavelet_32num_10logdb_3hz_250hz_welltrained_avg_freq_ODOR_RUN_states_SeqMemFreq.mat')
else
    %   load(['baseline_wavelet_32num_10logdb_3hz_250hz_' task '_avg_freq_ODOR_RUN_states'])
end
load(['baseline_info_wavelet_' num2str(NumVoices) 'num_10logdb_3hz_250hz_notched_artifact_reject.mat'])
load('BehaviorMatrix.mat')

%find position on maze with largest num of trials
[ middle_maze_idx_final ] = running_index( anim ,behavMatrix);

% condition 2 data: poke out - random selection of trials
if strcmp('p_i', lock)
    eventWindow = [-0.5 .5];
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


inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A&response_time==1;
inCorrLog    =  corrTrlLog==0;

% trial_log       = inCorrLog;
% idx = find(inSeqCorrLog);
% inSeqCorrLog(idx(23)+1:end) = 0;
trial_log       = inSeqCorrLog;
maze_run_trials = randperm(size(middle_maze_idx_final,1),sum(trial_log));
middle_maze_idx_final = middle_maze_idx_final(maze_run_trials);
LFP             = ExtractTrialData_SM(pokeInAlignedBehavMatrix, statMatrix(:,2));
LFP_trials_data = cell2mat(LFP(trial_log))';
for idx = 1:length(middle_maze_idx_final)
    pos_maz_data(idx, :) = statMatrix(middle_maze_idx_final(idx)-(fs/2):middle_maze_idx_final(idx)+(fs/2), 2);
end

% create wavelets
dt = 1/fs;
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

if strcmp('SeqMem', freq_states)
    freq_idx (:,1) = freq>6 & freq<9;
    freq_idx (:,2) = freq>9 & freq<12;
    freq_idx (:,3) = freq>11 & freq<19;
    freq_idx (:,4) = freq>19 & freq<36;
    freq_idx (:,5) = freq>36 & freq<55;
    freq_idx (:,6) = freq>60 & freq<100;
    freq_idx (:,7) = freq>126 & freq<250;
    cond_1_event_histogram     = zeros(size(freq_idx,2),size(cwt_power_exemp,2), sum(trial_log), chan_counter);
    cond_2_event_histogram     =  zeros(size(freq_idx,2),size(cwt_power_exemp,2),sum(trial_log), chan_counter);
else strcmp('yes', freq_states)
    freq_idx (:,1) = freq>3 & freq<5;
    freq_idx (:,2) = freq>4 & freq<12;
    freq_idx (:,3) = freq>12 & freq<25;
    freq_idx (:,4) = freq>25 & freq<55;
    freq_idx (:,5) = freq>60 & freq<100;
    freq_idx (:,6) = freq>126 & freq<250;
    cond_1_event_histogram     = zeros(size(freq_idx,2),size(cwt_power_exemp,2), sum(trial_log), chan_counter);
    cond_2_event_histogram     =  zeros(size(freq_idx,2),size(cwt_power_exemp,2),sum(trial_log), chan_counter);
end

% run wavelet exemplar
cwt             = cwtft({pos_maz_data(1, :),dt},...
    'scales',scales,'wavelet','morl');

cwt_power_exemp = 10*log10(abs(cwt.cfs).^2);

norm_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(pos_maz_data,1), chan_counter);

fs  = 1000;
flt1 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);
%%
for chan = 1:chan_counter
    
    % LFP data
    data = load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
    data_flt_1 =  filtfilt(flt1,data.statMatrix(:,2));
    for idx = 1:length(middle_maze_idx_final)
        pos_maz_data(idx, :) = data_flt_1(middle_maze_idx_final(idx)-(fs/2):middle_maze_idx_final(idx)+(fs/2));
    end
    
    % calculate logicals for all trials in condition 1
    parfor trial = 1:size(pos_maz_data,1)
        current_trial_logical = zeros(size(freq_idx,2), size(cwt_power_exemp,2));
        
        % wavelet spectrogram
        cwt = cwtft({pos_maz_data(trial, :),dt},...
            'scales',scales,'wavelet','morl');
        cwt_power = 10*log10(abs(cwt.cfs).^2);
        
        % normalize z-score relative to 1-hour session
        norm_freq = zeros(size(cwt_power));
        for a = 1:size(cwt_power,1) %loop thru freq
            for b = 1:size(cwt_power,2)%loop thru timepoints
                norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/chan_powr_std(a, chan);
            end
        end
        norm_freq_acrs_chan_cond_1(:,:,trial,chan) =  norm_freq;
        
        if strcmp('SeqMem', freq_states)
            
            theta1       = mean(cwt_power((freq>6 & freq<9), :));
            theta2       = mean(cwt_power((freq>9 & freq<12), :));
            alpha        = mean(cwt_power((freq>11 & freq<19), :));
            beta         = mean(cwt_power((freq>19 & freq<36), :));
            slow_gamma   = mean(cwt_power((freq>36 & freq<55), :));
            fast_gamma   = mean(cwt_power((freq>60 & freq<100), :));
            ripple       = mean(cwt_power((freq>126 & freq<250), :));
            current_trial = [theta1; theta2; alpha; beta; slow_gamma; fast_gamma; ripple];
            
        elseif strcmp('yes', freq_states)
            
            %  average into feq bands of interes
            delta         = mean(cwt_power((freq>3 & freq<5), :));
            theta         = mean(cwt_power((freq>4 & freq<12), :));
            beta          = mean(cwt_power((freq>12 & freq<25), :));
            slow_gamma    = mean(cwt_power((freq>25 & freq<55), :));
            fast_gamma    = mean(cwt_power((freq>60 & freq<100), :));
            ripple        = mean(cwt_power((freq>126 & freq<250), :));
            current_trial = [delta; theta; beta; slow_gamma; fast_gamma; ripple];
            
            
            % find indices above threshold
            for freqs = 1:size(current_trial,1) %loop thru freq
                for timepoints = 1:size(current_trial,2) %loop thru trial  timepoitns
                    if current_trial(freqs, timepoints)>freq_mean_chan(chan,freqs)+(thresh*freq_std_chan(chan,freqs));
                        current_trial_logical(freqs,timepoints) = 1;
                    end
                end
            end
            cond_1_event_histogram (:, :, trial, chan) = current_trial_logical;
            
        end
    end
    
    %
    %% condition 2: odor presentation
    %extract time-locked LFP data
    LFP = ExtractTrialData_SM(pokeInAlignedBehavMatrix, filtfilt(flt1,data.statMatrix(:,2)) );
    LFP_trials_data = cell2mat(LFP(trial_log))';
    
    %calculate logicals for all trials in condition 1
    for trial = 1:size(LFP_trials_data,1)
        current_trial_logical = zeros(size(freq_idx,2), size(cwt_power_exemp,2));
        
        % wavelet spectrogram
        cwt = cwtft({LFP_trials_data(trial, :),dt},...
            'scales',scales,'wavelet','morl');
        cwt_power = 10*log10(abs(cwt.cfs).^2);
        
        % average into feq bands of interes
        if strcmp('SeqMem', freq_states)
            
            theta1       = mean(cwt_power((freq>6 & freq<9), :));
            theta2       = mean(cwt_power((freq>9 & freq<12), :));
            alpha        = mean(cwt_power((freq>11 & freq<19), :));
            beta        = mean(cwt_power((freq>19 & freq<36), :));
            slow_gamma  = mean(cwt_power((freq>36 & freq<55), :));
            fast_gamma  = mean(cwt_power((freq>60 & freq<100), :));
            ripple      = mean(cwt_power((freq>126 & freq<250), :));
            current_trial = [theta1; theta2; alpha; beta; slow_gamma; fast_gamma; ripple];
            
        elseif strcmp('yes', freq_states)
            %  average into feq bands of interes
            delta         = mean(cwt_power((freq>3 & freq<5), :));
            theta         = mean(cwt_power((freq>4 & freq<12), :));
            beta          = mean(cwt_power((freq>12 & freq<25), :));
            slow_gamma    = mean(cwt_power((freq>25 & freq<55), :));
            fast_gamma    = mean(cwt_power((freq>60 & freq<100), :));
            ripple        = mean(cwt_power((freq>126 & freq<250), :));
            current_trial = [delta; theta; beta; slow_gamma; fast_gamma; ripple];
            
            % find indices above threshold
            for freqs = 1:size(current_trial,1) %loop thru freq
                for timepoints = 1:size(current_trial,2) %loop thru trial  timepoitns
                    if current_trial(freqs, timepoints)>freq_mean_chan(chan,freqs)+(thresh*freq_std_chan(chan,freqs));
                        current_trial_logical(freqs,timepoints) = 1;
                    end
                end
            end
            cond_2_event_histogram (:, :, trial, chan) = current_trial_logical;
        end
    end
end
% re-organizing running spectrogram data
% initialize matrix for all frequencies
mn_acrs_trials_run = zeros(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2), chan_counter);

% average across trials within chans
for elec = 1:size(norm_freq_acrs_chan_cond_1,4) %loop thru chans
    mn_acrs_trials_run (:, :, elec)  = nanmean(norm_freq_acrs_chan_cond_1(:, :, :, elec), 3);
end

% mean acrs chans
indiv_freq_run_trls_chans = nanmean(mn_acrs_trials_run,3);

% save data w. all trials for each anim
cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim\fig2_odorVsrun\run_matrices_AL_PM')
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;
run_spectrogram_AL_PM = norm_freq_acrs_chan_cond_1(:, :, :, direc_chan);
save(['anim' num2str(anim) '_AL_PM_' trial_selection], 'run_spectrogram_AL_PM');
%% plot running spectrogram
% initialize matrix for all frequencies
mn_acrs_trials = zeros(size(norm_freq_acrs_chan_cond_1,1),801, chan_counter);

% average across trials within chans
for elec = 1:size(norm_freq_acrs_chan_cond_1,4) %loop thru chans
    mn_acrs_trials(:, :, elec)  = mean(norm_freq_acrs_chan_cond_1(:, 101:end-100, :, elec), 3);
end
indiv_freq_inseq_trls_chans = mean(mn_acrs_trials,3);
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;
cond1 = mn_acrs_trials(:, :, direc_chan);

% plot 
for chan = 1:chan_counter
    figure
    tickmarks = 1:30:length(freq);
    imagesc(linspace(-.5+.0017,.5-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq), mn_acrs_trials(:, :, chan))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([min(min([ mn_acrs_trials(50:end,:,chan) ])) max(max(([ mn_acrs_trials(50:end,:,chan) ])))])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    saveas(gcf,['chan' num2str(chan) 'running_spectrogram.png'])
end

%% raster plot running trials:
if strcmp('SeqMem', freq_states)
    band_name  = {'theta1','theta2', 'alpha', 'beta','lo gamma','hi gamma','ripple'};
    band_range = {'6-9 hz', '9-12 hz', '11-19 hz', '19-36 hz', '36-55 hz', '60-100 hz','125-250 hz'};
else
    band_name  = {'delta','theta','beta','lo gamma','hi gamma','ripple'};
    band_range = {'3-5 hz', '4-12 hz', '12-25 hz', '25-55 hz', '60-100 hz', '125-250 hz'};
end


trial_data_1 = zeros(size(cond_1_event_histogram,3), size(cond_2_event_histogram,2), length(band_name));
trial_data_2 = zeros(size(cond_1_event_histogram,3), size(cond_2_event_histogram,2), length(band_name));
chan_trial_1_data = cell(1,chan_counter);
chan_trial_2_data = cell(1,chan_counter);

% save data per chan
for chan = 1:size(cond_2_event_histogram,4) % loop thru chans
    for bands = 1:size(band_name,2) % loop thru bands
        for trial = 1:size(cond_2_event_histogram, 3)
            trial_data_1(trial, :, bands) = cond_1_event_histogram(bands, :,trial, chan);
            trial_data_2(trial, :, bands) = cond_2_event_histogram(bands, :,trial, chan);

        end
    end
    chan_trial_1_data{chan} = trial_data_1;
    chan_trial_2_data{chan} = trial_data_2;
end

[anterolat_postermedial] = ca1_electodes(anim, chan_length)
direc_chan = anterolat_postermedial;
ca1_elec_length = 4;

%% plot third elec low and hi gamma for running
third_elec = direc_chan(3)
bnd = find(strcmp('lo gamma',band_name))
figure
subplot(2,2,1)
temp = chan_trial_1_data{third_elec}(:,:,bnd);
for trial = 1:size(cond_2_event_histogram, 3) % loop thru trials
    log_vector = logical(temp(trial,:)) ;
    vector_value = zeros(1,length(log_vector));
    vector_value(log_vector) = trial;
    vector_value(~log_vector) = -5;
    
    plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))), vector_value(51:end-50), '.g')
    xlim([-.5+.05 .5-.05])
    ylim([0 size(cond_2_event_histogram, 3)])
    xlabel('time from center of maze')
    ylabel('trials')
    hold on
end
title([band_name{bnd} ' band: ' band_range{bnd}])
set(gca, 'FontSize',12,'FontWeight','bold')

% smooth traces
win_size = 40;
sum_trace = sum(chan_trial_1_data{third_elec}(:,51:end-50,bnd),1);
sum_trace_smootha = conv(sum_trace, ones(1,win_size)/win_size, 'same');
clear sum_trace
bnd = find(strcmp('hi gamma',band_name))
sum_trace = sum(chan_trial_1_data{third_elec}(:,51:end-50,bnd),1);
sum_trace_smoothb=conv(sum_trace, ones(1,win_size)/win_size, 'same');
mx = max([sum_trace_smoothb sum_trace_smootha]);


subplot(2,2,3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smootha,'g','LineWidth',2.5)
xlim([-.5+.05 .5-.05])
ylim([0 mx])
set(gca, 'FontSize',12,'FontWeight','bold')

bnd = find(strcmp('hi gamma',band_name))
subplot(2,2,2)
temp = chan_trial_1_data{third_elec}(:,:,bnd);
for trial = 1:size(cond_2_event_histogram, 3) % loop thru trials
    log_vector = logical(temp(trial,:)) ;
    vector_value = zeros(1,length(log_vector));
    vector_value(log_vector) = trial;
    vector_value(~log_vector) = -5;
    
    plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))), vector_value(51:end-50), '.r')
    xlim([-.5+.05 .5-.05])
    ylim([0 size(cond_2_event_histogram, 3)])
    xlabel('time from center of maze')
    ylabel('trials')
    hold on
end
title([band_name{5} ' band: ' band_range{bnd}])
set(gca, 'FontSize',12,'FontWeight','bold')

sum_trace = sum(chan_trial_1_data{third_elec}(:,51:end-50,bnd),1);
sum_trace_smoothb=conv(sum_trace, ones(1,win_size)/win_size, 'same');
mx = max([sum_trace_smoothb sum_trace_smootha]);

subplot(2,2,4)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothb,  'r', 'LineWidth',2.5)
xlim([-.5+.05 .5-.05])
ylim([0 mx])
set(gca, 'FontSize',12,'FontWeight','bold')

figure
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smootha,  'g', 'LineWidth',2.5)
hold on
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothb,  'r', 'LineWidth',2.5)
xlim([-.5+.05 .5-.05])
ylim([0 mx])
legend({'low gamma', 'high gamma'})
title('running')
set(gca, 'FontSize',12,'FontWeight','bold')


%% plot third elec low and hi gamma for poke
third_elec = direc_chan(3)
bnd = find(strcmp('lo gamma',band_name))
figure
subplot(2,2,1)
temp = chan_trial_2_data{third_elec}(:,:,bnd);
for trial = 1:size(cond_2_event_histogram, 3) % loop thru trials
    log_vector = logical(temp(trial,:)) ;
    vector_value = zeros(1,length(log_vector));
    vector_value(log_vector) = trial;
    vector_value(~log_vector) = -5;
    plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))), vector_value(51:end-50), '.g')
    xlim([-.5+.05 .5-.05])
    ylim([0 size(cond_2_event_histogram, 3)])
    hold on
    plot(zeros(1,length(freq)+1), 0:length(freq),'Color','k','LineWidth',1)
    xlabel('time ')
    ylabel('trials')
    hold on
end
title([band_name{bnd} ' band: ' band_range{bnd}])
set(gca, 'FontSize',12,'FontWeight','bold')

% moving average
sum_trace = sum(chan_trial_2_data{third_elec}(:,51:end-50,bnd),1);
sum_trace_smootha =  conv(sum_trace, ones(1,win_size)/win_size, 'same');
clear sum_trace
subplot(2,2,3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smootha,'g','LineWidth',2.5)
bnd = find(strcmp('hi gamma',band_name))
sum_trace = sum(chan_trial_2_data{third_elec}(:,51:end-50,bnd),1);
sum_trace_smoothb = conv(sum_trace, ones(1,win_size)/win_size, 'same');
mx = max([sum_trace_smoothb sum_trace_smootha]);
xlim([-.5+.05 .5-.05])
ylim([0 mx])
set(gca, 'FontSize',12,'FontWeight','bold')


bnd = find(strcmp('hi gamma',band_name))
subplot(2,2,2)
temp = chan_trial_2_data{third_elec}(:,:,bnd);
 for trial = 1:size(cond_2_event_histogram, 3) % loop thru trials
            log_vector = logical(temp(trial,:)) ;
            vector_value = zeros(1,length(log_vector));
            vector_value(log_vector) = trial;
            vector_value(~log_vector) = -5;
            plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))), vector_value(51:end-50), '.r')
            xlim([-.5+.05 .5-.05])
            ylim([0 size(cond_2_event_histogram, 3)])
            hold on
            plot(zeros(1,length(freq)+1), 0:length(freq),'Color','k','LineWidth',1)
            xlabel('time ')
            ylabel('trials')
            hold on
        end
title([band_name{5} ' band: ' band_range{bnd}])
set(gca, 'FontSize',12,'FontWeight','bold')
subplot(2,2,4)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothb,'r','LineWidth',2.5)
xlim([-.5+.05 .5-.05])
ylim([0 mx])
set(gca, 'FontSize',12,'FontWeight','bold')

figure
sum_tracea = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('lo gamma',band_name))),1);
sum_traceb = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('hi gamma',band_name))),1);
sum_trace_smootha2 = conv(sum_tracea, ones(1,win_size)/win_size, 'same')
sum_trace_smoothb2 = conv(sum_traceb, ones(1,win_size)/win_size, 'same')
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smootha,  'g', 'LineWidth',2.5)
hold on
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothb,  'r','LineWidth',2.5)
plot(zeros(1,round(max([sum_trace_smootha2 sum_trace_smoothb2 ]))+1), 0:round(max([sum_trace_smootha2 sum_trace_smoothb2 ])),'Color','k','LineWidth',1)
xlim([-.5+.05 .5-.05])
ylim([0 round(max([sum_trace_smootha2 sum_trace_smoothb2 ]))])
legend({'low gamma', 'high gamma'})
title('InSeq Corr relative to Poke-In')
set(gca, 'FontSize',12,'FontWeight','bold')
%% low and hi gamma in all elecs
figure
for elec = 1:4
    third_elec = direc_chan(elec)

sum_tracea = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('lo gamma',band_name))),1);
sum_traceb = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('hi gamma',band_name))),1);
sum_tracec = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('ripple',band_name))),1);

sum_trace_smootha2 = conv(sum_tracea, ones(1,win_size)/win_size, 'same')
sum_trace_smoothb2 = conv(sum_traceb, ones(1,win_size)/win_size, 'same')
sum_trace_smoothc2 = conv(sum_tracec, ones(1,win_size)/win_size, 'same')

subplot(1,4,elec)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smootha2,  'g', 'LineWidth',3)
hold on
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothb2,  'r','LineWidth',3)
%plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothc2,  'b','LineWidth',3)
plot(zeros(1,ceil(max([sum_trace_smootha2 sum_trace_smoothb2 ]))+1), 0:ceil(max([sum_trace_smootha2 sum_trace_smoothb2 ])),'Color','k','LineWidth',1)

xlim([-.5+.05 .5-.05])
ylim([0 ceil(max([sum_trace_smootha2 sum_trace_smoothb2 ]))])
legend({'low gamma', 'high gamma'})
title(['InSeq Corr relative to Poke-In : chan' num2str(direc_chan(elec))])
set(gca, 'FontSize',12,'FontWeight','bold')
end


%%
figure
for elec = 1:4
third_elec = direc_chan(elec)
sum_tracea = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('lo gamma',band_name))),1);
sum_traceb = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('hi gamma',band_name))),1);
sum_tracec = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('ripple',band_name))),1);
sum_trace_smootha2 = conv(sum_tracea, ones(1,win_size)/win_size, 'same')
sum_trace_smoothb2 = conv(sum_traceb, ones(1,win_size)/win_size, 'same')
sum_trace_smoothc2 = conv(sum_tracec, ones(1,win_size)/win_size, 'same')

subplot(1,4,elec)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smootha2,  'g', 'LineWidth',3)
hold on
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothb2,  'r','LineWidth',3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smoothc2,  'b','LineWidth',3)
plot(zeros(1,ceil(max([sum_trace_smootha2 sum_trace_smoothb2 sum_trace_smoothc2]))+1), 0:ceil(max([sum_trace_smootha2 sum_trace_smoothb2 sum_trace_smoothc2])),'Color','k','LineWidth',1)

xlim([-.5+.05 .5-.05])
ylim([0 ceil(max([sum_trace_smootha2 sum_trace_smoothb2 sum_trace_smoothc2]))])

legend({'low gamma', 'high gamma', 'ripple'})
title(['InSeq Corr relative to Poke-In : chan' num2str(direc_chan(elec))])
set(gca, 'FontSize',12,'FontWeight','bold')
end
%%
figure
for elec = 1:4
third_elec = direc_chan(elec)
sum_tracea = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('theta1',band_name))),1);
sum_traceb = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('theta2',band_name))),1);
sum_tracec = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('alpha',band_name))),1);
sum_traced = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('beta',band_name))),1);
sum_tracee = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('lo gamma',band_name))),1);
sum_tracef = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('hi gamma',band_name))),1);
sum_traceg = sum(chan_trial_2_data{third_elec}(:,51:end-50,find(strcmp('ripple',band_name))),1);

sum_trace_smooth1 = conv(sum_tracea, ones(1,win_size)/win_size, 'same');
sum_trace_smooth2 = conv(sum_traceb, ones(1,win_size)/win_size, 'same');
sum_trace_smooth3 = conv(sum_tracec, ones(1,win_size)/win_size, 'same');
sum_trace_smooth4 = conv(sum_traced, ones(1,win_size)/win_size, 'same');
sum_trace_smooth5 = conv(sum_tracee, ones(1,win_size)/win_size, 'same');
sum_trace_smooth6 = conv(sum_tracef, ones(1,win_size)/win_size, 'same');
sum_trace_smooth7 = conv(sum_traceg, ones(1,win_size)/win_size, 'same');

subplot(1,4,elec)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smooth1,  '--k', 'LineWidth',3)
hold on
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smooth2,  'k','LineWidth',3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smooth3,  'y','LineWidth',3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smooth4,  'm','LineWidth',3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smooth5,  'g','LineWidth',3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smooth6,  'r','LineWidth',3)
plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum_trace_smooth7,  'b','LineWidth',3)

all_vals = [sum_trace_smooth1 sum_trace_smooth2 sum_trace_smooth3 sum_trace_smooth4 sum_trace_smooth5 sum_trace_smooth6 sum_trace_smooth7];
plot(zeros(1,ceil(max(all_vals))+1), 0:ceil(max(all_vals)),'Color','k','LineWidth',1)

xlim([-.5+.05 .5-.05])
ylim([0 ceil(max(all_vals))])
legend({'theta1', 'theta2', 'alpha', 'beta', 'low gamma', 'high gamma', 'ripple'})
title(['InSeq Corr relative to Poke-In : chan' num2str(direc_chan(elec))])
set(gca, 'FontSize',12,'FontWeight','bold')
end


%%
i = 0;
%subplot_val = [1 5 9 13 17 21; 2 6 10 14 18 22; 3 7 11 15 19 23; 4 8 12 16 20 24]
subplot_val = [1 5 9 13 17 21 25 29 33 37 41 45; 2 6 10 14 18 22 26 30 34 38 42 46; 3 7 11 15 19 23 27 31 35 39 43 47; 4 8 12 16 20 24 28 32 36 40 44 48]
figure
for elec = direc_chan
    i = 1+i;
    for bnd = 1:6
        if bnd == 1
            bnd_column = 1
        else
            bnd_column = (bnd*2)-1
        end
        subplot(12,ca1_elec_length,subplot_val(i, bnd_column))
        temp = chan_trial_1_data{elec}(:,:,bnd);
        for trial = 1:size(cond_2_event_histogram, 3) % loop thru trials
            log_vector = logical(temp(trial,:)) ;
            vector_value = zeros(1,length(log_vector));
            vector_value(log_vector) = trial;
            vector_value(~log_vector) = -5;
            
            plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))), vector_value(51:end-50), '.r')
            xlim([-.5+.05 .5-.05])
            ylim([0 size(cond_2_event_histogram, 3)])
            xlabel('time from center of maze')
            ylabel('trials')
            hold on
        end
                title([band_name{bnd} ' band: ' band_range{bnd}])
        set(gca, 'FontSize',8,'FontWeight','bold')
        subplot(12,ca1_elec_length,subplot_val(i, bnd_column+1))
                plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum(chan_trial_1_data{elec}(:,51:end-50,bnd),1),  'b')

        ylim([0 size(cond_2_event_histogram, 3)/2])
        set(gca, 'FontSize',8,'FontWeight','bold')
        
       
        clear temp
        
    end
    
end

%% subplots for odor presentation poke out
i = 0;
figure
%subplot_val = [1 5 9 13 17 21; 2 6 10 14 18 22; 3 7 11 15 19 23; 4 8 12 16 20 24]
subplot_val = [1 5 9 13 17 21 25 29 33 37 41 45; 2 6 10 14 18 22 26 30 34 38 42 46; 3 7 11 15 19 23 27 31 35 39 43 47; 4 8 12 16 20 24 28 32 36 40 44 48]

for elec = direc_chan
    i = 1+i
    
    for bnd = 1:6
        if bnd == 1
            bnd_column = 1
        else
            bnd_column = (bnd*2)-1
        end
        subplot(12,ca1_elec_length,subplot_val(i, bnd_column))
        temp = chan_trial_2_data{elec}(:,:,bnd);
        for trial = 1:size(cond_2_event_histogram, 3) % loop thru trials
            log_vector = logical(temp(trial,:)) ;
            vector_value = zeros(1,length(log_vector));
            vector_value(log_vector) = trial;
            vector_value(~log_vector) = -5;
            plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))), vector_value(51:end-50), '.r')
            xlim([-.5+.05 .5-.05])
            ylim([0 size(cond_2_event_histogram, 3)])
            hold on
            plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
            xlabel('time ')
            ylabel('trials')
            hold on
        end
        title([band_name{bnd} ' band: ' band_range{bnd}])
        set(gca, 'FontSize',8,'FontWeight','bold')
        subplot(12,ca1_elec_length,subplot_val(i, bnd_column+1))
        plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum(chan_trial_2_data{elec}(:,51:end-50,bnd),1),  'b')
        ylim([0 size(cond_2_event_histogram, 3)/2])
        set(gca, 'FontSize',8,'FontWeight','bold')
        
    end
    
end

%% subplots for odor presentation poke in
i = 0;
figure
%subplot_val = [1 5 9 13 17 21; 2 6 10 14 18 22; 3 7 11 15 19 23; 4 8 12 16 20 24]
subplot_val = [1 5 9 13 17 21 25 29 33 37 41 45; 2 6 10 14 18 22 26 30 34 38 42 46; 3 7 11 15 19 23 27 31 35 39 43 47; 4 8 12 16 20 24 28 32 36 40 44 48]

for elec = direc_chan
    i = 1+i
    
    for bnd = 1:6
        if bnd == 1
            bnd_column = 1
        else
            bnd_column = (bnd*2)-1
        end
        subplot(12,ca1_elec_length,subplot_val(i, bnd_column))
        temp = chan_trial_2_data{elec}(:,:,bnd);
        for trial = 1:size(cond_2_event_histogram, 3) % loop thru trials
            log_vector = logical(temp(trial,:)) ;
            vector_value = zeros(1,length(log_vector));
            vector_value(log_vector) = trial;
            vector_value(~log_vector) = -5;
            plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))), vector_value(51:end-50), '.r')
            xlim([-.5+.05 .5-.05])
            ylim([0 size(cond_2_event_histogram, 3)])
            hold on
            plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
            xlabel('time ')
            ylabel('trials')
            hold on
        end
        title([band_name{bnd} ' band: ' band_range{bnd}])
        set(gca, 'FontSize',8,'FontWeight','bold')
        subplot(12,ca1_elec_length,subplot_val(i, bnd_column+1))
        plot(linspace(-.5+.05,.5-.05, length(log_vector(51:end-50))),sum(chan_trial_2_data{elec}(:,51:end-50,bnd),1),  'b')
        ylim([0 size(cond_2_event_histogram, 3)])
        set(gca, 'FontSize',8,'FontWeight','bold')
        
    end
    
end


%% make raster plots cond 1
for chan = 1:chan_counter
    figure
    for bnd = 1:6
        subplot(length(band_name),1,bnd)
        % index data per chan and per band
        temp = chan_trial_2_data{chan}(:,:,bnd);
        
        for trial = 1:size(temp,1) % loop thru trials
            log_vector = logical(temp(trial,:)) ;
            vector_value = zeros(1,length(log_vector));
            vector_value(log_vector) = trial;
             vector_value(~log_vector) = -5;

            plot(linspace(-.5,.5, length(log_vector)), vector_value, '.r')
            xlabel('time ') 
            ylabel('trials')
            hold on
        end
        title([band_name{bnd} ' band'])
            ylim([1 size(temp,1)])
            set(gca, 'FontSize',8,'FontWeight','bold')

    end
    
    %enlarge to fullscreen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % save
    saveas(gcf,['run_states_chan' num2str(chan_length(chan)) '.png'])

end

%% make raster plots cond 2
for chan = 1:chan_counter
    figure
    for bnd = 1:6
        subplot(length(band_name),1,bnd)
        % index data per chan and per band
        temp = chan_trial_2_data{chan}(:,:,bnd);
        
        for trial = 1:size(temp,1) % loop thru trials
            log_vector = logical(temp(trial,:)) ;
            vector_value = zeros(1,length(log_vector));
            vector_value(log_vector) = trial;
            vector_value(~log_vector) = -5;

            plot(linspace(-.5,.5, length(log_vector)), vector_value, '.r')
            xlabel('time ')
            ylabel('trials')
            hold on
        end
        title([band_name{bnd} ' band'])
            ylim([1 size(temp,1)])
            set(gca, 'FontSize',8,'FontWeight','bold')

    end
    
    %enlarge to fullscreen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % save
    saveas(gcf,['run_states_chan' num2str(chan_length(chan)) 'B.png'])
close all
end


