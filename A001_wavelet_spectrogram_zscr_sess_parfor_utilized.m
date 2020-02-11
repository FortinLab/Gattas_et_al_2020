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

%%
%old way of getting times
% if strcmp('novel1', task)
%     x = [1.1 .8 1 1.1];
%  response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);
% 
% elseif strcmp('novel2', task)
%     x = [1 1 1 1.18];
%  response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);
% 
% elseif strcmp('welltrained', task)
%     x = [1.1 1.2 1 1.1 1.1 .8];
%  response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);
% end
%figure; plot(trial_num(inSeqLog), times(inSeqLog), 'r*') % look at response distribution
%hold on
%plot(trial_num(otSeqLog), times(otSeqLog), 'b*')
% legend('inseq', 'outseq')
% set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% title(['anim: '  num2str(anim) ' ' task])
%figure; plot( times, 'r*') % look at response distribution
% line([0 length(trial_num)], [x(anim), x(anim)], 'color', 'k', 'LineWidth', 2)
% ylim([0 3])




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
indiv_freq_inseq_trls_chans = indiv_freq_inseq_trls_chans(:, 101:end-100);
indiv_freq_otseq_trls_chans = indiv_freq_otseq_trls_chans(:, 101:end-100);
absolute_indiv_freq_inseq_trls_chans = absolute_indiv_freq_inseq_trls_chans(:, 101:end-100);
absolute_indiv_freq_otseq_trls_chans = absolute_indiv_freq_otseq_trls_chans(:, 101:end-100);

%% plot indiv freq
figure
mn = -1
mx = .6

% plot first condition
tickmarks = 1:30:length(freq);
subplot(2,1,1)
if strcmp('p_i', lock)
    imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq), indiv_freq_inseq_trls_chans)
elseif strcmp('p_o', lock) || strcmp('FrontReward', lock)
    imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq), indiv_freq_inseq_trls_chans)
end
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
%caxis([min(min([indiv_freq_inseq_trls_chans indiv_freq_otseq_trls_chans])) per*max(max(([indiv_freq_inseq_trls_chans indiv_freq_otseq_trls_chans])))])
caxis([mn mx])

xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)

% add appropriate figure title
if strcmp('corr_incorr_excld_A', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_1) ' corr trials excld A '])
elseif strcmp('corr_incorr_incld_A', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_1) ' corr trials incld A '])
    
elseif strcmp('inseq_outseq_con', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_1) ' incorr inseq trials exld A '])
elseif strcmp('inseqcorr_otseqincorr', trial_selection) || strcmp('inseqcorr_inseqincorr' , trial_selection) ||...
        strcmp('inseq_outseq_exp', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_1) ' corr inseq trials exld A '])
elseif strcmp('otseqcorr_otseqincorr', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_1) ' corr otseq trials exld A '])
end
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
hold on
line([0, 1.5], [120 120], 'color', 'r') % 19 hz
line([0, 1.5], [91 91], 'color', 'r') % 36 hz

% plot second condition
subplot(2,1,2)
if strcmp('p_i', lock)
    imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq), indiv_freq_otseq_trls_chans)
elseif strcmp('p_o', lock)|| strcmp('FrontReward', lock)
    imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq), indiv_freq_otseq_trls_chans)
end
set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
% add appropriate figure title
if strcmp('corr_incorr_excld_A', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_2) ' incorr trials '])
elseif strcmp('corr_incorr_incld_A', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_2) ' incorr trials incld A '])
elseif strcmp('inseq_outseq_exp', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_2) ' corr otseq trials exld A '])
elseif strcmp('inseq_outseq_con', trial_selection) || strcmp('inseqcorr_otseqincorr', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_2) ' incorr otseq trials exld A '])
elseif strcmp('inseqcorr_inseqincorr', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_2) ' incorr inseq trials exld A '])
elseif strcmp('otseqcorr_otseqincorr', trial_selection)
    title(['spectrogram across chans for ' num2str(total_tiral_log_2) ' incorr otseq trials exld A '])
end
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

hold on
line([0, 1.5], [120 120], 'color', 'r') % 19 hz
line([0, 1.5], [91 91], 'color', 'r') % 36 hz


%% plot indiv elecs
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;
mx=.8

hold on
line([0, 1.5], [120 120], 'color', 'r') % 19 hz
line([0, 1.5], [91 91], 'color', 'r') % 36 hz

for elec = direc_chan
        figure
        subplot(2,1,1)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,:,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,:,elec))
        end
        
        set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn mx])
        xlabel('time')
        ylabel('freq')
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
hold on
line([0, 1.5], [120 120], 'color', 'r') % 19 hz
line([0, 1.5], [91 91], 'color', 'r') % 36 hz

        % add appropriate figure title
        if strcmp('corr_incorr_excld_A', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_1(elec)) ' corr trials excld A '])
        elseif strcmp('corr_incorr_incld_A', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_1(elec)) ' corr trials incld A '])
            
        elseif strcmp('inseq_outseq_con', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_1(elec)) ' incorr inseq trials exld A '])
        elseif strcmp('inseqcorr_otseqincorr', trial_selection) || strcmp('inseqcorr_inseqincorr' , trial_selection) ||...
                strcmp('inseq_outseq_exp', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: '  num2str(chan_trials_log_1(elec)) ' corr inseq trials exld A '])
        elseif strcmp('otseqcorr_otseqincorr', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: '  num2str(chan_trials_log_1(elec)) ' corr otseq trials exld A '])
        end
        set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    
        
        subplot(2,1,2)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,elec))
        end
        set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn mx])
        xlabel('time')
        ylabel('freq')
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
hold on
line([0, 1.5], [120 120], 'color', 'r') % 19 hz
line([0, 1.5], [91 91], 'color', 'r') % 36 hz

        % add appropriate figure title
        if strcmp('corr_incorr_excld_A', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_2(elec)) ' incorr trials '])
        elseif strcmp('corr_incorr_incld_A', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_2(elec)) ' incorr trials incld A '])
        elseif strcmp('inseq_outseq_exp', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_2(elec)) ' corr otseq trials exld A '])
        elseif strcmp('inseq_outseq_con', trial_selection) || strcmp('inseqcorr_otseqincorr', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_2(elec)) ' incorr otseq trials exld A '])
        elseif strcmp('inseqcorr_inseqincorr', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_2(elec)) ' incorr inseq trials exld A '])
        elseif strcmp('otseqcorr_otseqincorr', trial_selection)
            title(['chan ' num2str(chan_length(elec)) ' spectrogram: ' num2str(chan_trials_log_2(elec)) ' incorr otseq trials exld A '])
        end
        set(gca, 'FontSize', 14, 'FontWeight', 'bold')
      % saveas(gcf, [trial_selection 'elec' num2str(elec) '.png']) 
    end
    



%% spatial map indications for each anim
[medial_elecs, lateral_elecs, i, i2,...
 subplot_row_lateral, subplot_row_medial , ...
 subplot_column_lateral, subplot_column_medial] = electrode_map( anim );

mn = -.4
mx =1.5

% plot first condition
tickmarks = 1:30:size(indiv_freq_inseq_trls_chans,1);
i_counter = 0;
figure
for elec = lateral_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_lateral(anim), subplot_column_lateral(anim), i2(i_counter))
    if strcmp('p_i', lock)
        imagesc(linspace(-.5+.0017,1.5-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,101:end-100,find(chan_length==elec)))
    elseif strcmp('p_o', lock)
        imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,101:end-100,find(chan_length==elec)))
    end
    set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title(['chan ' num2str( lateral_elecs( i_counter)) ' : ' num2str(chan_trials_log_1(find(chan_length==elec))) ' trials'])

end

i_counter = 0;
figure
for elec = medial_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_medial(anim), subplot_column_medial(anim), i(i_counter))
    if strcmp('p_i', lock)
        
        imagesc(linspace(-.5+.0017,1.5-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,101:end-100,find(chan_length==elec)))
    elseif strcmp('p_o', lock)
        imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,101:end-100,find(chan_length==elec)))
    end
     set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    
    
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title(['chan ' num2str( medial_elecs( i_counter)) ' : ' num2str(chan_trials_log_1(find(chan_length==elec))) ' trials'])

end



%% second condition electrode map plots

i_counter = 0;
figure
for elec = lateral_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_lateral(anim), subplot_column_lateral(anim), i2(i_counter))
    if strcmp('p_i', lock)
        imagesc(linspace(-.5+.0017,1.5-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,find(chan_length==elec)))
    elseif strcmp('p_o', lock)
        imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,find(chan_length==elec)))
    end
     set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title(['chan ' num2str( lateral_elecs( i_counter)) ' : ' num2str(chan_trials_log_2(find(chan_length==elec))) ' trials'])

end

figure
i_counter = 0;

for elec = medial_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_medial(anim), subplot_column(anim), i(i_counter))
    if strcmp('p_i', lock)
        
        imagesc(linspace(-.5+.0017,1.5-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,find(chan_length==elec)))
    elseif strcmp('p_o', lock)
        imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,find(chan_length==elec)))
    end
     set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title(['chan ' num2str( medial_elecs( i_counter)) ' : ' num2str(chan_trials_log_2(find(chan_length==elec))) ' trials'])

end


%%
% 
% % unit activity
% chan_name = 'SuperChris-2-12-09_MS_T'
% chan_units      = [];
% unit_activity   = cell(chan_counter,1);
% for chan = 1:chan_counter
%     load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
%     [tetStart, tetEnd] = regexp(statMatrixColIDs, 'T([0-9]*)-U([0-9]*)');
%     column_idx = find(~cellfun(@isempty,tetStart));
%     chan_units = [chan_units; [chan_length(chan) length(column_idx)]];
%     if ~ isempty(column_idx)
%     unit_activity{chan} = [sum(statMatrix(:, column_idx))];
%     end
% end