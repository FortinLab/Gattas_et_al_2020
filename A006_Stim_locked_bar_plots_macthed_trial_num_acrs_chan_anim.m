clear all; close all
anim = 5
lock = 'p_o' %'p_o' % p_i: poke-in % p_o: poke-out
task = 'welltrained' %novel1, novel2, welltrained
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

% load behavioral and neural data
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load('baseline_info_wavelet_32num_10logdb_3hz_250hz_notched_artifact_reject.mat')
load('chan_artifact_thresh.mat')
load('BehaviorMatrix.mat')
load([chan_name '4.mat']) % load exmp chan

% electrodes along CA1 axis for each animal for group average
cd('D:\Gattas\ephys_data_final\group_plots')
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

% variables
fs         = 1000; % sampling rate

% define eventWindow/lock
if strcmp('p_i', lock)
    eventWindow = [-0.5 1.5];
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
% anim 3 has a 47 min outlier trials = nan
if anim==3 && strcmp('novel2', task) || anim==4 && strcmp('novel2', task)  
    times(find(times==max(times))) = nan;
end
[ cutoff ] = get_response_time_cutoff( anim, task, times,inSeqLog ,otSeqLog, trial_num)
 response_time= [pokeInAlignedBehavMatrix.PokeDuration]>cutoff;

% old approach
% figure; plot(trial_num(inSeqLog), times(inSeqLog), 'r*') % look at response distribution
% hold on
% plot(trial_num(otSeqLog), times(otSeqLog), 'b*')
% if strcmp('novel1', task)
%     x = [1.1 .8 1 1.1];
%     response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);
%     
% elseif strcmp('novel2', task)
%     x = [1 1 1 1.18];
%     response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);
%     
% elseif strcmp('welltrained', task)
%     x = [1.1 1.2 1 1.1 1.1 .8];
%     response_time= [pokeInAlignedBehavMatrix.PokeDuration]>x(anim);
% end
% line([0 length(trial_num)], [x(anim), x(anim)], 'color', 'k', 'LineWidth', 2)
% line([0 length(trial_num)], [1.2, 1.2], 'color', 'm', 'LineWidth', 2)
% 
% ylim([0 3])
% legend('inseq', 'outseq')
% set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% title(['anim: '  num2str(anim) ' ' task])
%figure; plot( times, 'r*') % look at response distribution
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

% save figure and inspect trial counts
saveas(gcf, 'responseDistrib.eps')
saveas(gcf, 'responseDistrib.jpeg')
[sum(inSeqCorrLog) sum(otSeqCorrLog) sum(inSeqWrongLog) sum(otSeqWrongLog)]

%% wavelet params
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
norm_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_1), length(direc_chan));
norm_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_2), length(direc_chan));
norm_freq_acrs_chan_cond_3 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_3), length(direc_chan));
norm_freq_acrs_chan_cond_4 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), sum(trial_log_4), length(direc_chan));
trial_log_1_artifact = cell(length(direc_chan),1);
trial_log_2_artifact = cell(length(direc_chan),1);
trial_log_3_artifact = cell(length(direc_chan),1);
trial_log_4_artifact = cell(length(direc_chan),1);


tic
counter = 0;
for chan =direc_chan % loop thru chans
    counter = counter+1;
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
        cwt_power = 10*log10(abs(cwt.cfs).^2);
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
        norm_freq_acrs_chan_cond_1(:, :, trial, counter) =  norm_freq;
        
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
        cwt_power = 10*log10(abs(cwt.cfs).^2);
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
        norm_freq_acrs_chan_cond_2(:, :, trial, counter) =  norm_freq;
        
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
        cwt_power = 10*log10(abs(cwt.cfs).^2);
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
        norm_freq_acrs_chan_cond_3(:, :, trial, counter) =  norm_freq;
        
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
        cwt_power = 10*log10(abs(cwt.cfs).^2);
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
        norm_freq_acrs_chan_cond_4(:, :, trial, counter) =  norm_freq;
        
    end
    
    trial_log_4_artifact{chan}=  trial_log_4_artifact_temp';
    
end
toc


% save non-downsampled trial data per animal
cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
save(['anim' num2str(anim) '_spectrogram_data_' lock '_' task], 'norm_freq_acrs_chan_cond_1', 'norm_freq_acrs_chan_cond_2', 'norm_freq_acrs_chan_cond_3','norm_freq_acrs_chan_cond_4')
beep
%% find cond with min val
counter_4_1 = 0;
counter_4_2 = 0;
counter_4_3 = 0;
counter_4_4 = 0;

counter_3_1 = 0;
counter_3_2 = 0;
counter_3_3 = 0;
counter_3_4 = 0;

counter_2_1 = 0;
counter_2_2 = 0;
counter_2_3 = 0;
counter_2_4 = 0;

counter_1_1 = 0;
counter_1_2 = 0;
counter_1_3 = 0;
counter_1_4 = 0;

trial_counter = max([sum(trial_log_1) sum(trial_log_2) sum(trial_log_3) sum(trial_log_4)]);
for chan = 1:4
    for trial = 1:trial_counter %trial_counter
        % cond 4
        if  trial<=sum(trial_log_4) & ~isnan(norm_freq_acrs_chan_cond_4(:,:, trial, chan))
            if chan ==1
                counter_4_1 = counter_4_1+1;
                chan_1_cond_4{counter_4_1} = norm_freq_acrs_chan_cond_4(:,:, trial, chan);
            elseif chan ==2
                counter_4_2 = counter_4_2+1;
                chan_2_cond_4{counter_4_2} = norm_freq_acrs_chan_cond_4(:,:, trial, chan);
            elseif chan ==3
               counter_4_3 = counter_4_3+1;

                chan_3_cond_4{counter_4_3} = norm_freq_acrs_chan_cond_4(:,:, trial, chan);
            else
                counter_4_4 = counter_4_4+1;
                chan_4_cond_4{counter_4_4} = norm_freq_acrs_chan_cond_4(:,:, trial, chan);
            end
        end
        
        % cond 3
        if  trial<=sum(trial_log_3) & ~isnan(norm_freq_acrs_chan_cond_3(:,:, trial, chan))
            if chan ==1
                counter_3_1 = counter_3_1+1;
                chan_1_cond_3{counter_3_1} = norm_freq_acrs_chan_cond_3(:,:, trial, chan);
            elseif chan ==2
                counter_3_2 = counter_3_2+1;
                chan_2_cond_3{counter_3_2} = norm_freq_acrs_chan_cond_3(:,:, trial, chan);
            elseif chan ==3
                counter_3_3 = counter_3_3+1;
                chan_3_cond_3{counter_3_3} = norm_freq_acrs_chan_cond_3(:,:, trial, chan);
            else
                counter_3_4 = counter_3_4+1;
                chan_4_cond_3{counter_3_4} = norm_freq_acrs_chan_cond_3(:,:, trial, chan);
            end
        end        
                
        % cond 2
        if  trial<=sum(trial_log_2) & ~isnan(norm_freq_acrs_chan_cond_2(:,:, trial, chan))
            if chan ==1
                 counter_2_1 = counter_2_1+1;
                %cond_chan_
                chan_1_cond_2{chan,counter_2_1} = norm_freq_acrs_chan_cond_2(:,:, trial, chan);
            elseif chan ==2
                counter_2_2 = counter_2_2+1;
                chan_2_cond_2{counter_2_2} = norm_freq_acrs_chan_cond_2(:,:, trial, chan);
            elseif chan ==3
                counter_2_3 = counter_2_3+1;
                chan_3_cond_2{counter_2_3} = norm_freq_acrs_chan_cond_2(:,:, trial, chan);
            else
                counter_2_4 = counter_2_4+1;
                chan_4_cond_2{counter_2_4} = norm_freq_acrs_chan_cond_2(:,:, trial, chan);
            end
        end
        
        % cond 1
        if  trial<=sum(trial_log_1) & ~isnan(norm_freq_acrs_chan_cond_1(:,:, trial, chan))            
            if chan ==1
                counter_1_1 = counter_1_1+1;
                %cond_chan_
                chan_1_cond_1{counter_1_1} = norm_freq_acrs_chan_cond_1(:,:, trial, chan);
            elseif chan ==2
                counter_1_2 = counter_1_2+1;                
                chan_2_cond_1{counter_1_2} = norm_freq_acrs_chan_cond_1(:,:, trial, chan);
            elseif chan ==3
                counter_1_3 = counter_1_3+1;                
                chan_3_cond_1{counter_1_3} = norm_freq_acrs_chan_cond_1(:,:, trial, chan);
            else
                counter_1_4 = counter_1_4+1;                
                chan_4_cond_1{counter_1_4} = norm_freq_acrs_chan_cond_1(:,:, trial, chan);
            end
        end        
    end       
end

% min trial num
if anim ==1
    cond_1_num = [size(chan_1_cond_1,2) size(chan_2_cond_1,2) size(chan_3_cond_1,2) size(chan_4_cond_1,2) ]
    cond_2_num = [size(chan_1_cond_2,2) size(chan_2_cond_2,2) size(chan_3_cond_2,2) size(chan_4_cond_2,2) ]
    cond_4_num = [size(chan_1_cond_4,2) size(chan_2_cond_4,2) size(chan_3_cond_4,2) size(chan_4_cond_4,2) ]
else
    cond_1_num = [size(chan_1_cond_1,2) size(chan_2_cond_1,2) size(chan_3_cond_1,2) size(chan_4_cond_1,2) ]
    cond_2_num = [size(chan_1_cond_2,2) size(chan_2_cond_2,2) size(chan_3_cond_2,2) size(chan_4_cond_2,2) ]
    cond_3_num = [size(chan_1_cond_3,2) size(chan_2_cond_3,2) size(chan_3_cond_3,2) size(chan_4_cond_3,2) ]
    cond_4_num = [size(chan_1_cond_4,2) size(chan_2_cond_4,2) size(chan_3_cond_4,2) size(chan_4_cond_4,2) ]
end
beep
%% save matched num
if strcmp('novel1', task)
    if strcmp('p_o', lock)
        cond1 = [57 58 78 69]
        cond2 = [1   0 0 3]
        cond3 = [2   1 11 9]
        cond4 = [13  8 8 10]
    else
        cond1 = []
        cond2 = []
        cond3 = []
        cond4 = []
    end
elseif strcmp('novel2', task)
        if strcmp('p_o', lock)
        cond1 = []
        cond2 = []
        cond3 = []
        cond4 = []
    else
        cond1 = []
        cond2 = []
        cond3 = []
        cond4 = []
    end
elseif strcmp('welltrained', task)
    if strcmp('p_o', lock)
        cond1 = [9 9 9 8 8 8 ]
        cond2 = [8 8 8 9 9 9 ]
        cond3 = [0 17 1 6 17 10]
        cond4 = [3 5 8 15 14 6]
    else
        cond1 = [8 9 8 8 8 8]
        cond2 = [9 8 8 8 8 8]
        cond3 = [0 15 1 7 16 10]
        cond4 = [2 4 9 15 13 6]
    end
end

%%
for trial = 1:cond1(anim)
    cond_1(:,:,trial,1)= chan_1_cond_1{trial};
    cond_1(:,:,trial,2)= chan_2_cond_1{trial};
    cond_1(:,:,trial,3)= chan_3_cond_1{trial};
    cond_1(:,:,trial, 4)= chan_4_cond_1{trial};
end
for trial = 1:cond2(anim)
    cond_2(:,:,trial,1) = chan_1_cond_2{trial};
    cond_2(:,:,trial,2) = chan_2_cond_2{trial};
    cond_2(:,:,trial,3) = chan_3_cond_2{trial};
    cond_2(:,:,trial,4) = chan_4_cond_2{trial};
end

for trial = 1:cond3(anim)
    cond_3(:,:,trial,1) = chan_1_cond_3{trial};
    cond_3(:,:,trial,2) = chan_2_cond_3{trial};
    cond_3(:,:,trial,3) = chan_3_cond_3{trial};
    cond_3(:,:,trial,4) = chan_4_cond_3{trial};
end
for trial = 1:cond4(anim)
    cond_4(:,:,trial,1) = chan_1_cond_4{trial};
    cond_4(:,:,trial,2) = chan_2_cond_4{trial};
    cond_4(:,:,trial,3) = chan_3_cond_4{trial};
    cond_4(:,:,trial,4) = chan_4_cond_4{trial};
end
%%
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
if strcmp('welltrained', task)
    if anim ==1
        save(['anim' num2str(anim) '_exact_match_' lock], 'cond_1', 'cond_2', 'cond_4')
    else
        save(['anim' num2str(anim) '_exact_match_' lock], 'cond_1', 'cond_2', 'cond_3','cond_4')
    end
else
    save(['anim' num2str(anim) '_exact_match_' lock '_' task], 'cond_1', 'cond_2', 'cond_3','cond_4')
    
end