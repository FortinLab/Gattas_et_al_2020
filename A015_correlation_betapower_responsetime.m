% correl between beta power and response time
clear all; close all;clc
lock ='p_o';
task ='welltrained' ;

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

%%
close all;clc
corr_mtx = zeros(5,4);
p_mtx    = zeros(5,4);

for anim = 5
    % get InSeq poke durations
    cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])
    load('BehaviorMatrix.mat')
    eventWindow = [-0.5 0.5];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut')
    times          = [pokeInAlignedBehavMatrix.PokeDuration];
    inSeqLog       = [pokeInAlignedBehavMatrix.TranspositionDistance]==0;
    odor_excld_A = [pokeInAlignedBehavMatrix.Odor]~=1;
    corrTrlLog   = [pokeInAlignedBehavMatrix.Performance]==1;
    trial_num = 1:length(inSeqLog);
    otSeqLog = [pokeInAlignedBehavMatrix.TranspositionDistance]~=0;

    [ cutoff ] = get_response_time_cutoff( anim, task, times,inSeqLog ,otSeqLog, trial_num);
    response_time= [pokeInAlignedBehavMatrix.PokeDuration]>cutoff;
    inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A&response_time==1;
    cond_responses = times(inSeqCorrLog);

    % get data/beta power
    cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
    load(['anim' num2str(anim) '_spectrogram_data_' lock '_' task])
    figure;hold on
    % for inseq corr cond
    for chan = 1:4 % loop thru chans
        inseq_corr_beta_pwr = norm_freq_acrs_chan_cond_1((freq>20 & freq<40), 250:500, :, chan);
        
        inseq_corr_beta_pwr_avg = zeros(size(inseq_corr_beta_pwr,3), 1);
        for trials = 1:size(inseq_corr_beta_pwr,3)
            inseq_corr_beta_pwr_avg(trials,1) = mean(mean(inseq_corr_beta_pwr(:, :, trials), 1),2);
        end
        
        cond1_responses = cond_responses(cond_responses < prctile(cond_responses, 95)& cond_responses > prctile(cond_responses, 5));
        inseq_corr_beta_pwr_avg = inseq_corr_beta_pwr_avg(cond_responses < prctile(cond_responses, 95)& cond_responses > prctile(cond_responses, 5));
        [r,p] = corrcoef(inseq_corr_beta_pwr_avg(~isnan(inseq_corr_beta_pwr_avg)), cond1_responses(~isnan(inseq_corr_beta_pwr_avg)));
        corr_mtx(anim,chan) = r(2);
        p_mtx(anim,chan)    = p(2);
        plot(cond1_responses, inseq_corr_beta_pwr_avg, '*')
        xlabel('poke duration')
        ylabel('beta power')
     
    end
    title(num2str(anim))
   % legend({'elec 1', 'elec 2', 'elec 3', 'elec 4'})
   ylim([-1.5 2])
   xlim([[1.2 1.7]])
end

%%
% for otseq incorr cond

for chan = 1:chan_counter
inseq_corr_beta_pwr = norm_freq_acrs_chan_cond_2((freq>19 & freq<36),  101:end-100, :, chan);

inseq_corr_beta_pwr_avg = zeros(size(inseq_corr_beta_pwr,3), 1);
for trials = 1:size(inseq_corr_beta_pwr,3)
    inseq_corr_beta_pwr_avg(trials,1) = mean(mean(inseq_corr_beta_pwr(:,  :, trials), 1),2);
end

responses =  [pokeInAlignedBehavMatrix.PokeDuration];
cond_responses = responses(trial_log_2);
cond2_responses = cond_responses; % wont remove trials bc too small a number to begin with
%cond2_responses = cond_responses(cond_responses < prctile(cond_responses, 95)& cond_responses > prctile(cond_responses, 5));
%inseq_corr_beta_pwr_avg = inseq_corr_beta_pwr_avg(cond_responses < prctile(cond_responses, 95)& cond_responses > prctile(cond_responses, 5));

[r,p] = corrcoef(inseq_corr_beta_pwr_avg(~isnan(inseq_corr_beta_pwr_avg)), cond2_responses(~isnan(inseq_corr_beta_pwr_avg)));
if length(r)>1
    corr_val(chan) = r(2);
    sig(chan) = p(2);
    figure; plot(cond2_responses, inseq_corr_beta_pwr_avg, '*')
    title(['chan' num2str(chan) ' corrcoeff ' num2str(r(2)) ' p-value ' num2str(p(2))])
else
    corr_val(chan) = nan;
    sig(chan) = nan;
end

end
[nanmean(corr_val(sig<=0.05))   std(corr_val(sig<=0.05), 1)  length(corr_val(sig<=0.05))/length(sig)]
