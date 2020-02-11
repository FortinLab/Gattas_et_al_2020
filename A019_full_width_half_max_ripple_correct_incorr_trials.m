% first run A001
mn_acrs_trials_inseq (:, :, elec)  = nanmean(norm_freq_acrs_chan_cond_1(:, :, :, elec), 3);
mn_acrs_trials_outseq (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_2(:, :, :, elec), 3);
cond_1 = mn_acrs_trials_inseq;
cond_2 = mn_acrs_trials_outseq;
%% calc FWHM at group elec 
win_size = 40;
strt = 100
stp  = 500


temp = nanmean(cond_1,3);
ripple_trace = mean(temp(freq>150 & freq<250, strt:stp),1);
low_gamma_trace = mean(temp(freq>25 & freq<55, strt:stp),1);
hi_gamma_trace = mean(temp(freq>60 & freq<100, strt:stp),1);

smooth_ripple_trace = conv(ripple_trace, ones(1,win_size)/win_size, 'same');
smooth_low_gamma_trace = conv(low_gamma_trace, ones(1,win_size)/win_size, 'same');
smooth_hi_gamma_trace = conv(hi_gamma_trace, ones(1,win_size)/win_size, 'same');


smooth_ripple_trace1 = smooth_ripple_trace;
[fullwidth_at_halfmax,I_peak, half_max_idx_a, half_max_idx_b] = fwhm_SG(smooth_ripple_trace)
%figure;plot(smooth_ripple_trace);
% hold on
% plot(I_peak, smooth_ripple_trace(I_peak), '*r', 'LineWidth',5)
% plot(half_max_idx_a, smooth_ripple_trace(half_max_idx_a), '*r', 'LineWidth',5)
% plot(half_max_idx_b, smooth_ripple_trace(half_max_idx_b), '*r', 'LineWidth',5)
clear temp smooth_ripple_trace

temp = nanmean(cond_2,3);
ripple_trace = mean(temp(freq>150 & freq<250, strt:stp),1);
smooth_ripple_trace = conv(ripple_trace, ones(1,win_size)/win_size, 'same');
smooth_ripple_trace2 = smooth_ripple_trace;
low_gamma_trace2 = mean(temp(freq>25 & freq<55, strt:stp),1);
hi_gamma_trace2 = mean(temp(freq>60 & freq<100, strt:stp),1);
smooth_low_gamma_trace2 = conv(low_gamma_trace2, ones(1,win_size)/win_size, 'same');
smooth_hi_gamma_trace2 = conv(hi_gamma_trace2, ones(1,win_size)/win_size, 'same');

[fullwidth_at_halfmax,I_peak, half_max_idx_a, half_max_idx_b] = fwhm_SG(smooth_ripple_trace)
% figure;plot(smooth_ripple_trace);
% hold on
% plot(I_peak, smooth_ripple_trace(I_peak), '*r', 'LineWidth',5)

figure;
plot(smooth_ripple_trace1, 'g', 'LineWidth',5)
hold on
plot(smooth_ripple_trace2, 'k', 'LineWidth',5)
plot(low_gamma_trace, 'm', 'LineWidth',5)
plot(low_gamma_trace2, 'b', 'LineWidth',5)
plot(hi_gamma_trace, 'r', 'LineWidth',5)
plot(hi_gamma_trace2, 'y', 'LineWidth',5)


legend({'corr ripple', 'incorr', 'corr lo gamma', 'incorr', 'incorr hi gamma', 'incorr'})
title (['anim ' num2str(anim)])
xlim([0 stp])
set(gca, 'FontSize',14,'FontWeight','bold')
%% check win size
fullwidth_at_halfmax_chan= zeros(chan_counter,2)
for chan = 1:chan_counter
    figure
    temp1 = cond_1 (:, :, chan);
    ripple_trace1 = mean(temp1(freq>150 & freq<250, strt:stp),1);
    smooth_ripple_trace1 = conv(ripple_trace1, ones(1,win_size)/win_size, 'same');
    subplot(2,1,1)
    plot(smooth_ripple_trace1)
    [fullwidth_at_halfmax,I_peak, half_max_idx_a, half_max_idx_b] = fwhm_SG(smooth_ripple_trace1)
    hold on
    plot(I_peak, smooth_ripple_trace1(I_peak), '*r', 'LineWidth',5)
    plot(half_max_idx_a, smooth_ripple_trace1(half_max_idx_a), '*r', 'LineWidth',5)
    plot(half_max_idx_b, smooth_ripple_trace1(half_max_idx_b), '*r', 'LineWidth',5)
    if I_peak>100 && I_peak<300
        fullwidth_at_halfmax_chan(chan,1) = fullwidth_at_halfmax;
    else
        fullwidth_at_halfmax_chan(chan,1) = nan;
    end
    title( ['chan ' num2str(chan) '  FWWHM value: ' num2str(fullwidth_at_halfmax_chan(chan,1))])
    clear  fullwidth_at_halfmax
    
    temp2 = cond_2 (:, :, chan);
    ripple_trace2 = mean(temp2(freq>150 & freq<250, strt:stp),1);
    smooth_ripple_trace2 = conv(ripple_trace2, ones(1,win_size)/win_size, 'same');
    subplot(2,1,2)
    plot(smooth_ripple_trace2)
    [fullwidth_at_halfmax,I_peak, half_max_idx_a, half_max_idx_b] = fwhm_SG(smooth_ripple_trace2)
    hold on
    plot(I_peak, smooth_ripple_trace2(I_peak), '*r', 'LineWidth',5)
    plot(half_max_idx_a, smooth_ripple_trace2(half_max_idx_a), '*r', 'LineWidth',5)
    plot(half_max_idx_b, smooth_ripple_trace2(half_max_idx_b), '*r', 'LineWidth',5)
    if I_peak>100 && I_peak<300
        fullwidth_at_halfmax_chan(chan,2) =fullwidth_at_halfmax;
    else
        fullwidth_at_halfmax_chan(chan,2) = nan;
    end
    title([num2str(chan) '  FWWHM value: ' num2str(fullwidth_at_halfmax_chan(chan,2))])
    clear  fullwidth_at_halfmax
end
cond1_mean = nanmean(fullwidth_at_halfmax_chan(:,1))
cond2_mean = nanmean(fullwidth_at_halfmax_chan(:,2))


%% calc FWHM foro each elec 
fullwidth_at_halfmax_chan = zeros(chan_counter,2);
temp1 = []
for  chan = [1 2 4:chan_counter]
    temp1 = cond_1 (:, :, chan);
    ripple_trace1 = mean(temp1(freq>150 & freq<250,strt:stp),1);
    smooth_ripple_trace1 = conv(ripple_trace1, ones(1,win_size)/win_size, 'same');
[fullwidth_at_halfmax,I_peak, half_max_idx_a, half_max_idx_b] = fwhm_SG(smooth_ripple_trace1)
        if I_peak>150 && I_peak<300
          fullwidth_at_halfmax_chan(chan,1) = fullwidth_at_halfmax;    
        else
           fullwidth_at_halfmax_chan(chan,1) = nan;
        end

        temp2 = cond_2 (:, :, chan);
        ripple_trace2 = mean(temp2(freq>150 & freq<250,strt:stp),1);
        smooth_ripple_trace2 = conv(ripple_trace2, ones(1,win_size)/win_size, 'same')
        [fullwidth_at_halfmax,I_peak, half_max_idx_a, half_max_idx_b] = fwhm_SG(smooth_ripple_trace2)
        if I_peak>150 && I_peak<300
            fullwidth_at_halfmax_chan(chan,2) = fullwidth_at_halfmax;
        else
            fullwidth_at_halfmax_chan(chan,2) = nan;
        end
end

cond1_mean = nanmean(fullwidth_at_halfmax_chan(:,1))
cond2_mean = nanmean(fullwidth_at_halfmax_chan(:,2))

%% save values for each animal
%anim_FWHM_ripple = zeros(6,2)
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
load('anim_FWHM_ripple.mat')
anim_FWHM_ripple(anim,:) =[cond1_mean cond2_mean] 
save('anim_FWHM_ripple', 'anim_FWHM_ripple')