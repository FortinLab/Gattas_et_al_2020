% plot elecs medial to lateral for indiv animals
% steps: run A001, then run this, saves output matrix and figure, then generate group plots
%%
per = 1;
tickmarks = 1:30:length(freq);
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;
ca1_elec_lanim_1_response_timesanim_1_response_timesength = 4;
cd('D:\Gattas\ephys_data_final\group_plots')

% plot inseq correct

mx = 1
mn = -1

i   = 1;
figure
if strcmp('inseq_outseq_exp', trial_selection) || strcmp('inseqcorr_inseqincorr', trial_selection) || strcmp('otseqcorr_otseqincorr', trial_selection)
    for elec = direc_chan
        subplot(2,length(direc_chan),i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,101:end-100,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,101:end-100,elec))
        end
        colorbar
        caxis([mn per*mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title (['chan '  num2str(chan_length(elec)) ' iseq +' ])
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
    end
    

    %% plot outseq correct
    for elec = direc_chan
        subplot(2,length(direc_chan),i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,101:end-100,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,101:end-100,elec))
        end
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn per*mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title (['chan '   num2str(chan_length(elec)) ' outseq +' ])
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
    end
    set(gca, 'FontSize', 8, 'FontWeight', 'bold')
    if strcmp('p_i',lock)
        fig_nam = ['Animal_' num2str(anim) '_lateral_to_media_Channels_' fignam '_POKE_IN.png']
    else
        fig_nam = ['Animal_' num2str(anim) '_lateral_to_media_Channels_' fignam '_POKE_OUT.png']
    end
    saveas(gcf,fig_nam)
    
    %% plot inseq incorrect
elseif strcmp('inseq_outseq_con', trial_selection)
    i   = 1;
    % figure
    for elec = direc_chan
        subplot(2,length(direc_chan),i)
        if strcmp('p_i', lock)
            imagesc(linspace(-.5+.0017,minimum_RT_outseq-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,:,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,:,elec))
        end
        colorbar
        caxis([mn per*mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title (['chan '   num2str(chan_length(elec))  ])
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1
    end
    
    
    %%  plot outseq INcorrect
    for elec = direc_chan
        subplot(2,length(direc_chan),i)
        if strcmp('p_i', lock)
            imagesc(linspace(-.5+.0017,minimum_RT_outseq-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_outseq(:,:,elec))
        end
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn per*mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title (['chan '   num2str(chan_length(elec)) ])
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
    end
    
    set(gca, 'FontSize', 8, 'FontWeight', 'bold')
    if strcmp('p_i',lock)
        fig_nam = ['Animal_' num2str(anim) '_lateral_to_media_Channels_' fignam '_POKE_IN.png']
    else
        fig_nam = ['Animal_' num2str(anim) '_lateral_to_media_Channels_' fignam '_POKE_OUT.png']
    end
    saveas(gcf,fig_nam)
    
    %%insep corr - outseq incorr
elseif strcmp('inseqcorr_otseqincorr', trial_selection)
    for elec = direc_chan
        subplot(1,length(direc_chan),i)
        imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  mn_acrs_trials_inseq(:,:,elec)- mn_acrs_trials_outseq(:,:,elec))        
        colorbar
        colormap jet
       % caxis([mn per*mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title (['chan '   num2str(chan_length(elec)) ])
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
    end
    cd('C:\Users\FortinLab\Desktop\Example Data Sets\welltrained_group_data')
    cond1= mn_acrs_trials_inseq(:, :, direc_chan)-mn_acrs_trials_outseq(:, :, direc_chan);
    filename= ['anim' num2str(anim) '_AL_PM_' trial_selection '_' lock '_matched_trials_DiffRT'];
    save(filename, 'cond1')
    
    % running spectrogram
elseif strcmp('run', trial_selection)
       % freq_range= freq<36;
        %new_freq = freq(freq_range);
      %  tickmarks =1:15:length(new_freq);
        for elec = direc_chan
        subplot(1,length(direc_chan),i)
        imagesc(linspace(-.5+.0017,.5-.0017,size(mn_acrs_trials_run,2)), 1:length(freq),  mn_acrs_trials_run(:,101:end-100,elec))
        colorbar
        colormap jet
        caxis([mn mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title (['chan '  num2str(chan_length(elec))])
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
         
        end
end
%% save figures along axis
cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim\fig2_runVsOdor')
hgexport(gcf, ['anim' num2str(anim) '_AL_PM_' trial_selection '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
hgexport(gcf, ['anim' num2str(anim) '_AL_PM_' trial_selection '.eps'], hgexport('factorystyle'), 'Format', 'eps')

%% save figures one at a time
i = 0;
for elec = direc_chan
    figure
    i = i+1;
    imagesc(linspace(-.5+.0017,.5-.0017,size(mn_acrs_trials_run,2)), 1:length(freq),  mn_acrs_trials_run(:,101:end-100,elec))
    colorbar
    colormap jet
    caxis([mn mx])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title (['chan '  num2str(chan_length(elec))])
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    hgexport(gcf, ['anim' num2str(anim) '_AL_PM_' trial_selection '_elec_' num2str(i) '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
    hgexport(gcf, ['anim' num2str(anim) '_AL_PM_' trial_selection '_elec_' num2str(i) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
    
    
end
%% save data
if strcmp('run', trial_selection)
    cond1 = mn_acrs_trials_run(:, :, direc_chan);
    filename= ['anim' num2str(anim) '_AL_PM_' trial_selection];
    save(filename, 'cond1')
else
    cond1 = mn_acrs_trials_inseq(:, :, direc_chan);
    cond2 = mn_acrs_trials_outseq(:, :, direc_chan);
    if strcmp('yes',match)
        filename= ['anim' num2str(anim) '_AL_PM_' trial_selection '_' lock '_matched_trials_' task];
    else
        filename= ['anim' num2str(anim) '_AL_PM_' trial_selection '_' lock '_nonmatched_trials_' task];
    end
    save(filename, 'cond1', 'cond2')
end

