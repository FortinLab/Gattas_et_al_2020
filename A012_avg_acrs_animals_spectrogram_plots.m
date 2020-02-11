%% load and save data from all anims
cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
trial_selection = 'run'
%%
animal = 1:5;
if strcmp('run',trial_selection)
        cond_1_all_anim = zeros(size(cond1));
    
    %cond_1_all_anim = zeros(size(mn_acrs_trials,1), size(mn_acrs_trials,2));
elseif strcmp('novel1',task)||strcmp('novel2',task)
    cond_1_all_anim = zeros(size(cond1));
    cond_2_all_anim = zeros(size(cond1));
    e1_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
    e2_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
    e3_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
    e4_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
else
    if strcmp('yes',match)
        load('anim1_AL_PM_inseq_outseq_exp_p_o_matched_trials.mat')
    else
        load(['anim1_AL_PM_inseq_outseq_exp_' lock '_nonmatched_trials_' task '.mat'])
    end
    
    cond_1_all_anim = zeros(size(cond1));
    cond_2_all_anim = zeros(size(cond1));
    e1_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
    e2_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
    e3_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
    e4_c2 = zeros(size(cond1,1), size(cond1,2), length(animal));
end

e1_c1 = zeros(size(cond1,1), size(cond1,2), length(animal));
e2_c1 = zeros(size(cond1,1), size(cond1,2), length(animal));
e3_c1 = zeros(size(cond1,1), size(cond1,2), length(animal));
e4_c1 = zeros(size(cond1,1), size(cond1,2), length(animal));


% init 4 electrodes for each condition

for a = 1:length(animal)
    if strcmp('run',trial_selection)
        load(['anim' num2str(animal(a)) '_AL_PM_' trial_selection])
    elseif strcmp('novel1',task)
        load(['anim' num2str(animal(a)) '_AL_PM_' trial_selection '_' lock '_nonmatched_trials_' task])
    elseif strcmp('novel2',task)
      load(['anim' num2str(animal(a)) '_AL_PM_' trial_selection '_' lock '_nonmatched_trials_' task])

    else
        if strcmp('yes',match)
            load(['anim' num2str(animal(a)) '_AL_PM_' trial_selection '_' lock '_matched_trials'])
        else
            load(['anim' num2str(animal(a)) '_AL_PM_' trial_selection '_' lock '_nonmatched_trials_' task '.mat'])
        end
    end

    % cond 1
    e1_c1(:,:,a) = cond1(:,:,1); %elec1cond1
    e2_c1(:,:, a) = cond1(:,:,2);
    e3_c1(:,:,a) = cond1(:,:,3);
    e4_c1(:,:, a) = cond1(:,:,4);
    
    if ~strcmp('run',trial_selection)
    %cond 2
    e1_c2(:,:, a) = cond2(:,:,1);
    e2_c2(:,:, a) = cond2(:,:,2);
    e3_c2(:,:, a) = cond2(:,:,3);
    e4_c2(:,:, a) = cond2(:,:,4);
    end
end

% avg acrs anims
cond_1_all_anim(:,:, 1)     = nanmean(e1_c1,3);
cond_1_all_anim_std(:,:, 1) = nanstd(e1_c1,1,3);

cond_1_all_anim(:,:, 2)     = nanmean(e2_c1,3);
cond_1_all_anim_std(:,:, 2) = nanstd(e2_c1,1,3);

cond_1_all_anim(:,:, 3)     = nanmean(e3_c1,3);
cond_1_all_anim_std(:,:, 3) = nanstd(e3_c1,1,3);

cond_1_all_anim(:,:, 4)     = nanmean(e4_c1,3);
cond_1_all_anim_std(:,:, 4) = nanstd(e4_c1,1,3);

    if ~strcmp('run',trial_selection)

cond_2_all_anim(:,:, 1)     = nanmean(e1_c2,3);
cond_2_all_anim_std(:,:, 1) = nanstd(e1_c2,1,3);

cond_2_all_anim(:,:, 2)     = nanmean(e2_c2,3);
cond_2_all_anim_std(:,:, 2) = nanstd(e2_c2,1,3);

cond_2_all_anim(:,:, 3)     = nanmean(e3_c2,3);
cond_2_all_anim_std(:,:, 3) = nanstd(e3_c2,1,3);

cond_2_all_anim(:,:, 4)     = nanmean(e4_c2,3);
cond_2_all_anim_std(:,:, 4) = nanstd(e4_c2,1,3);
    end
%
tickmarks = 1:30:size(cond_1_all_anim,1);
i   = 1;
mn = -3
mx = 3
figure
%

if strcmp('inseq_outseq_exp', trial_selection)|| strcmp('inseqcorr_inseqincorr', trial_selection)|| strcmp('otseqcorr_otseqincorr', trial_selection)
    % plot inseq correct
    for elec = 1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_1_all_anim(:,101:end-100,elec))./(cond_1_all_anim_std(:,101:end-100,elec)))
     
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_1_all_anim(:,101:end-100,elec))./(cond_1_all_anim_std(:,101:end-100,elec)))

        end
        
        h=colorbar
        ylabel(h,'group mean / group std')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        caxis([mn mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        colormap jet
        i = i+1;
    end
    
    % plot outseq correct
    for elec =  1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_2_all_anim(:,101:end-100,elec))./(cond_2_all_anim_std(:,101:end-100,elec)))
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_2_all_anim(:,101:end-100,elec))./(cond_2_all_anim_std(:,101:end-100,elec)))
        end
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        h=colorbar
        ylabel(h,'group mean / group std')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        
        %caxis([round(2*mn)  round(2*mx)])
        caxis([mn mx])
        colormap jet
        
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        i = i+1;
    end
elseif strcmp('inseq_outseq_con', trial_selection)
    fig_title = 'INSEQ INCORR VS. OUTSEQ INCORR'
    
    % plot inseq incorrect
    for elec = 1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  (cond_1_all_anim(:,:,elec))./(cond_1_all_anim_std(:,:,elec)))
            mx = max(max(max([cond_1_all_anim cond_2_all_anim])));
            mn = min(min(min([cond_1_all_anim cond_2_all_anim])));
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  (cond_1_all_anim(:,:,elec))./(cond_1_all_anim_std(:,:,elec)))
            mx = max(max(max([cond_1_all_anim cond_2_all_anim])));
            mn = min(min(min([cond_1_all_anim cond_2_all_anim])));
        end
        h=colorbar
        ylabel(h,'group mean / group std')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        caxis([mn mx])
        colormap jet
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        i = i+1;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
        
    end
    
    % plot outseq correct
    for elec =  1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  (cond_2_all_anim(:,:,elec))./(cond_2_all_anim_std(:,:,elec)))
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  (cond_2_all_anim(:,:,elec))./(cond_2_all_anim_std(:,:,elec)))
        end
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'FontSize', 12, 'FontWeight', 'bold')
        h=colorbar
        ylabel(h,'group mean / group std')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        caxis([mn mx])
        colormap jet
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
        
    end
elseif strcmp('inseqcorr_outseqincorr', trial_selection)
     % plot inseq incorrect
    for elec = 1:4
        subplot(1,4,elec)
            imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  (cond_1_all_anim(:,:,elec))./(cond_1_all_anim_std(:,:,elec)))
            mx = max(max(max([cond_1_all_anim ])));
            mn = min(min(min([cond_1_all_anim ])));
        
        h=colorbar
        ylabel(h,'group mean / group std')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        caxis([-4 5])
        colormap jet
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        i = i+1;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    end
elseif strcmp('run', trial_selection)
    for elec = 1:4
        subplot(1,4,elec)
            imagesc(linspace(-.5+.0017,.5-.0017,size(cond1,2)), 1:length(freq),  (cond_1_all_anim(:,101:end-100,elec))./(cond_1_all_anim_std(:,101:end-100,elec)))
%             mx = max(max(max([cond_1_all_anim ])));
%             mn = min(min(min([cond_1_all_anim ])));
        
        h=colorbar
        ylabel(h,'group mean / group std')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        caxis([mn mx])
        colormap jet
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        i = i+1;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    end
end
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

%% make plot with out T ratio - only group average
tickmarks = 1:30:size(cond_2_all_anim,1);
i   = 1;
mn = -.8
mx = .8

figure
 for elec = 1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_1_all_anim(:,101:end-100,elec)))
     
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_1_all_anim(:,101:end-100,elec)))

        end
        
        h=colorbar
        ylabel(h,'group mean')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        caxis([mn mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        colormap jet
        i = i+1;
    end
    
    % plot outseq correct
    for elec =  1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_2_all_anim(:,101:end-100,elec)))
        elseif strcmp('p_o', lock)
            imagesc(linspace(eventWindow(1)+.0017,eventWindow(2)-.0017,size(cond_1_all_anim,2)), 1:length(freq),  (cond_2_all_anim(:,101:end-100,elec)))
        end
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        h=colorbar
        ylabel(h,'group mean ')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        
        %caxis([round(2*mn)  round(2*mx)])
        caxis([mn mx])
        colormap jet
        
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        i = i+1;
    end

%% group mean not a t ratio for running
figure
    for elec = 1:4
        subplot(1,4,elec)
            imagesc(linspace(-.5+.0017,.5-.0017,size(cond1,2)), 1:length(freq),  (cond_1_all_anim(:,101:end-100,elec)))
       
        h=colorbar
        ylabel(h,'group mean')
        set(h, 'FontSize', 12, 'FontWeight', 'bold')
        caxis([mn mx])
        colormap jet
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'FontSize', 12, 'FontWeight', 'bold')
        i = i+1;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    end
%%
per =1;
tickmarks = 1:30:size(cond_1_all_anim_std,1);
i   = 1;
figure
if strcmp('inseq_outseq_exp', trial_selection)
    fig_title = 'INSEQ CORR VS. OUTSEQ CORR'
    % plot inseq correct
    for elec = 1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(-.5+.0017,minimum_RT_outseq-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  (cond_1_all_anim_std(:,:,elec)))
            mx = max(max(max([cond_2_all_anim_std cond_1_all_anim_std])));
            mn = min(min(min([cond_2_all_anim_std cond_1_all_anim_std])));
        elseif strcmp('p_o', lock)
            imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  (cond_1_all_anim_std(:,:,elec)))
            mx = max(max(max([cond_1_all_anim_std cond_2_all_anim_std])));
            mn = min(min(min([cond_1_all_anim_std cond_2_all_anim_std])));
        end
        
        h=colorbar
        ylabel(h,'SD')
        caxis([mx mn])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
    end
    
    % plot outseq correct
    for elec =  1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(-.5+.0017,minimum_RT_outseq-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  cond_2_all_anim_std(:,:,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  cond_2_all_anim_std(:,:,elec))
        end
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mx mn])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
    end
elseif strcmp('inseq_outseq_con', trial_selection)
    fig_title = 'INSEQ INCORR VS. OUTSEQ INCORR'
    
    % plot inseq incorrect
    for elec = 1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(-.5+.0017,minimum_RT_outseq-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  cond_1_all_anim(:,:,elec))
            mx = max(max(max([cond_1_all_anim cond_2_all_anim])));
            mn = min(min(min([cond_1_all_anim cond_2_all_anim])));
        elseif strcmp('p_o', lock)
            imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  cond_1_all_anim(:,:,elec))
            mx = max(max(max([cond_1_all_anim_std cond_1_all_anim_std])));
            mn = min(min(min([cond_1_all_anim_std cond_1_all_anim_std])));
        end
        colorbar
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))    
        i = i+1;
    end
    
    % plot outseq correct
    for elec =  1:4
        subplot(2,4,i)
        if strcmp('p_i', lock)
            imagesc(linspace(-.5+.0017,minimum_RT_outseq-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  cond_2_all_anim_std(:,:,elec))
        elseif strcmp('p_o', lock)
            imagesc(linspace(-1+.0017,1-.0017,size(indiv_freq_inseq_trls_chans,2)), 1:length(freq),  cond_2_all_anim_std(:,:,elec))
        end
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn per*mx])
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        i = i+1;
    end
    
end
supertitle(fig_title)
set(gca, 'FontSize', 12, 'FontWeight', 'bold')