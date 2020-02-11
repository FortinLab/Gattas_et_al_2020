clear all;close all;clc

anim = 6
task = 'welltrained' 
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim) '\ripples'])
load('ripple_times.mat')
%%
if anim == 1
    %novel 1
    if strcmp('novel1',task)
        chan_name = 'SuperChris-Novel1-IntBadCutLFP_T';
        
        %novel 2
    elseif strcmp('novel2',task)
        chan_name = 'SuperChris-Novel2-IntBadCut.plx_T';
        
        %well trained
    else
        
        chan_name = 'SuperChris-2-12-09_SG_final_T';
        
    end
    chan_length = [1:10 12:16 18:23];
    
elseif anim ==2
    %novel1
    if strcmp('novel1',task)
        chan_name = 'Stella-Novel1-IntBadCutLFP_T';
        chan_length = [2:10 12:23];
        
        %novel2
    elseif strcmp('novel2',task)
        chan_name = 'Stella-Novel2-IntBadCutLFP_T';
        chan_length = [2:10 12:23];
        
        
        %well trained
    else        chan_name = 'Stella-2-12-2009_SG_final_T';
        chan_length = [2:10 12:16 18:23];
    end
    
elseif anim==3
    %novel 1
    if strcmp('novel1',task)
        chan_name = 'Barat-Novel1-IntBadCutLFP_T';
        %novel 2
    elseif strcmp('novel2',task)
        chan_name = 'Barat-Novel2-IntBadCutLFP_T';
        
        %well trained
    else
        chan_name = 'Barat-11-06-2008Skips_mrg_SG_final_T';
    end
    
    chan_length  = [1 3:10 12:21 23];
elseif anim==4
    if strcmp('novel1',task)
        chan_name = 'Buchanan-Novel1-IntBadCutLFP_T';
    elseif strcmp('novel2',task)
        chan_name = 'Buchanan-Novel2-IntBadCutLFP_T';
    else
        chan_name ='Buchanan4-20-withskips_mrg_SG_final_T';
    end
    chan_length = [1 2 4:10 12:13 15:23];
elseif anim ==5
    chan_name ='Mitt_July18_5odorswithSkips_SG_final_T';
    chan_length = [1:10 12:23];
    
    
elseif anim ==6
    chan_name = 'SAS01_SessiongGE44_mrg_GEcut1_T';
    chan_length = [1:2 4:10 12:18 21:24];
end
chan_counter = length(chan_length);

%% spatial map indications for each anim
[medial_elecs, lateral_elecs, i, i2,...
 subplot_row_lateral, subplot_row_medial , ...
 subplot_column_lateral, subplot_column_medial] = electrode_map( anim );

%%
% lateral elecs
figure
i_counter = 0;
for chan = lateral_elecs
    i_counter = i_counter+1;
    ripples = mean(ripple_mtx{find(chan==chan_length)},1);
    subplot(subplot_row_lateral(anim), subplot_column_lateral(anim), i2(i_counter))
    plot(ripples(700:1300))
    xlim([0 length(700:1300)])
    ylim([-.17 .2])
    title(['chan ' num2str(chan)])
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
end

% medial elecs
figure
i_counter = 0;
for chan = medial_elecs
    i_counter = i_counter+1;
    ripples = mean(ripple_mtx{find(chan==chan_length)},1);
    subplot(subplot_row_medial(anim), subplot_column_medial(anim), i(i_counter))
    plot(ripples(700:1300))
    xlim([0 length(700:1300)])
    ylim([-.17 .2])
    title(['chan ' num2str(chan)])
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
end
