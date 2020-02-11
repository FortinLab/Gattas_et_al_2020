% resample
fig_name = 'fig4'
if strcmp('fig4', fig_name)
    freq_range = freq>20 & freq<40;  % experimental cond
    start_time = 250;
    end_time   = 500;
   % mx = .4;
 %   mn = 0;
    % fig 4
    cond1_chan1 =  nanmean(squeeze(nanmean(cond1_ch1_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond1_chan2 =  nanmean(squeeze(nanmean(cond1_ch2_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond1_chan3 =  nanmean(squeeze(nanmean(cond1_ch3_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond1_chan4 =  nanmean(squeeze(nanmean(cond1_ch4_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond1 = nanmean([cond1_chan1' cond1_chan2' cond1_chan3' cond1_chan4'],2);
    
    cond2_chan1 =  nanmean(squeeze(nanmean(cond2_ch1_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond2_chan2 =  nanmean(squeeze(nanmean(cond2_ch2_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond2_chan3 =  nanmean(squeeze(nanmean(cond2_ch3_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond2_chan4 =  nanmean(squeeze(nanmean(cond2_ch4_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond2 = nanmean([cond2_chan1' cond2_chan2' cond2_chan3' cond2_chan4'],2);
    
    cond3_chan1 =  nanmean(squeeze(nanmean(cond3_ch1_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond3_chan2 =  nanmean(squeeze(nanmean(cond3_ch2_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond3_chan3 =  nanmean(squeeze(nanmean(cond3_ch3_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond3_chan4 =  nanmean(squeeze(nanmean(cond3_ch4_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond3 = nanmean([cond3_chan1' cond3_chan2' cond3_chan3' cond3_chan4'],2);
    
    cond4_chan1 =  nanmean(squeeze(nanmean(cond4_ch1_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond4_chan2 =  nanmean(squeeze(nanmean(cond4_ch2_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond4_chan3 =  nanmean(squeeze(nanmean(cond4_ch3_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond4_chan4 =  nanmean(squeeze(nanmean(cond4_ch4_temp_temp(freq_range, start_time:end_time,:),1)),1);
    cond4 = nanmean([cond4_chan1' cond4_chan2' cond4_chan3' cond4_chan4'],2);
    
elseif strcmp('fig5', fig_name)
    % fig 5
    mx =.2;
    mn = -.2;
    cond1 = zeros(size(cond1_ch1_temp_temp,3),1); cond2 = zeros(size(cond2_ch1_temp_temp,3),1);
    cond3 = zeros(size(cond3_ch1_temp_temp,3),1); cond4 = zeros(size(cond4_ch1_temp_temp,3),1);
    temp1 = []; temp2 =[]; temp3 =[]; temp4 =[];
    
    for iTrl = 1:size(cond1_ch1_temp_temp,3)
        temp1 = cond1_ch1_temp_temp(:,:,iTrl);
        temp2 = cond1_ch2_temp_temp(:,:,iTrl);
        temp3 = cond1_ch3_temp_temp(:,:,iTrl);
        temp4 = cond1_ch4_temp_temp(:,:,iTrl);
        
        cond1(iTrl)=  nanmean([nanmean(temp1(logical(zmap_index))) nanmean(temp2(logical(zmap_index))) ...
            nanmean(temp3(logical(zmap_index)))  nanmean(temp4(logical(zmap_index))) ]);
    end
    
    
    temp1 = []; temp2 =[]; temp3 =[]; temp4 =[];
    for iTrl = 1:size(cond2_ch1_temp_temp,3)
        temp1 = cond2_ch1_temp_temp(:,:,iTrl);
        temp2 = cond2_ch2_temp_temp(:,:,iTrl);
        temp3 = cond2_ch3_temp_temp(:,:,iTrl);
        temp4 = cond2_ch4_temp_temp(:,:,iTrl);
        
        cond2(iTrl)=  nanmean([nanmean(temp1(logical(zmap_index))) nanmean(temp2(logical(zmap_index))) ...
            nanmean(temp3(logical(zmap_index)))  nanmean(temp4(logical(zmap_index))) ]);
    end
    
    temp1 = []; temp2 =[]; temp3 =[]; temp4 =[];
    for iTrl = 1:size(cond3_ch1_temp_temp,3)
        temp1 = cond3_ch1_temp_temp(:,:,iTrl);
        temp2 = cond3_ch2_temp_temp(:,:,iTrl);
        temp3 = cond3_ch3_temp_temp(:,:,iTrl);
        temp4 = cond3_ch4_temp_temp(:,:,iTrl);
        
        cond3(iTrl)=  nanmean([nanmean(temp1(logical(zmap_index))) nanmean(temp2(logical(zmap_index))) ...
            nanmean(temp3(logical(zmap_index)))  nanmean(temp4(logical(zmap_index))) ]);
    end
    
    temp1 = []; temp2 =[]; temp3 =[]; temp4 =[];
    for iTrl = 1:size(cond4_ch1_temp_temp,3)
        temp1 = cond4_ch1_temp_temp(:,:,iTrl);
        temp2 = cond4_ch2_temp_temp(:,:,iTrl);
        temp3 = cond4_ch3_temp_temp(:,:,iTrl);
        temp4 = cond4_ch4_temp_temp(:,:,iTrl);
        
        cond4(iTrl)=  nanmean([nanmean(temp1(logical(zmap_index))) nanmean(temp2(logical(zmap_index))) ...
            nanmean(temp3(logical(zmap_index)))  nanmean(temp4(logical(zmap_index))) ]);
    end
end

cond1= cond1(~isnan(cond1));
cond2= cond2(~isnan(cond2));
cond3= cond3(~isnan(cond3));
cond4= cond4(~isnan(cond4));

trialCount = [length(cond1) length(cond2) length(cond3) length(cond4)];
low_bound = min(trialCount);

% draw fixed num of trials randomly from 3 conds, and generated a sampled
% dist
nperm  = 1000;
cond1_sampled_dis = zeros(low_bound,nperm);
cond2_sampled_dis = zeros(low_bound,nperm);
cond3_sampled_dis = zeros(low_bound,nperm);

% generate dist
for iPerm = 1:nperm
    cond1_sampled_dis(:,iPerm) = cond1(randi([1 length(cond1)],1,low_bound));
    cond2_sampled_dis(:,iPerm) = cond2(randi([1 length(cond2)],1,low_bound));
    cond3_sampled_dis(:,iPerm) = cond3(randi([1 length(cond3)],1,low_bound));
end
cond1_sampled_dis = mean(cond1_sampled_dis,2);
cond2_sampled_dis = mean(cond2_sampled_dis,2);
cond3_sampled_dis = mean(cond3_sampled_dis,2);
cond4_sampled_dis = cond4;

% plot all 4 trial types
figure

bar_vectora1 = [nanmean(cond1_sampled_dis) nanmean(cond4_sampled_dis) nanmean(cond2_sampled_dis) nanmean(cond3_sampled_dis)  ];
bar(bar_vectora1, 'k')
hold on
bar_vectora1sem = [nanstd(cond1_sampled_dis)/sqrt(length(cond1_sampled_dis)) nanstd(cond4_sampled_dis)/sqrt(length(cond4_sampled_dis)) ...
    nanstd(cond2_sampled_dis)/sqrt(length(cond2_sampled_dis)) nanstd(cond3_sampled_dis)/sqrt(length(cond3_sampled_dis))];
errorbar(1:4,bar_vectora1,bar_vectora1sem, 'k')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'outseq -' 'outseq +' 'inseq -' },'XTickLabelRotation',45)
%set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])

% run stats
[ pvals, p_fdr, p_masked] = TrialTypeStats( cond1_sampled_dis, cond2_sampled_dis, cond3_sampled_dis, cond4_sampled_dis )
cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim\fig4_beta')
hgexport(gcf, ['group_betapower_4TrialTypes_pooled' lock '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, ['group_betapower_4TrialTypes_pooled' lock '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')

%% fig 4 plot match/correct
plot_cond = 'correct'
    mx = .6;
    mn = -.1;
if strcmp('match', plot_cond)
    cond1_plt = nanmean([cond1_sampled_dis cond4_sampled_dis],2);
    cond2_plt = nanmean([cond2_sampled_dis cond3_sampled_dis],2);
    legend_entry = {'match' 'mismatch'};
else strcmp('correct', plot_cond)
    cond1_plt = nanmean([cond1_sampled_dis cond2_sampled_dis],2);
    cond2_plt = nanmean([cond3_sampled_dis cond4_sampled_dis],2);
    legend_entry = {'correct' 'incorrect'};
end

% traces for indiv animals
figure
bar([nanmean(cond1_plt) nanmean(cond2_plt)], 'k')
hold on
errorbar(1:2, [nanmean(cond1_plt) nanmean(cond2_plt)],...
    [nanstd(cond1_plt)/sqrt(length(cond1_plt)) nanstd(cond2_plt)/sqrt(length(cond2_plt))], 'kx', 'LineWidth',2)%ylim([mn mx])
set(gca, 'XTick', 1:2, 'XTickLabel', legend_entry,'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])


% run stats per anim
reps  = 1000;
adata = cond1_plt(~isnan(cond1_plt));
bdata = cond2_plt(~isnan(cond2_plt));
p1= permutation_unpaired(adata, bdata, reps)


cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim\fig4_beta')
hgexport(gcf, ['group_betapower_' plot_cond '_pooled_' lock '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, ['group_betapower_' plot_cond '_pooled_' lock '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
