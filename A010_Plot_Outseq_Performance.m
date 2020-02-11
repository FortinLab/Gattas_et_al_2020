inSeqLog     = [pokeInAlignedBehavMatrix.TranspositionDistance]==0;
otSeqLog     = [pokeInAlignedBehavMatrix.TranspositionDistance]~=0;
odor_excld_A = [pokeInAlignedBehavMatrix.Odor]~=1;
corrTrlLog   = [pokeInAlignedBehavMatrix.Performance]==1;

inSeqLog_exldA = inSeqLog==1 & odor_excld_A==1
otSeqLog_exldA = otSeqLog==1 & odor_excld_A==1

inSeqLog_exldA_corr = inSeqLog_exldA==1 & corrTrlLog==1
otSeqLog_exldA_corr = otSeqLog_exldA==1 & corrTrlLog==1

inseq_perf = sum(inSeqLog_exldA_corr)/sum(inSeqLog_exldA)
outseq_perf = sum(otSeqLog_exldA_corr)/sum(otSeqLog_exldA)

figure;
bar([inseq_perf outseq_perf])
ylim([0 1])
set(gca, 'XTickLabel',{'inseq ','outseq '}, 'XTickLabelRotation',45,  'FontSize',14,'FontWeight','Bold')
ylabel('percet correct')
title(['animal ' num2str(anim) ': ' task ])
%% all anims

anim_perf(anim,:) = [inseq_perf outseq_perf]

figure;
bar([anim_perf(:,2)])
title('animal performance on outseq trials')
ylabel('% outseq corr')
xlabel('animal #')
set(gca, 'FontSize',14,'FontWeight','Bold')
