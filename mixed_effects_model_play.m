load flu
flu2 = stack(flu,2:10,'NewDataVarName','FluRate',...
    'IndVarName','Region');
flu2.Date = nominal(flu2.Date);
flu2 = dataset2table(flu2);
plot(flu2.WtdILI,flu2.FluRate,'ro')
xlabel('WtdILI')
ylabel('Flu Rate')
lme = fitlme(flu2,'FluRate ~ 1 + WtdILI + (1|Date)')

altlme = fitlme(flu2,'FluRate ~ 1 + WtdILI + (1|Date) + (WtdILI-1|Date)',...
'Exclude',[98,107])

lme = fitlme(data,'betapower ~ 1 + inseqcorr + inseqcorr + inseqcorr + inseqcorr + (1|anim)'+ inseqcorr|anim + ...+outseq|anim)
