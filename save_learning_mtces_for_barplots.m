function [ ~ ] = save_learning_mtces_for_barplots( task, cond_1_all_animals)

cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')

% Save learning matrices to make bar plots
if strcmp('novel1',task)
    cond1=[];
    cond1 = cond_1_all_animals;
    save('novel1_InSeq_learning','cond1')
    
elseif strcmp('novel2',task)
    cond2=[];
    cond2 = cond_1_all_animals;
    save('novel2_InSeq_learning','cond2')
    
elseif strcmp('welltrained',task)
    cond3=[];
    cond3 = cond_1_all_animals;
    save('welltrained_InSeq_learning','cond3')
    
end

