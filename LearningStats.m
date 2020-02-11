function [ pvals,p_fdr,p_masked ] = LearningStats( reps,cond1, cond2, cond3 )
 % run stats per anim
    adata = cond1;
    bdata = cond2;
    p1     = permutation_unpaired(adata, bdata, reps)
    
    adata = cond1;
    bdata = cond3;
    p2    = permutation_unpaired(adata, bdata, reps)
    
    adata = cond3;
    bdata = cond2;
    p3    = permutation_unpaired(adata, bdata, reps)
    
    alpha = .05;
    pvals = [p1 p2 p3];
    [p_fdr, p_masked] = fdr( pvals, alpha)
    disp(['novel 2 vs novel1 : ' num2str(pvals(1))])
    disp(['welltrained vs novel1 : ' num2str(pvals(2))])
    disp(['welltrained vs novel2 : ' num2str(pvals(3))])


end

