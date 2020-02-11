function [ pvals,p_masked,p_fdr ] = ElecsAlongAxisStats( reps,cond_1_mean )

  adata = cond_1_mean(1,:); % elec 1 vs 2. 
  bdata = cond_1_mean(2,:);
  p1     = permutation_unpaired(adata, bdata, reps)
  
  adata = cond_1_mean(1,:);
  bdata = cond_1_mean(3,:);
  p2    = permutation_unpaired(adata, bdata, reps)
  
  adata = cond_1_mean(1,:);
  bdata = cond_1_mean(4,:);
  p3    = permutation_unpaired(adata, bdata, reps)
  
  adata = cond_1_mean(2,:);
  bdata = cond_1_mean(3,:);
  p4    = permutation_unpaired(adata, bdata, reps)
  
  adata = cond_1_mean(2,:);
  bdata = cond_1_mean(4,:);
  p5    = permutation_unpaired(adata, bdata, reps)
  
  adata = cond_1_mean(3,:);
  bdata = cond_1_mean(4,:);
  p6    = permutation_unpaired(adata, bdata, reps)
    
  alpha = .05
  pvals = [p1 p2 p3 p4 p5 p6];
  [p_fdr, p_masked] = fdr( pvals, alpha)

end

