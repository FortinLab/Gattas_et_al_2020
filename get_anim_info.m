function [chan_name, chan_length ] = get_anim_info( anim, task )
%This function gives output of chan_name and length for each animal and task

if anim ==1
    if strcmp('novel1',task)
        chan_name = 'SuperChris-2-17-09NewOdorsSkips2_T';
    elseif strcmp('novel2',task)
        chan_name = 'SuperChris-Novel2-IntBadCut_T';
    else
        chan_name = 'SuperChris-2-12-09_SG_final_T';
    end
        chan_length = [1:10 12:16 18:23];

elseif anim ==2
    if strcmp('novel1',task)
        chan_name = 'Stella-2-17-09-newodors2_T';
            chan_length = [2:10 12:23];

    elseif strcmp('novel2',task)
        chan_name = 'Stella-Novel2-IntBadCut_T'; % chan 13 file is corrupt
            chan_length = [2:10 12 14:23];

    else
        chan_name = 'Stella-2-12-2009_SG_final_T';
            chan_length = [2:10 12:23];
    end
    
elseif anim==3
    if strcmp('novel1',task)
        chan_name = 'Barat-11-08-2008_NovelOdors_mrg_T';
        chan_length = [1:10 12:21 23]; % didnt use chan 22
    elseif strcmp('novel2',task)
        chan_name = 'Barat-Novel2-IntBadCut_T'; 
        chan_length = [1:10 12:21 23]; % didnt use chan 22
    else
        chan_name = 'Barat-11-06-2008Skips_mrg_SG_final_T';
        chan_length = [1 3:10 12:21 23];
    end
    
elseif anim==4
    if strcmp('novel1',task)
        chan_name ='Buchanan4-21-NewOdors2ndtry_T';
    elseif strcmp('novel2',task)
        chan_name = 'Buchanan-Novel2-IntBadCut_T';
    else
        chan_name ='Buchanan4-20-withskips_mrg_SG_final_T';
    end
    chan_length = [1 2 4:10 12:13 15:23];
    
elseif anim==5
    chan_name ='Mitt_July18_5odorswithSkips_SG_final_T';
    chan_length = [1:10 12:23];
elseif anim ==6
    if strcmp('novel1',task)
        chan_name ='';
    elseif strcmp('novel2',task)
        chan_name = '';
    else
        chan_name = 'SAS01_SessiongGE44_mrg_GEcut1_T';
    end
        chan_length = [1:2 4:10 12:18 21:24];
end
end

