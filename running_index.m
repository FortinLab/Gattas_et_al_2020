function [ middle_maze_idx_final ] = running_index( anim ,behavMatrix)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
info = [];
middle_maze_idx = [];
for pos = 600:700
    % find all time indices for position
    if anim ==6
        middle_maze_idx =[middle_maze_idx; find(behavMatrix(:, 19)==pos)];
    else
        middle_maze_idx =[middle_maze_idx; find(behavMatrix(:, 17)==pos)];
    end
    
    if isempty(middle_maze_idx)
        continue
    end
    
end

% find logical vector for middle of the maze indices that are 5 seconds
% apart
unq_pos = zeros(1,length(middle_maze_idx));
for a = 1:length(middle_maze_idx)
   if a ==1
       unq_pos(a)= 1;
   elseif middle_maze_idx(a)- middle_maze_idx(a-1)>5000
       unq_pos(a)= 1;
   end
end

run_trials = sum(unq_pos);
middle_maze_idx_final = middle_maze_idx(logical(unq_pos));

end

%old code
% info = [];
% for pos = 500:800
%     % find all time indices for position
%     if anim ==6
%         middle_maze_idx = find(behavMatrix(:, 19)==pos);
%     else
%         middle_maze_idx = find(behavMatrix(:, 17)==pos);
%     end
%     
%     if isempty(middle_maze_idx)
%         continue
%     end
%     unq_pos = zeros(1,length(middle_maze_idx));
% for a = 1:length(middle_maze_idx)
%    if a ==1
%        unq_pos(a)= 1;
%    elseif middle_maze_idx(a)- middle_maze_idx(a-1)>5000
%        unq_pos(a)= 1;
%    end
% end
% 
%     run_trials = sum(unq_pos);
%     info = [info; pos run_trials];
%     clear middle_maze_idx
% end
% 
% % find x-position value w/ max trial num
% row_idx  = find(info(:,2)==max(info(:,2)));
% pos_final = info(row_idx(1),1);
% 
% if anim ==6
%     middle_maze_idx = find(behavMatrix(:, 19)==pos_final);
% else
%     middle_maze_idx = find(behavMatrix(:, 17)==pos_final);
% end
% unq_pos = zeros(1,length(middle_maze_idx))
% for a = 1:length(middle_maze_idx)
%    if a ==1
%        unq_pos(a)= 1;
%    elseif middle_maze_idx(a)- middle_maze_idx(a-1)>5000
%        unq_pos(a)= 1;
%    end
% end
% middle_maze_idx_final = middle_maze_idx(logical(unq_pos));

