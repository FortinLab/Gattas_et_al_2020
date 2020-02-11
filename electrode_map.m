function [ medial_elecs, lateral_elecs, i, i2,...
    subplot_row_lateral, subplot_row_medial , ...
    subplot_column_lateral, subplot_column_medial] = electrode_map( anim )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if anim ==1
    medial_elecs = [15 16 12 14 18 13 19 23 22 21 20];
    lateral_elecs = [1 2 4 3 5 6 9 7 10 8];
    i  = [2 4 6 8 10 12 15 16 18 22 24];
    i2 = [1 3 5 8 10 17 21 24 26 29]; %lateral
elseif anim ==2
    medial_elecs = [12 5 8 6 2 10 7 3 9 4];
    lateral_elecs = [22 23  15 14 21 17 16 20 18 19]; % removed from 3rd member of lateral 13
    i  = [4  8 14 16 18 19 21 23 26 28];
    i2 = [1 3 9 10 11 14 15 17 19 23]; %lateral 7 third member of lateral

elseif anim ==3
    medial_elecs = [10 5 8 6 4 7 9 3 1 ];
    lateral_elecs = [17 15 20 14 18 21 13 16 19 23 12 ];
    i2 = [2 7 8 11 13 15 16 18 19 20 22]; %lateral
    i  = [2 6 8 12 16 17 19 21 25];
elseif anim == 4
    medial_elecs  = [9 1 10 8 6 2 4 7 5];
    lateral_elecs = [18 17 15 19 16 13 23 20 21 12 22];
    i2 = [3 8 10 11 13 15 17 18 24 26 29]; %lateral
    i  = [4 6 9 11 13 15 16 26 30];
elseif anim ==5
    medial_elecs  = [3 4 10 2 8 9 1 6 7 5];
    lateral_elecs = [17 18 19 16 23 20 14 15 21 13 22 12];
    i2 = [3 5 6 8 10 12 13 15 17 20 24 28]; %lateral
    i  = [3 5 6 7 9 11 13 15 16 19];
elseif anim ==6
    medial_elecs  = [6 8 7 9 10 12 13 17 14 16 18 15];
    lateral_elecs = [5 2 4 1 24 23 22 21];
    i2 = [3 5:8  10:11 13];
    i  = [1:2  5:14];
end
subplot_row_lateral    = [5 5 5 5 5 5];
subplot_row_medial     = [5 5 5 4 4 5];

subplot_column_lateral = [6 5 5 6 6 3];
subplot_column_medial  = [5 6 5 8 5 3];
end

