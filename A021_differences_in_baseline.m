clear all;close all;clc
anim = 6
task = 'welltrained'
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load('baseline_info_wavelet_32num_20logdb_3hz_250hz_notched_artifact_reject.mat')
[ anterolat_postermedial] = ca1_electodes( anim, chan_length )
direc_chan = anterolat_postermedial;
fs = 1000
dt = 1/fs;
NumVoices = 32;
a0 = 2^(1/NumVoices);
wavCenterFreq = 6/(2*pi);
minfreq = 3;
maxfreq = 250;
minscale = wavCenterFreq/(maxfreq*dt);
maxscale = wavCenterFreq/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;
freq = wavCenterFreq./(fs*scales.*dt);
freq_range= freq>19 & freq<36; 
%anim_mean = cell(1,6);


anim_mean{anim}(1:4) = [ mean(chan_powr_mn(freq_range, anterolat_postermedial),1) ]

chan_1 = mean([anim_mean{1}(1) anim_mean{2}(1) anim_mean{3}(1) anim_mean{4}(1) anim_mean{5}(1)])
chan_2 = mean([anim_mean{1}(2) anim_mean{2}(2) anim_mean{3}(2) anim_mean{4}(2) anim_mean{5}(2)])
chan_3 = mean([anim_mean{1}(3) anim_mean{2}(3) anim_mean{3}(3) anim_mean{4}(3) anim_mean{5}(3)])
chan_4 = mean([anim_mean{1}(4) anim_mean{2}(4) anim_mean{3}(4) anim_mean{4}(4) anim_mean{5}(4)])

chan_1s = std([anim_mean{1}(1) anim_mean{2}(1) anim_mean{3}(1) anim_mean{4}(1) anim_mean{5}(1)])
chan_2s = std([anim_mean{1}(2) anim_mean{2}(2) anim_mean{3}(2) anim_mean{4}(2) anim_mean{5}(2)])
chan_3s = std([anim_mean{1}(3) anim_mean{2}(3) anim_mean{3}(3) anim_mean{4}(3) anim_mean{5}(3)])
chan_4s = std([anim_mean{1}(4) anim_mean{2}(4) anim_mean{3}(4) anim_mean{4}(4) anim_mean{5}(4)])

figure
bar ([chan_1 chan_2 chan_3 chan_4])
hold on
errorbar(1:4,[chan_1 chan_2 chan_3 chan_4],[chan_1s chan_2s chan_3s chan_4s]/sqrt(5), 'rx')
