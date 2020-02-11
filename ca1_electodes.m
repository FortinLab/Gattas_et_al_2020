function [ anterolat_postermedial ] = ca1_electodes( anim , chan_length)
%   this function takes in the input of animal number, and gives output of 4
%   channels chosen along the CA1 axis.
if anim ==1
    anterolat_postermedial= [find(chan_length==1) find(chan_length==3) find(chan_length==9) find(chan_length==7)];
elseif anim ==2
    anterolat_postermedial =  [find(chan_length==22) find(chan_length==15) find(chan_length==8) find(chan_length==3)];
elseif anim ==3
    anterolat_postermedial  =   [find(chan_length==17) find(chan_length==19) find(chan_length==7) find(chan_length==1)];
elseif anim ==4
    anterolat_postermedial =  [find(chan_length==18) find(chan_length==20) find(chan_length==8) find(chan_length==5)];
elseif anim ==5
    anterolat_postermedial =  [find(chan_length==17) find(chan_length==21) find(chan_length==1) find(chan_length==5)];
elseif anim ==6
    anterolat_postermedial =  [find(chan_length==1) find(chan_length==10) find(chan_length==14) find(chan_length==15)];
end

end

