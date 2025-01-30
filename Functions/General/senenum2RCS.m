function [row,col] = senenum2RCS(sene)
idx=sene/12;
idx_2=round(rem(idx,1)/0.0833);
if idx<=1
    row=2; if idx_2==0; col=12; else col=idx_2; end
elseif idx<=2
    row=3; if idx_2==0; col=1; else col=13-idx_2; end
elseif idx<=3
    row=4; if idx_2==0; col=12; else col=idx_2; end
elseif idx<=4
    row=5; if idx_2==0; col=1; else col=13-idx_2; end
elseif idx<=5
    row=6; if idx_2==0; col=12; else col=idx_2; end
elseif idx<=6
    row=7; if idx_2==0; col=1; else col=13-idx_2; end
end
end