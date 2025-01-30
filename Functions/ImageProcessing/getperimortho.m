function cporthos=getperimortho(perimset,perimstep,cpidx)
perilength=size(perimset,1);
leftidx=mod(cpidx-perimstep,perilength); leftidx(leftidx==0)=perilength;
leftcoors=perimset(leftidx,:);
rightidx=mod(cpidx+perimstep,perilength); rightidx(rightidx==0)=perilength;
rightcoors=perimset(rightidx,:);
pointdiff=rightcoors-leftcoors;
%cptangents=rad2deg(atan2(pointdiff(:,2),pointdiff(:,1)));
cptangents=rad2deg(atan(pointdiff(:,2)./pointdiff(:,1)));
cporthos=cptangents+90; %ranges from 0-180
end