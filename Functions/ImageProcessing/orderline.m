function order=orderline(coorset,epcoorset)
numpoints=size(coorset,1);
order=zeros(numpoints,1);
relcoor=zeros(size(coorset));
%%% find endpoint index and set it to first in the order %%%%%%%%%%%%%%%%%%
epidx=find(coorset(:,1)==epcoorset(1) & coorset(:,2)==epcoorset(2));
order(1)=epidx;
%%% set coordinates relative to endpoint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relcoor(:,1)=coorset(:,1)-epcoorset(1);
relcoor(:,2)=coorset(:,2)-epcoorset(2);
%%% remove endpoint from consideration in remaining relative coordinates %%
relcoor(epidx,:)=[10000,10000];
%%% find coordinate adjacent to endpoint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adjidx=find(abs(relcoor(:,1))<=1 & abs(relcoor(:,2))<=1);
%%% repeat to complete the ordering of contiguous coordinates %%%%%%%%%%%%%
for pp=2:numpoints-1
    order(pp)=adjidx;
    x=relcoor(adjidx,1); y=relcoor(adjidx,2);
    relcoor(:,1)=relcoor(:,1)-x;
    relcoor(:,2)=relcoor(:,2)-y;
    relcoor(adjidx,:)= [10000,10000];
    adjidx=find(abs(relcoor(:,1))<=1 & abs(relcoor(:,2))<=1);
end
order(numpoints)=adjidx;
end