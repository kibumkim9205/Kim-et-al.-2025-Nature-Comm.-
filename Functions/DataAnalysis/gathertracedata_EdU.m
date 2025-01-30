function [tracedata,tracestats,motherstats,IFdata,EdUdata]=gathertracedata_EdU(datadir,shot,motheroption,daughteroption,IFoption,EdUoption)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy');
if IFoption
    load([datadir,'IF_',shot,'.mat'],'IFdata');
    load([datadir,'EdU_',shot,'.mat'],'EdUdata');
else
    IFdata=ones(size(genealogy))*NaN;
    EdUdata=ones(size(genealogy))*NaN;
end
%%% gate by genealogy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gatemother=zeros(size(genealogy));
mothers=genealogy(~isnan(genealogy));
cellid=(1:length(genealogy))';
nonmothers=cellid(~ismember(cellid,mothers));
if motheroption==0 %no gating
    gatemother(:)=1;
elseif motheroption==1 %only traces that end in mitosis
    gatemother(mothers)=1;
elseif motheroption==2 %only traces that don't end in mitosis
    gatemother(nonmothers)=1;
end
if daughteroption==0
    gatedaughter=ones(size(genealogy)); %no gating
elseif daughteroption==1
    gatedaughter=~isnan(genealogy); %only traces that start with mitosis
elseif daughteroption==2
    gatedaughter=isnan(genealogy); %only traces that don't start with mitosis
end
gategenealogy=gatemother & gatedaughter;
%%% gate by presence of IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IFoption
    gateIF=~isnan(IFdata(:,1));
else
    gateIF=ones(size(genealogy));
end
if EdUoption
    gateEdU=~isnan(EdUdata(:,1));
else
    gateEdU=ones(size(genealogy));
end
%%% combine gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=gategenealogy & gateIF & gateEdU;
samplecellsID=find(samplecells);
%%% get full durations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
tracedata=linkancestry(tracedata,tracestats,samplecellsID);
%%% record mother stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if daughteroption==1
    motherstats=getmotherstatsonly(tracedata,tracestats,samplecellsID);
else
    motherstats=ones(size(genealogy))*NaN;
end
%%% remove gated out cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracedata=tracedata(samplecells,:,:);
tracestats=tracestats(samplecells,:);
motherstats=motherstats(samplecells,:);
IFdata=IFdata(samplecells,:);
EdUdata=EdUdata(samplecells,:);
end