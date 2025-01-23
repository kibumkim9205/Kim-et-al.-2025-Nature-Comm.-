%% Initializes and clears the workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% Setting up row and column numbers for each condition %%%%%%%%%%%%%%%%%%%
option=2; %1:MCF-10A (E2F), 2: MCF-10A (Cdc25A), MCF-10A (ACTB);
siteV=1:32; 
Row_1=1; Row_2=2; Row_3=3;
if option==1 %MCF-10A (E2F) 
    conditions={
        'tit1E2F',Row_1,1,siteV; %
        'tit2E2F',Row_1,2,siteV; %
        'tit3E2F',Row_1,3,siteV; %
        'tit4E2F',Row_1,4,siteV; %    
        'tit5E2F',Row_1,5,siteV; %  
        };
elseif option==2 %MCF-10A (cdc25A) *Order: Col1-Col4-Col2-Col3-Col5!!!
    conditions={
        'tit1cdc25a',Row_2,1,siteV; %
        'tit2cdc25a',Row_2,2,siteV; %
        'tit3cdc25a',Row_2,3,siteV; %
        'tit4cdc25a',Row_2,4,siteV; %    
        'tit5cdc25a',Row_2,5,siteV; %  
        };
elseif option==3 %MCF-10A (ACTB)
    conditions={
        'tit1ACTB',Row_3,1,siteV; %
        'tit2ACTB',Row_3,2,siteV; %
        'tit3ACTB',Row_3,3,siteV; %
        'tit4ACTB',Row_3,4,siteV; %    
        'tit5ACTB',Row_3,5,siteV; %  
        };
end
%% GENERAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='C:\';
experimentpath='\';
resultdir=[projectpath,experimentpath,'\'];
datadir=[projectpath,experimentpath,'\'];

if ~exist(resultdir,'dir')
    mkdir(resultdir)
end
%% Cell cycle gating
close all;clc;
numCell=3000;
[Data label]=CellCycleGating(conditions,datadir,option,numCell);
cellPhase={'G1','S','G2'};
figu(0.5,1);
%% CHOOSE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% IFdata=[xcoor,ycoor,nuc_area,sigA1,sigA2,sigA3arealow,sigA3sumlow,sigA3countlow,sigA3areahigh,sigA3sumhigh,sigA3counthigh];

label=conditions(:,1);
cellPhase={'G1','S','G2'};
saveoption=0;

option2=2; % 1: with normalization, 2: without normalization

protein=9; protein2=11;
cellS=3; % 1: G1; 2: S; 3: G2; 4: all
DisplayOption='Panel'; %Scatter, plot, Panel, boxplot

counter=0;
for i=cellS
    counter=counter+1;
    data{counter}=Stain_HistComparison_HW_GOI(conditions,datadir,resultdir,protein,protein2,DisplayOption,i,option,option2);
end

rawdata = data{1,1};
condition = size(rawdata, 2); %get number columns
normdata = cell(1, condition); %make new matrix of 9 columns except control (condition - 1)
meanControl = mean(rawdata{1, 1}); %get the mean of column 1, control
for i = 1:condition %ivide the current column by the first column
normdata{1, i} = mean(rawdata{1, i}) ./ meanControl;
end

if saveoption
    save([resultdir 'pFOXM1_',cellLine,'_20240206.mat'],'label','cellPhase','data','Data');
end
% figu(0.4,1)
%% Violin plot (Single cell)
close all;
clear Data;
numCond=size(data{1},2);
count=0;

for i=1:numCond
    count=count+1;
    Data{i}=data{1}{count};
    meanVal(i)=nanmean(Data{i});
    cellNum(i)=3000;
end

fig_position = [200 200 600 400];
figure('Position', fig_position);
plot_single_cells_violin_shape([1:numCond],Data,0.2,2000,[0.4 0.4 0.4],'dots',2);

violin(Data,'facecolor',[0 1 0;1 1 0;1 0 0;0 1 0;0 1 0;1 1 0;1 0 0;0 1 0;1 1 0;1 1 0;1 1 0;1 1 0],'edgecolor','b',...
   'bw',3.5,...
   'mc','r',...
   'medc','b--')

set(gca,'Xtick',[1:5],'XTickLabel',{'tit1(GM-GFS)','tit2','tit3','tit4','tit5(Full-Growth)'},'fontsize',10,'tickdir','out','linewidth',2);
ylabel('count (Cdc25a)','FontSize',12); ylim([0 10]);

figu(0.4,0.3);
