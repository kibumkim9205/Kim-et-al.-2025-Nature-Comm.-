%% Initializes and clears the workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
%% GENERAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='\';
experimentpath='\';
datadir=[projectpath,experimentpath,'Data\'];
resultdir=[projectpath,experimentpath,'Results\'];
if ~exist(resultdir,'dir')
    mkdir(resultdir)
end
%%
rows=1:4;
cols=1:12;
sites=1:32;

plot=0; %1: Check histogram; 0: collect data;

numrows=length(rows);
numcols=length(cols);
numsites=length(sites); 
shots=numrows*numcols; % Purpose: designate numrows, cols, sites, shots

for idx=1:shots
    colidx=mod(ceil(idx),numcols); %mod: remainder afer division; ceil: round toward positive infinity
    if colidx==0
        colidx=numcols; % Purpose: 
    end
    col=cols(colidx);
    rowidx=ceil(idx/(numcols));
    row=rows(rowidx);

    Alldata=[]; % Purpose: make 'All data' vector
    for site=sites
        %Nikon
        [rowstring,colstring,sitestring]=wellnum2strRCS_3(row,col,site);
        shot=[rowstring,colstring,'_',sitestring];
        if exist([datadir,'IF_',shot,'.mat'])
            load([datadir,'IF_',shot,'.mat'],'IFdata');
            Alldata=[Alldata;IFdata];
        else
            continue
        end
    end
    if isempty(Alldata)
        continue
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);  % Purpose: designate Hoechst, EdU value in Alldata 
    EdUval=Alldata(:,5);
    EdUval(EdUval<1)=1; % log(1)=0, so convert it to 1
    EdUval=log2(EdUval); %log2(EdUval)=2^EdUval (produce concies graphs by preventing units from becoming too large

     if row==1:2
        threshEdU=9;
    elseif row==3:4
        threshEdU=8;  
     end
    
    if plot==1
        if idx==1
            figure;
        end
        subplot(2,12,idx); histogram(Hoechstval);
        subplot(2,12,idx+12); histogram(EdUval); figu(0.3,0.6);
        numCell_nogating_Data(row,col)=size(Hoechstval,1);
        vline(thr5eshEdU);
        ylabel(['Row: ' num2str(row) ' Col: ' num2str(col)])
        figu(0.5,1)
    end

    %Rb=Alldata(:,5);%./Alldata(:,protein);
    %posval=Rb>=0; Rb(Rb<1)=1; Rb=log2(Rb);

    %Cas9=Alldata(:,6);%./Alldata(:,protein);
    %posval_2=Cas9>=0; Cas9(Cas9<1)=1; Cas9=log2(Cas9);

    %gating=posval & posval_2;% & CDK
    %EdUval=EdUval(gating);
    %Rb=Rb(gating);
    %Cas9=Cas9(gating);
    %Hoechstval=Hoechstval(gating);

    Hoechst_gating=Hoechstval>500000;
    EdU_gating=EdUval>1 & EdUval<16;

    gating=Hoechst_gating & EdU_gating;
    numCell=sum(gating);
    EdUval=EdUval(gating);
    Hoechstval=Hoechstval(gating);


    numS=sum(EdUval>threshEdU);
    perS=numS/numCell*100;

    %%%%%% visualize Hoechst vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure;dscatter(Hoechstval,EdUval);
    figure;histogram(Cas9,100,'Normalization','probability'); vline(threshCas9); xlim([8 13]);
    figure;histogram(EdUval,100,'Normalization','probability'); vline(threshEdU); xlim([3 14.5]);
    %}
    numCellData(row,col)=numCell;

    HoechstRaw{row,col}=Hoechstval;
    EdURaw{row,col}=EdUval;

    SData(row,col)=perS;
    %%%%%% visualize Hoechst vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure;histogram(Hoechstval);
    figure;histogram(EdUval);
    figure;dscatter(Hoechstval,EdUval);
    figure;dscatter(CDK2,EdUval);
    figure;dscatter(Rb,EdUval);
    figure;dscatter(CDK4_stain,EdUval);
    figure;histogram(EdUval,100,'Normalization','probability'); vline(threshEdU); xlim([3 14.5]);
    %}
end
SData=fliplr(SData);
%%
conditions={'RPE1+Adagrasib','RPE1+Adagrasib','RPE1+Sotorasib','RPE1+Sotorasib'};
%%
save([resultdir 'Data_'],'HoechstRaw','EdURaw','numCellData','SData','conditions');
%%
row=5; col=12;
figure;histogram(EdURaw{row,col});
vline(10);

%%
row=3; col=12;
figure;
dscatter(HoechstRaw{row,col},EdURaw{row,col})
hline(10);

%%
yMax=20000;
figure;imagesc(numCellData(1:end,1:end),[0 yMax]); colorbar; title('Number of Cells');
set(gca,'Xtick',[1:12],'XTickLabel',[1:12],'fontsize',12,'tickdir','out','linewidth',2);
set(gca,'Ytick',[1:8],'YTickLabel',{'RPE1(Adagrasib)','RPE1(Adagrasib)','RPE1(Sotorasib)','RPE1(Sotorasib)','MCF-10A(Adagrasib)','REP1(Adagrasib)','MCF-10A(RM042)','RPE1(RM042)'},'fontsize',12,'tickdir','out','linewidth',2);
colormap('cyan_yellow');
%colormap('blue_yellow')
figu(0.5,0.6);
%%
gating=numCellData>1000;
%%
figure;imagesc(gating); colorbar; title('Gating');
set(gca,'Xtick',[1:12],'XTickLabel',[1:12],'fontsize',12,'tickdir','out','linewidth',2);
set(gca,'Ytick',[1:7],'YTickLabel',{'A549','H1373','H1792','H358','SW1573','MCF-10A','RPE1'},'fontsize',12,'tickdir','out','linewidth',2);
figu(0.5,0.6);
%%
yMax=40;
figure;imagesc(SData,[0 yMax]); colorbar; title('%S-phase intensity plate');
set(gca,'Xtick',[1:12],'XTickLabel',[0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30],'fontsize',12,'tickdir','out','linewidth',2);
set(gca,'Ytick',[1:8],'YTickLabel',{'A549(RM028)','H1792(RM028)','MCF-10A(RM028)','RPE1(RM028)','A549(RM042)','H1792(RM042)','MCF-10A(RM042)','RPE1(RM042)'},'fontsize',12,'tickdir','out','linewidth',2);
colormap('cyan_yellow');    
figu(0.5,0.6);
%colormap('blue_yellow')
%%
normOption=0;

conc=[0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000];

name={'Adagrasib, Sotorasib'};
color={'b','b','g','k','r','b','b','g'};
count=1;

figure;
for count=1:2;
    i=count;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;


    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
    ylabel('S-Phase (%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
end
figu(0.4,1)
%%
normOption=0;

conc=[0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000];

name={'Adagrasib, Sotorasib'};
color={'b','b','b','b','r','b','b','g'};
count=1;

figure;
for count=3:4;
    i=count;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;


    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
    ylabel('S-Phase (%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
end
figu(0.4,1)
%% A549+RM028
normOption=1;

conc=[0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000];

name={'RM028'};
color={'k'};
count=0;

figure;
for count=1;
    i=1;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
        SData_norm2{count}=response;
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;


    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    ylabel('S-Phase(%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
end
figu(0.4,1)
%%
combinedata='C:\1. Project\1. KRAS NSLC\Combine\IC50 RM028 NSCLC\Data\';
save([combinedata 'Data_RM28 A549-1'],'SData_norm2','SData');
save([resultdir 'Data_RM28 A549-1'],'SData_norm2','SData');
% save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');
%% H1373+RM028
normOption=1;

conc=[0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000];

name={'RM028'};
color={'k','r','b','g','k','r','b','g'};
count=0;

figure;
for count=1;
    i=2;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
        SData_norm2{count}=response;
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;


    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    ylabel('S-Phase(%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
end
figu(0.4,1)
%%
combinedata='C:\1. Project\1. KRAS NSLC\Combine\IC50 RM028 NSCLC\Data\';
save([combinedata 'Data_RM28 H1373-1'],'SData_norm2','SData');
save([resultdir 'Data_RM28 H1373-1'],'SData_norm2','SData');
% save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');
%% H1792+RM028
normOption=1;

conc=[0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000];

name={'RM028'};
color={'k','b','b','g','k','r','b','g'};
count=0;

figure;
for count=1;
    i=3;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
        SData_norm2{count}=response;
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;


    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    ylabel('S-Phase(%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
end
figu(0.4,1)
%%
combinedata='C:\1. Project\1. KRAS NSLC\Combine\IC50 RM028 NSCLC\Data\';
save([combinedata 'Data_RM28 H1792-1'],'SData_norm2','SData');
save([resultdir 'Data_RM28 H1792-1'],'SData_norm2','SData');
% save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');
%% H358+RM028
normOption=1;

conc=[0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000];

name={'RM028'};
color={'k','r','r','g','k','r','b','g'};
count=0;

figure;
for count=1;
    i=4;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
        SData_norm2{count}=response;
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;


    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    ylabel('S-Phase(%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
end
figu(0.4,1)
%%
combinedata='C:\1. Project\1. KRAS NSLC\Combine\IC50 RM028 NSCLC\Data\';
save([combinedata 'Data_RM28 H358-1'],'SData_norm2','SData');
save([resultdir 'Data_RM28 H358-1'],'SData_norm2','SData');
% save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');
%% MCF10A+RM028
normOption=1;

conc=[0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300];

name={'RM028'};
color={'k','r','b','b','r','r','b','g'};
count=0;

figure;
for count=1;
    i=5;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
        SData_norm2{count}=response;
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;
    
    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    ylabel('S-Phase(%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
end
figu(0.4,1)
%%
combinedata='C:\1. Project\1. KRAS NSLC\Combine\IC50 RM028 NSCLC\Data\';
save([combinedata 'Data_RM28 MCF10A-1'],'SData_norm2','SData');
save([resultdir 'Data_RM28 MCF10A-1'],'SData_norm2','SData');
% save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');
%% RPE1+RM028
normOption=1;

conc=[0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300];

name={'RM028'};
color={'k','r','b','g','k','b','b','g'};
count=0;

figure;
for count=1;
    i=6;
    if normOption
        SData_norm=SData(i,:);%-minS;
        response=SData_norm/SData_norm(1);
        SData_norm2{count}=response;
    else
        response=SData(i,:);
    end
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    y1=response;


    subplot(1,7,count);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    % ylim([0.0 1.2]);
    ylabel('S-Phase(%)','FontSize',12);xlabel('drug','FontSize',12);
    set(gca,'Xtick',[conc],'XTickLabel',[0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300],'fontsize',11,'tickdir','out','linewidth',1);
    title([conditions(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
    xlabel('Dose (nM)'); %ylabel('CDK2 activity');
end
figu(0.4,1)
%%
combinedata='C:\1. Project\1. KRAS NSLC\Combine\IC50 RM028 NSCLC\Data\';
save([combinedata 'Data_RM28 RPE1-1'],'SData_norm2','SData');
save([resultdir 'Data_RM28 RPE1-1'],'SData_norm2','SData');
% save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');
%% IC50 with individual dots
dose=[0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000]; % uM
count=0;
for i=1 % # conditions
    count=count+1;
    plotV=SData1(1:12,i);
    doseresponse=plotV
    %doseresponse=doseresponse/doseresponse(1);
    numMax=length(doseresponse);
    dose_gated=dose(gating(1:8,i));
    response{count}=reshape(doseresponse,[1 numMax]);
    
    %if ismember(i,1:2:30)
    %    figure; hold on;
    %end
    
    subplot(1,1,count); 
    y1=response(1,:)
    minDose=min(dose);
    maxDose=max(dose);
    xvals=(dose);
    fitpara=fit(xvals',doseresponse,'exp1'); % the two input (xvals, Response) should have the same size, add ' if needed to convert
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints, fitpara(xpoints),'Color',[0 0 0],'linewidth',2)

%     x_label_names={'Palbociclib', 'Ribociclib', 'Abemaciclib', 'Palbociclib', 'Ribociclib', 'Abemaciclib','Palbociclib', 'Ribociclib', 'Abemaciclib'}
%     xlabel(x_label_names{i});
    hold on
    halfConcentration(i)=(log(2)/abs((fitpara.b))); %convert to concentration
    xlim([dose(1) dose(end)]);
    ylim([0.0 1.2]);
    scatter(xvals',doseresponse,'Marker','o','MarkerFaceColor','k'); %plot(fitpara);
    title([conditions{i} ': IC50=' num2str(halfConcentration(i))],'fontsize',12);
        xlabel('Dose (nM)','FontSize',12);
        %ylabel('Response','FontSize',12);
        ylabel('S-phase (%)','FontSize',12);
    %if ismember(i,2:2:30)
    %    count=0;
    %    ylim([0 1.2]); %xlim([0 end]);
    %    ylabel('Response','FontSize',12);
       
        figu(0.3,0.3);
        
    %end
end
%% 
% option=1;

normOption=1;
name={'RM028'};
% numCond=size(rows,2);
color={'k','r','b','g','k','r','b','g'};
count=0;
% for i=2:2:5 %numCond
i=1;
    count=1;
      %count=count+i;
%     count=i/2;
    if normOption
        minS=min(nanmean(SData(i:i+1,:)));
        maxS=max(nanmean(SData(i:i+1,:)))-minS;
        SData_norm=nanmean(SData(i:i+1,:))-minS;
        response=SData_norm/SData_norm(1,1);
        SData_norm2{count}=response;
    else
        response=SData(count,1:11);
    end
    subplot(2,2,count);
    if count==1
        conc=[0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000];
    elseif count==2
        conc=[0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000];
        
    end
    y1=response(1,:);
    minDose=min(conc);
    maxDose=max(conc);
    xvals=(conc);
    fitpara=fit(xvals',y1','exp1');
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,fitpara(xpoints),'Color',color{i},'LineWidth',2)
    hold on;
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    scatter(xvals',y1','Marker','o','MarkerFaceColor',color{i});
    xlim([conc(1) conc(end)]);
    ylim([0.0 1.2]);
    ylabel('Response','FontSize',12);xlabel('drug','FontSize',12);
    title([name(count), 'IC50=' num2str(halfLife(i))],'fontsize',10);
    figu(0.3,0.3);
% end
xlabel('Dose (nM)'); %ylabel('CDK2 activity');

%%
combinedata='C:\1. Project\2. Breast cancer\Combine\IC50 Palbo drug-naive AT3OVA\Data\';
save([combinedata 'Data_20230427_3'],'SData_norm2','SData');
save([resultdir 'Data_20230427_3'],'SData_norm2','SData');
% save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');
%%
dose=[10 30 100 300 1000 3000 6000 10000]; % nM
count=0;


    count = count + 1; 
    response = SData (1:8,i);
    response = respons';
    figure;hold on
    
%% IC50 with individual dots

dose=[10 30 100 300 1000 3000 6000 10000]; % nM
count=0;


for i=2 % # conditions
    count=count+1;
    plotV=SData(1:8,i);
    
    doseresponse=plotV(1:8,i);
    doseresponse=doseresponse/doseresponse(1);
    numMax=length(doseresponse);
    
    response{count}=reshape(doseresponse,[1 numMax]);
    
%     if ismember(i,1:2:30)
%         figure; hold on;
%     end
    
    subplot(1,2,count); hold on;
    xvals=(dose);
    fitpara=fit(xvals',doseresponse,'exp1'); % the two input (xvals, Response) should have the same size, add ' if needed to convert
    halfConcentration(i)=(log(2)/abs((fitpara.b))); %convert to concentration
    
    scatter(xvals',doseresponse,'Marker','o','MarkerFaceColor','k'); plot(fitpara);
    title([conditions{i} ': IC50=' num2str(halfConcentration(i))],'fontsize',12);
    if ismember(i,2:2:30)
        count=0;
        ylim([0 1.2]); %xlim([0 end]);
        ylabel('Response','FontSize',12);
        xlabel('Palbociclib','FontSize',12);
        figu(0.3,0.8);
        
    end
end
%%
%% Normalization
clear hillCoeff EC50
count=0;
%dose=[0.001 0.01 0.1 0.5 1 3 5 8 10 50 100 500];
conc=[0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30]; 
% nM
 %
    count=count+1;

        responses=SData1(:,2);
        responses=responses/responses(1);
   
%     responses{count}=reshape(doseresponse,[1 8]);
    figure; hold on; 
    [hillCoeff ec50(count)]=doseResponse(conc,responses');
    ylim([0 1.2])

xlabel('Dose (nM)'); ylabel('S-phase cells (%)');
%%
% dose=[0.01 0.03 0.1 0.3 1 3 6 10];
% for i=1:3
%     response=mean(com_data{i});
%     figure; hold on;
%     [hillCoeff ec50(i)]=doseResponse(dose,response);
%     if i==1
%         ylim([0 65])
%     end
% end
% xlabel('Dose (nM)'); ylabel('S-phase cells (%)');

%% IC50

conc=[0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30];
figure; count=0;
for i=[1:12]   
  %  if ismember(i,6:11)
  %      doseresponse=mean(SData(1:8,i:i+1)')/mean(SData(1,i:i+1)'); %normalization by control condition
  %  else
        doseresponse=mean(SData1(1:12,i)')/mean(SData1(1,i)');    
  % end
    
    %MAXpoint=find(doseresponse==max(doseresponse));
    %MINpoint=find(doseresponse==min(doseresponse));
    %decreaseS=doseresponse(MAXpoint:MINpoint);        
    %minNorm=min(decreaseS);
    %decreaseS=decreaseS-minNorm;
    %decreaseS=decreaseS/max(decreaseS);
    
    count=count+1;
    subplot(1,6,count); hold on;
    xvals=(conc);
    fitpara=fit(xvals',doseresponse','exp1');
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    
    plot(fitpara,xvals',doseresponse');hold on;
    
    %     set(gca,'Xtick',[0:100:lastFrame-fisrtFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
    %     set(gca,'Ytick',[0.3:0.3:1.5],'YTickLabel',[0.3:0.3:1.5],'fontsize',24,'tickdir','in','linewidth',5)
%     ylim([0 1.2]); xlim([0 end]);
    ylabel('Response','FontSize',20);xlabel('Palbociclib','FontSize',20);title(['IC50=' num2str(halfLife(i))],'fontsize',20);
    figu(0.5,1);
end


%% Normalization
clear hillCoeff ec50
count=0;
dose=[10 30 100 300 1000 3000 10000 30000]; % nM
for i=[7:8] %
    count=count+1;
    plotV=SData(1:8,i);
    doseresponse=plotV(gating(1:8,i));
    doseresponse=doseresponse/doseresponse(1);
    numMax=length(doseresponse);
    dose_gated=dose(gating(1:8,i));
    response{count}=reshape(doseresponse,[1 numMax]);
    figure; hold on;
    
    [hillCoeff EC50V(count)]=ec50(dose_gated',response{count}');
    title(conditions{i});
    ylim([0 1.2]);
end
xlabel('Dose (nM)'); ylabel('S-phase cells (%)');
%%
conc=[10 30 100 300 1000 3000 6000 10000];
count=0;
% Response=[];
for i=[1,3,5]    
     ismember(i,1:10)
        Response=mean(SData(1:8,i:i+1)')/mean(SData(1,i:i+1)'); %normalization by control condition
     
    Response{i}=Response;
    count=count+1;    
    xvals=(conc);
    fitpara=fit(xvals',Response','exp1');
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time    
end
%%
save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');


%%
figure;imagesc(RbData(4,2:end),[8 10]); colorbar; title('Rb intensity');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Xtick',[1:10],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
%set(gca,'Ytick',[1:8],'YTickLabel',{'Bortezomib','Ixazomib','Celastrol','Carfilzomib','Epoxomicin','ONX0914','Delanzomib','Oprozomib'},'fontsize',12,'tickdir','out','linewidth',2);



%%
count=1;
for i=[2:7]
    count=count+1;
    xvals=[0.00001 0.001 5 10 50];
    val=plotV(i,2:6);
    % fitpara=fit(xvals',val','exp2');
    %     halfLife(i)=log(2)/abs((fitpara.b)); %convert to time
    figure; hold on;
    [hillCoeff ec50(count)]=doseResponse(xvals,val)
    
    
    %plot(fitpara,xvals',val');hold on;
end


