%% gene correlation analysis 1
close all
clear all
clc
%%
direct='/data/Technion_analysis/zebrafish/sc_100410/SLC32A1/after';% SLC17A6A % SLC32A1
cd(direct)
% nupt=table2array(readtable('/data/Technion_analysis/zebrafish/sc_100410/list of neuropeptides.xlsx'));
% currnupt=upper(string([nupt(:,1);(nupt(:,3))]));
% currnupt(cellfun('isempty',currnupt)) = []; % remove empty strings
%% load
disp('zf_v')
load([direct,'/Sorted.mat'],'data_orig_all_sorted','sample_sorted','bar_ar_sorted','T_cells_tmp','s_xy_sorted') % vdata
load([direct,'/Orignal.mat'], 'geneid_all')%vgenes
v_gen1=geneid_all;
%% find specific gene subset 
genef="HCRT";
% gf1=strcmpi(geneid_all,genef(1))
subdata1=data_orig_all_sorted(:,:);
%% gene correlation analysis 1
direct='/data/Technion_analysis/zebrafish/sc_100410/SLC17A6A/after';% SLC17A6A % SLC32A1
cd(direct)
% nupt=table2array(readtable('/data/Technion_analysis/zebrafish/sc_100410/list of neuropeptides.xlsx'));
% currnupt=upper(string([nupt(:,1);(nupt(:,3))]));
% currnupt(cellfun('isempty',currnupt)) = []; % remove empty strings
%% load
disp('zf_v')
load([direct,'/Sorted.mat'],'data_orig_all_sorted','sample_sorted','bar_ar_sorted','T_cells_tmp','s_xy_sorted') % vdata
load([direct,'/Orignal.mat'], 'geneid_all')%vgenes
v_gen2=geneid_all;
%% find specific gene subset 
% genef="HCRT";
% gf2=strcmpi(geneid_all,genef(1))
subdata2=data_orig_all_sorted(:,:);
%%
[inBoth,xi,xf] = intersect(v_gen1, v_gen2);
m_data=subdata1(xi,:);
mr_data=subdata2(xf,:);
all_data=[m_data,mr_data];
ex_genes=upper({'MT-','zgc:','si:'});
excluded=startsWith(inBoth,ex_genes); % to be excluded genes 
all_data(excluded,:)=[];
inBoth(excluded)=[];
in = find(sum(all_data>0,2)>50);% a gene that is expressed in more than 50 cells
all_data=all_data(in,:);
z_genes=inBoth(in); 
gf=strcmpi(z_genes,genef(1));
fdata=all_data(gf,:);
% exclude those from corealtion:  

%% corr:
 [RHO,~] = corr(fdata',all_data','Type','Spearman');
 %% tops
%  [topuv,topui]=sort(RHO,'ascend');
%  topup=topup(1:10);
    
  [topdv,topdi]=sort(RHO,'descend');
  topten=topdi(1:10);
  z_top=z_genes(topten);
  r_top=topdv(1:10);
  save('/data/Technion_analysis/zebrafish/sc_100410/HCRT_analysis/top10c.mat','z_top','r_top')

%% heatmap ofexpression sorted by RHO  
figure;
set(gcf,'color','w');
inall=[];
%marker = gname;%,'GSTP1','PTGDS','AQP1'};%{'STMN2','SNAP25','CLDN5','FLT1','EPAS1','ID3','GFAP','HES5','SLC1A3','ID4','SLC4A4','MFGE8','SLC1A2','SLC6A11','MPZ','PLP1','cldn19','NINJ2','SLC6A9','SOX9','OLIG2','OLIG1','SOX10','TRAF4','TAGLN','S100B'};
% marker=["SLC17A6A";"SLC17A7A";"GAD2";"SLC32A1";gname]; 
marker=z_top;
for i=1:length(marker)
in=find(strcmpi(marker(i),z_genes));
inall=[inall;in];
end
alltmp=1:size(all_data,2);

datamarkers=all_data(inall,:);
    Zpcay = linkage(datamarkers,'ward','correlation');
    Dx = pdist(datamarkers,'correlation');
%     leafordery = optimalleaforder(Zpcay,Dx);
%    datamarkers= datamarkers(leafordery,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
gr_tmp_mark=z_genes(inall);
% xline(alltmp(diff(ic)))
% hold on
% idtmp=alltmp(logical(diff(ic)));
% idtmp=[1,idtmp];
% % idtmp(end)=[];
% xline(idtmp)
set(gca,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)
% set(gca,'xtick',idtmp,'xticklabel',uc, 'fontsize', 10)
colormap('summer');
%% colocalization:" calculate ratio of coexpression with HCRT 
expth=5;%threshold
cogen=["HCRT","CRHB","NPVF","QRFP","ANXA13L","HMX3A","LHX9"];
coscore=zeros(length(cogen),1);
for i=1:length(cogen)
gf=find(z_genes==cogen(i));
thdata=all_data(gf,:)>expth;

if i==1

basex=thdata;
end
coscore(i)=sum(thdata & basex)/sum(basex);
end
% thdata=codata>expth;
% coscore=sum(thdata & basex)/basex)';
% draw bar 
X = categorical(cogen);
X = reordercats(X,cogen);
Y = coscore;

figure('color','w')
bx=bar(X,Y);
labels1 = string(bx(1).YData);
xtips1 = bx(1).XEndPoints;
ytips1 = bx(1).YEndPoints;
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylabel(['% Colocalization with HCRT , EXP-TH:',num2str(expth)])

