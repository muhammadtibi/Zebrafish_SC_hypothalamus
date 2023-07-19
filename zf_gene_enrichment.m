%% Combined: load rerun (/data/Technion_analysis/zebrafish/zf_scratch_nocls_combined.m)
close all
clear all
clc
%%
disp('zf_sc')

direct='/data/Technion_analysis/zebrafish/sc_100410/comparative/COMBINED';% SLC17A6A % SLC32A1
cd(direct)
load([direct,'/rerun.mat'],'data_orig_all','sample','geneid_all','cellid','all_flags','mapped_xy') % vdata

set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureVisible','on');% off / on

cd('/data/Technion_analysis/zebrafish/sc_100410/HCRT_analysis')

%% get draw
genesx=geneid_all;
ball_c=extractBefore(cellid,'_10X');
[bu_c,bu_x,bu_i]=unique(ball_c);
bu_c=string(natsort(bu_c));
ball_c=string(ball_c);
xidx=1:length(bu_c);
glutc=spring(21);
gadc=winter(28);
bgad=contains(bu_c,"GABA");
bglut=contains(bu_c,"Glut");
combinec=zeros(length(bu_c),3);
combinec(bgad,:)=gadc;
combinec(bglut,:)=glutc;
all_data=data_orig_all;
ex_genes=upper({'jun','fos','atf3','egr','agr2','hsp70.1','plk3','btg','ubb','gadd45','MT-','rrad', 'npas4', 'dnajb', 'nr4a1', 'hsps1', 'dusp1', 'dusp5', 'jdp2b','zgc:','si:'});
exgi=startsWith(geneid_all,ex_genes);
all_data(exgi,:)=[]; % take out ex_genes
genesx=geneid_all(~exgi); % take out ex_genes

u_cx=regexprep(bu_c,'zf_','');
u_cx=regexprep(u_cx,'_','-');
%% HCRT flag 
% creat T_cells_tmp binary (HCRT/NO) flag!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Hcrtx=contains(u_cx,"hcrt");
T_cells_tmp= ball_c==bu_c(Hcrtx);%un/c
% or for all clusters
% T_cells_tmp=bu_i; % un/c
%% or glut.gaba ct flags
ggx=contains(ball_c,"GABA");
ball_c(~ggx)=[];
T_cells_tmp= bu_i(ggx);%un/c
all_data(:,~ggx)=[];
%% calculate mean_data 
mean_data=zeros(length(genesx),length(ball_c));

for cc=1:length(bu_c)
    cc
    gxy = all_data(:,strcmpi(bu_c(cc),ball_c));
    mean_data(:,cc)=mean(gxy,2);
end
%%  enrichment cluster gene analysis :sc
% exlude genes: 

[rs,is]=sort(T_cells_tmp,'ascend');
T_cells_tmp=T_cells_tmp(is);
all_data=all_data(:,is);
[ind_gr_tmp_mark,~,gr_center] = markertablefeatures(T_cells_tmp,all_data,3);% top_genes by me
figure;
set(gcf,'color','w');
inall=[];
%marker = gname;%,'GSTP1','PTGDS','AQP1'};%{'STMN2','SNAP25','CLDN5','FLT1','EPAS1','ID3','GFAP','HES5','SLC1A3','ID4','SLC4A4','MFGE8','SLC1A2','SLC6A11','MPZ','PLP1','cldn19','NINJ2','SLC6A9','SOX9','OLIG2','OLIG1','SOX10','TRAF4','TAGLN','S100B'};
% marker=["SLC17A6A";"SLC17A7A";"GAD2";"SLC32A1";gname];
marker=genesx(ind_gr_tmp_mark);
for i=1:length(marker)
    in=find(strcmpi(marker(i),genesx));
    inall=[inall;in];
end
alltmp=1:size(all_data,2);
idx=1:length(u_cx);
datamarkers=all_data(inall,:);
% datamarkers(:,1:28)=datamarkers(:,2:2:length(u_cx));
% datamarkers(:,29:end)=datamarkers(:,1:2:length(u_cx))
% clnames(1:21)=u_cx(2:2:length(u_cx));
% clnames(22:end)=u_cx(2:2:length(u_cx));
Zpcay = linkage(datamarkers,'ward','correlation');
Dx = pdist(datamarkers,'correlation');
% leafordery = optimalleaforder(Zpcay,Dx);
% datamarkers= datamarkers(leafordery,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
% datamarkers_cn(2:)
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
% gr_tmp_mark=geneid_all(inall);
% xline(alltmp(diff(ic)))
hold on
idtmp=alltmp(logical(diff(rs)));
% idtmp=[idtmp];
% idtmp(end)=[];
xline(idtmp,'k','LineWidth',2)
set(gca,'ytick',[1:length(marker)],'yticklabel',marker, 'fontsize', 10)
set(gca,'xtick',gr_center,'xticklabel',1:length(gr_center), 'fontsize', 10)
colormap('summer');
%% choose the topgenes 
[rs,is]=sort(T_cells_tmp,'ascend');
T_cells_tmp=T_cells_tmp(is);
all_data=all_data(:,is);
[ind_gr_tmp_mark,vi0] = markertablefeatures_tmp(T_cells_tmp,all_data,10,1);% top_genes by me
top10g=genesx(ind_gr_tmp_mark(:,2));
en_score=vi0(:,2);
save('/data/Technion_analysis/zebrafish/sc_100410/HCRT_analysis/top10g.mat','top10g','en_score')
%% hbar genes ( u can start from here)
load('/data/Technion_analysis/zebrafish/sc_100410/HCRT_analysis/top10c.mat')
% load('/data/Technion_analysis/zebrafish/sc_100410/HCRT_analysis/top10g.mat')
% [inBoth,xi,xf]=intersect(z_top,top10g);
% r_top=r_top(xi);   
% en_score=en_score(xf); 
r_top(z_top=="HCRT")=[];
z_top(z_top=="HCRT")=[];
%%
X = categorical(z_top);
X = reordercats(X,z_top);
Y = [r_top];
figure('color','w')
barh(X,Y)
xlabel('Rho score with HCRT')
%% hbar genes ( u can start from here)
% load('/data/Technion_analysis/zebrafish/sc_100410/HCRT_analysis/top10c.mat')
load('/data/Technion_analysis/zebrafish/sc_100410/HCRT_analysis/top10g.mat')
% [inBoth,xi,xf]=intersect(z_top,top10g);
% r_top=r_top(xi);   
% en_score=en_score(xf); 
en_score(top10g=="HCRT")=[];
top10g(top10g=="HCRT")=[];

X = categorical(top10g);
X = reordercats(X,top10g);

X = reordercats(X,top10g);
Y = [en_score];
figure('color','w')
barh(X,Y)
xlabel('Enrichment score with HCRT cluster')
