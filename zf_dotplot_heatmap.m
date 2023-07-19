%% zf_dotplot: gene epxressioon 
clear all
close all
clc
cond=1;
%%
% colormapx=turbo(21+28); % creat all ct  colormap  : g:7:6
if cond==1
    cmap= spring(21); %slc17a6
    str1="zf_Glut_";
    str2="Glut-";
    cd('//data/Technion_analysis/zebrafish/sc_100410/SLC17A6A/after');
load('Sorted.mat','cellid_sorted','s_xy_sorted','T_cells_tmp','data_orig_all_sorted');
load('Orignal.mat', 'geneid_all')
elseif cond==2
    cmap= winter(28); %gad2
    str1="zf_GABA";
    str2="GABA-";
    cd('/data/Technion_analysis/zebrafish/sc_100410/SLC32A1/after');
load('Sorted.mat','cellid_sorted','s_xy_sorted','T_cells_tmp','data_orig_all_sorted');
load('Orignal.mat', 'geneid_all')

else % all
    cmap= [winter(28);spring(21)]; % all filp
    str1="zf_";
    str2="";
    cd('/data/Technion_analysis/zebrafish/sc_100410/comparative/COMBINED');
    load('rerun.mat','data_orig_all','mapped_xy','cellid','geneid_all');
    s_xy_sorted=mapped_xy;
    data_orig_all_sorted=data_orig_all;
    cellid_sorted=cellid;
    T_cells_tmp=str2double(string(extractBefore(cellid,'-')));
T_cells_tmp(5163:end)=T_cells_tmp(5163:end)+28;% UPDATE!!!!!!!!!!!
end




%% update names
cellid_sortedx=extractBefore(cellid_sorted,'_10X');
u_nm=regexprep(cellid_sortedx,str1,'');
u_cx=unique(u_nm,'stable');
ut=unique(T_cells_tmp,'stable');
if 0% no name
    u_cx=strcat(str2,string(ut));
elseif cond~=3
    u_cnum=strcat(str2,string(ut),'-');
    u_cx=strcat(u_cnum,u_cx)';
end
u_cx=regexprep(u_cx,'_','-')';
% ut=unique(T_cells_tmp);
% cmap=spring(length(u_cx));
mapped_xy=s_xy_sorted;
n_data=data_orig_all_sorted;
%% caculate mean and 95th percinitle and mean of each sample mper gene per cell type
mean_data1=zeros(length(geneid_all),length(ut));
p1_95=zeros(length(geneid_all),length(ut));
for cc=1:length(ut)
    cc
    gxy = n_data(:,T_cells_tmp==ut(cc));
    %gm1=find(bx_f(gxy,1)==1);
    mean_data1(:,cc)=mean(gxy,2);
    y_percentage(:,cc)=(sum(gxy>0,2)*100/size(gxy,2))+0.001;
    p1_95(:,cc)=prctile(gxy,95,2);
    c_size(cc)=size(gxy,2);
end
norm_c_size=c_size/max(c_size);
%% for autognenes
% gname=extractBefore(cluster_name,'-');
% nangen=ismissing(gname);
% gname(nangen)=cluster_name(nangen);
% gname=unique(gname,'stable');
%% for marker specific 
% NPT: c/uc
nupt=table2array(readtable('/data/Technion_analysis/zebrafish/sc_100410/np_tf_list.xlsx'));
gname=unique(upper(string([nupt(:,1);nupt(:,2);(nupt(:,3))])),'stable');% change to add genes 
gname(cellfun('isempty',gname)) = []; % remove empty strings
gname=[flip(gname)];
% TF: c/uc 
genet=readtable('/data/Technion_analysis/zebrafish/sc_100410/TFs.xlsx');
gname=upper(string(table2array(genet(:,1))));
gname(cellfun('isempty',gname)) = []; % remove empty strings
gname=[gname;"FOXP4";"LHX6"];
%% Dotplot of markers by exp
% marker = {'SLC17A6A','SLC17A7A','GAD2','SLC32A1','Size'};%,'GSTP1','PTGDS','AQP1'};%{'STMN2','SNAP25','CLDN5','FLT1','EPAS1','ID3','GFAP','HES5','SLC1A3','ID4','SLC4A4','MFGE8','SLC1A2','SLC6A11','MPZ','PLP1','cldn19','NINJ2','SLC6A9','SOX9','OLIG2','OLIG1','SOX10','TRAF4','TAGLN','S100B'};
marker = [gname];%,'GSTP1','PTGDS','AQP1'};%{'STMN2','SNAP25','CLDN5','FLT1','EPAS1','ID3','GFAP','HES5','SLC1A3','ID4','SLC4A4','MFGE8','SLC1A2','SLC6A11','MPZ','PLP1','cldn19','NINJ2','SLC6A9','SOX9','OLIG2','OLIG1','SOX10','TRAF4','TAGLN','S100B'};
% marker=[{'SLC17A6A'};{'SLC17A7A'};{'GAD2'};{'SLC32A1'};marker;{'Size'}];
marker=[marker;{'Number of cells'}];
marker=flip(marker);

gene_data=zeros(length(marker),size(mean_data1,2)); % c/uc

%%
figure('Position', [0 0 1000 400]);
set(gcf,'color','w');
markerx=[marker;"ADCYAP1B";"SLC17A6A";"SLC32A1";"GAD2"];
for i=1:length(markerx)
    if i>1
        in=find(strcmpi(markerx(i),geneid_all));
        if ~isempty(in)
            
            mean_data_gene=mean_data1(in,:);
            gene_data(i,:)=mean_data_gene;
            y_percentage_gene=y_percentage(in,:);

            tmpthhigh= max(mean_data_gene);
            color=mean_data_gene./tmpthhigh;
            scatter(1:length(ut),i*ones(length(ut),1),(color+0.001)*100/2,cmap,'filled');hold on;
        end
    end
    if i==1
        scatter(1:length(ut),i*ones(length(ut),1),(norm_c_size)*100/2,'k','filled');hold on;
    end

end


yticks(1:length(markerx))
yticklabels(lower(markerx))
xticks(1:length(ut))
% cluster_name=regexprep(cluster_name,'_','-');
xticklabels(ut) % uc/c
% xticklabels(lower(u_cx)) % uc/c

% xticklabels(1:length(cluster_name))
gene_datax=gene_data(1:length(marker),1:size(mean_data1,2));
[~,mxi]=max(gene_datax');
[~,is]=sort(mxi);
marker=marker(is);
%
ax=gca;
ax.FontSize = 6;
set(gca, 'XAxisLocation', 'top')
xtickangle(90);
ylim([1 length(markerx)+0.5])
grid on 
axis tight;
%save2png('fig1_Dotplot_expression',gcf,700)
% export_fig fig1_Dotplot_expression.pdf
%% heatmap mean   
figure;
set(gcf,'color','w');
inall=[];
%marker = gname;%,'GSTP1','PTGDS','AQP1'};%{'STMN2','SNAP25','CLDN5','FLT1','EPAS1','ID3','GFAP','HES5','SLC1A3','ID4','SLC4A4','MFGE8','SLC1A2','SLC6A11','MPZ','PLP1','cldn19','NINJ2','SLC6A9','SOX9','OLIG2','OLIG1','SOX10','TRAF4','TAGLN','S100B'};
% marker=["SLC17A6A";"SLC17A7A";"GAD2";"SLC32A1";gname]; 
marker=[gname;"LHX9"];
for i=1:length(marker)
in=find(strcmpi(marker(i),geneid));
inall=[inall;in];
end
alltmp=1:size(data,2);

datamarkers=data(inall,:);
    Zpcay = linkage(datamarkers,'ward','correlation');
    Dx = pdist(datamarkers,'correlation');
    leafordery = optimalleaforder(Zpcay,Dx);
   datamarkers= datamarkers(leafordery,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
gr_tmp_mark=geneid(inall);
% xline(alltmp(diff(ic)))
hold on
idtmp=alltmp(logical(diff(ic)));
idtmp=[1,idtmp];
% idtmp(end)=[];
xline(idtmp)
set(gca,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark(leafordery), 'fontsize', 10)
set(gca,'xtick',idtmp,'xticklabel',uc, 'fontsize', 10)
colormap('summer');
