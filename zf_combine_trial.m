
%% gf_vs_m
clear all
close all
clc
nrmz=5000;
%% load mouse-new (my)
disp('mm-new')
load('/bigdata/all_mouse_dataset/loom_n.mat','mn_matrix','mn_genes','mn_region','mn_marker','mn_location','mn_cl','mn_cln','mn_id','mn_dev','mn_ts','mn_nt','mn_t3');
m_nt=mn_nt;
m_data=mn_matrix';
m_name=mn_cln;
m_genes=mn_genes;
m_sample=ones(length(mn_id),1);%mn_id;
m_region=mn_region;
m_location=mn_location;
m_marker=mn_marker;
m_rank=string(mn_t3);
%% sampling 
sp=1;
if sp==1
disp('sample mm')
n_m=unique(m_name,'stable');
y_all=[];
for rs=1:length(n_m)
    ig=find(m_name==n_m(rs));% for a cell cluster
    if  length(ig)>200
        rs
        y = randsample(ig,200); % here we take only XX cells of the real cell typoe cluster
        y_all=[y_all;y];
    else
        y_all=[y_all;ig];
    end

end
m_name=m_name(y_all);
m_region=m_region(y_all);
m_sample=m_sample(y_all);
m_rank=m_rank(y_all);
m_marker=m_marker(y_all);
m_location=m_location(y_all);
m_nt=m_nt(y_all);
m_data=m_data(:,y_all);
end

%% a6 zf : 3
disp('zf')
disp('a6')
% nrmz=20000;
cd('/data/Technion_analysis/zebrafish/sc_100410/SLC17A6A')
load('Sorted.mat', 'data_orig_all_sorted','T_cells_tmp','sample_sorted','all_flags_sorted','s_xy_sorted','cellid_sorted')
load('Orignal.mat', 'geneid_all')


% > find  all_flags_sorted
% usbat=natsort(unique(sample_sorted,'stable'));
% ur=[2 2 2 1 1 1 2 1 2 1]; % r=1 / u=2
% bh=['h','b','b','h','h','h','b','b','b','h'];
all_flags_sorted=string(zeros(length(sample_sorted),1));
% for i= 1: length(usbat)
% 
%     flags_i=contains(sample_sorted,usbat(i));
%     all_flags_sorted(flags_i)=ur(i);
% end
% < flag


a6data=data_orig_all_sorted;%round(data_orig_all_sorted./repmat(sum(data_orig_all_sorted),length(data_orig_all_sorted(:,1)),1)*nrmz);
% [~, txt, ~] = xlsread('SLC17A6_Clusters.xlsx',1,'B:B');
% a6_name=txt([2:end],:);
txt=readtable('clusters names glut.xlsx');
a6_name=table2array(txt(:,end));% length is number of clusters
% add glut2 before each cell type name
a6_name=insertBefore(a6_name,1,'zf_Glut_');
sa6=sample_sorted;
ia6=cellid_sorted;


% group simiilar clusters
[u_c,~,X] = unique(a6_name,'stable');% X is replica index
ca6=strings(length(T_cells_tmp),1);
T_cells_x=zeros(size(T_cells_tmp));
for cl=1:length(a6_name)
    ca6(T_cells_tmp==cl)=string(u_c(X(cl)));
    T_cells_x(T_cells_tmp==cl)=X(cl);
end
% remove exclude names cluster
out_idx=ca6=="zf_Glut_exclude";
sa6(out_idx)=[];
a6data(:,out_idx)=[];
T_cells_x(out_idx)=[];
ca6(out_idx)=[];
all_flags_sorted(out_idx,:)=[];
s_xy_sorted(out_idx,:)=[];
ia6(out_idx)=[];

ba6_genes=geneid_all;
fa6=all_flags_sorted;
X=unique(X);
try
    X(X==find(u_c=="zf_Glut_exclude"))=[];
end
u_c(u_c=="zf_Glut_exclude")=[];
u6cell=u_c;
% rename variables and save new sorted
n_data=a6data;
n_celltype=ca6;
n_sample=sa6;
n_cellid=ia6;

save('n_sorted.mat','n_celltype','n_sample','X','out_idx','T_cells_x','n_data','all_flags_sorted','n_cellid','s_xy_sorted','u_c','-v7.3');
clear data_orig_aull_sorted geneid_all c T_cells_tmp sample_sorted
%% or load glut 
load('/data/Technion_analysis/zebrafish/sc_100410/SLC17A6A/n_sorted.mat','n_celltype','n_sample','X','out_idx','T_cells_x','n_data','all_flags_sorted','n_cellid','s_xy_sorted','u_c');

%% gad zf : 1
disp('gad2')

% nrmz=20000;

cd('/data/Technion_analysis/zebrafish/sc_100410/SLC32A1')
load('Sorted.mat', 'data_orig_all_sorted','T_cells_tmp','sample_sorted','all_flags_sorted','s_xy_sorted','cellid_sorted')
load('Orignal.mat', 'geneid_all')
gaddata=data_orig_all_sorted;%round(data_orig_all_sorted./repmat(sum(data_orig_all_sorted),length(data_orig_all_sorted(:,1)),1)*nrmz);
txt=readtable('clusters names gaba.xlsx');
gad_name=table2array(txt(:,end));
% add gaba before each cell type name
gad_name=insertBefore(gad_name,1,'zf_GABA_');

sgad=sample_sorted;

igad=cellid_sorted;

% group simiilar clusters
[u_c,~,X] = unique(gad_name,'stable');
cgad=strings(length(T_cells_tmp),1);
T_cells_x=zeros(size(T_cells_tmp));
for cl=1:length(gad_name)
    cgad(T_cells_tmp==cl)=string(u_c(X(cl)));
    T_cells_x(T_cells_tmp==cl)=X(cl);
end
% remove exclude names cluster
out_idx=cgad=="zf_GABA_exclude";
sgad(out_idx)=[];
gaddata(:,out_idx)=[];
cgad(out_idx)=[];
T_cells_x(out_idx)=[];
all_flags_sorted(out_idx,:)=[];
s_xy_sorted(out_idx,:)=[];
igad(out_idx)=[];

bgad_genes=geneid_all;
fgad=all_flags_sorted;
X=unique(X);
X(X==find(u_c=="zf_GABA_exclude"))=[];
u_c(u_c=="zf_GABA_exclude")=[];
ugadcell=u_c;
% rename variables and save new sorted
n_data=gaddata;
n_celltype=cgad;
n_sample=sgad;
n_cellid=igad;

save('n_sorted.mat','n_celltype','n_sample','X','out_idx','T_cells_x','n_data','all_flags_sorted','n_cellid','geneid_all','s_xy_sorted','u_c','-v7.3');

clear data_orig_all_sorted  c T_cells_tmp sample_sorted
%% or load gaba 
load('/data/Technion_analysis/zebrafish/sc_100410/SLC32A1/n_sorted.mat','n_celltype','n_sample','X','out_idx','T_cells_x','n_data','all_flags_sorted','n_cellid','s_xy_sorted','u_c');

%% (just run for decoy names) concatante all nourons for bdata  g76 (optional)
gads=1;
disp('conc ng&7&6')
if gads==1 % gad
cd('/data/Technion_analysis/zebrafish/sc_100410/comparative/GABA')
b_genes=bgad_genes;
ball_data=[gaddata];
ball_c=[cgad];
ball_s=[sgad];
ball_f=[fgad];
bu_c=[ugadcell];
elseif gads==0 % glut
cd('/data/Technion_analysis/zebrafish/sc_100410/comparative/GLUT')
b_genes=ba6_genes;
ball_data=[a6data];
ball_c=[ca6];
ball_s=[sa6];
ball_f=[fa6];
bu_c=[u6cell];
else % combined 
cd('/data/Technion_analysis/zebrafish/sc_100410/')

% [inBoth,xi,xf] = intersect(ba7_genes, ba6_genes) ;
% a6data=a6data(xf,:);
% a7data=a7data(xi,:);
% sa6=sa6(xf,:);
% sa7=sa7(xi,:);
geneid67=ba6_genes;
a76data=a6data;
c76=ca6;
s76=sa6;
f76=fa6;
u76cell=u6cell;
[inBoth,xi,xf] = intersect(geneid67, bgad_genes) ;
a76data=a76data(xi,:);
gaddata=gaddata(xf,:);

% % s67=s67(xi,:);
% % sgad=sgad(xf,:);

b_genes=inBoth;
ball_data=[gaddata,a76data];
ball_c=[cgad;c76];
ball_s=[sgad;s76];
ball_f=[fgad;f76];
bu_c=[ugadcell;u76cell];
end
save('b_g76.mat','ball_f','ball_s','ball_c','ball_data','b_genes','bu_c','-v7.3')
%% update zf genes to account for a & b
disp('conc m&b')
zf_genes=string(b_genes);
for i=1: length(m_genes)
    i
    xg= startsWith(b_genes,upper(m_genes(i)));
    zf_genes(xg)=upper(m_genes(i));
end
[Cx, iax,iac] = unique(upper(zf_genes),'stable','row');

% and combine similar genes (and matrix rows)
zf_genes=Cx;

dats=zeros(length(Cx),size(ball_data,2));
for gpx=1:size(ball_data,2)
    dats(:,gpx)=accumarray(iac,ball_data(:,gpx));
end
ball_data=dats;
%% random sampling gf
if sp==1
disp('sample gf')
y_all=[];
n_g=unique(ball_c,'stable');
for rs=1:length(n_g)
    ig=find(ball_c==n_g(rs));% for a cell cluster
    if  length(ig)>200
        rs
        y = randsample(ig,200); % here we take only XX cells of the real cell typoe cluster
        y_all=[y_all;y];
    else
        y_all=[y_all;ig];
    end

end
g_name=ball_c(y_all);
g_ca=ball_c(y_all);
g_sample=ball_s(y_all);
g_data=ball_data(:,y_all);
end

%% concatante MM & B mutual genes - for mouse bat
disp('conc m&g')

[inBoth,xi,xf] = intersect(upper(m_genes), zf_genes) ;
m_data=m_data(xi,:);
m_data=double(round(m_data./repmat(sum(m_data),length(m_data(:,1)),1)*nrmz));
g_data=g_data(xf,:);
g_data=double(round(g_data./repmat(sum(g_data),length(g_data(:,1)),1)*nrmz));
geneid=inBoth;


 %% neurotasmitor flag
disp('m g76 - nt flag')
if gads==1
gluts=0;
else
gluts=1;
end
% ClusterName=m_name;% CLUSTER NAME
ntgd=contains(m_nt,'GABA');
ntach=contains(m_nt,'Acetylcholine');
ntd=contains(m_nt,'Dopamine');
ntv=contains(m_nt,'Glutamate');
ntz=contains(m_name,"INH");
nty=contains(m_name,"GLU");
nnn=ntgd & ntv;
% nnn=~ntv & ~ ntgd;
ntf=zeros(length(m_nt),1);% 0=> NN
ntf(ntgd | ntach | ntd & ~nty)=gads;% gad a& ach
ntf(ntv & ~ntz)=gluts; % glut

% nti(nnn)=3;
ntf(nnn)=0;
flag_m=ntf;

% only tmp
% nti=ones(length(m_name),1);
% flag_i=zeros(length(m_name),1);
%% create mouse m76 flag
disp('m g76')

% ntgd= m_rank=='Cerebellumneurons';
% ntach= m_rank=='Hindbrainneurons';
% ntv= m_rank=='Spinalcordneurons';
% nttp= m_rank=='Telencephalon projecting neurons';
% ntti= m_rank=='Telencephalon interneurons';
% nttx= m_rank=='Cholinergic, monoaminergic and peptidergic neurons';
% ntty= m_rank=='Di- and mesencephalon neurons';
nttz= contains(m_region,"Hypothalamus");% !!!!!!!!!!!!!!!!!!!!!

% nti=ones(length(m_rank),1);
% nti(ntgd | ntach | ntv)=0;% 0 = filter out

%
nti=zeros(length(m_rank),1);
nti(nttz)=1;% 0 = filter in 
% choose INH /GLU or uncomment 
% ntf=contains(m_name,"INH");
% flag_m=ntf;

ntx=logical(nti) & logical(ntf);
%% create gf flag
disp('gf g76')

ntgd=contains(g_ca,'GABA');
ntgl=contains(g_ca,'Glut');
% nnn=~ntv & ~ ntgd;
flag_gf=zeros(length(g_ca),1);% 0=> NN
flag_gf(ntgd)=gads;% gad a& ach
flag_gf(ntgl)=gluts;% gad a& ach

%% concatante
disp('conc all')
all_name=[m_name(logical(ntx));g_name(logical(flag_gf))];
all_data=[m_data(:,logical(ntx)),g_data(:,logical(flag_gf))];
sampleid=[m_sample(logical(ntx));g_sample(logical(flag_gf))];
n_m = unique(m_name(logical(ntx)),'stable');
n_g = unique(g_name(logical(flag_gf)),'stable');
%% decoy ids & flags
disp('dec flg')
cellid=cell(length(all_name),1);
flag_mgf=2.*ones(length(all_name),1);% gf=2 | b=1
flag_mgf(1:sum(ntx))=1;
flag_g76=[flag_m(ntx);flag_gf(logical(flag_gf))]; % 0 1 2
flag_mgfg76=[flag_m(ntx);flag_gf(logical(flag_gf))+3];% +3 to make indexes 0 1 2  | 3 4
n_mg=[n_m;n_g];
flag_rgn=[m_region(ntx);strings(sum(logical(flag_gf)),1)];
flag_loc=[m_location(ntx);strings(sum(logical(flag_gf)),1)];
flag_gen=[m_marker(ntx);strings(sum(logical(flag_gf)),1)];

%% save
disp('save')

save("n_mzf_10x.mat","all_name","all_data","geneid","flag_mgf","flag_g76","flag_mgfg76","cellid","sampleid",'n_mg','n_m','n_g','flag_rgn','flag_loc','flag_gen','-v7.3')
