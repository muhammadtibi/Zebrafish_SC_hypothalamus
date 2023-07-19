
%% gf_vs_m
clear all
close all
clc
nrmz=5000;
sp=1;
gads=2; % 2 is all

%% load mouse-new (my)
disp('mm')
load('/bigdata/all_mouse_dataset/loom_n.mat','mn_matrix','mn_genes','mn_region','mn_marker','mn_location','mn_cl','mn_cln','mn_id','mn_dev','mn_ts','mn_nt','mn_t3');
m_nt=mn_nt;
m_data=mn_matrix';
mx_data=m_data;
m_name=mn_cln;
m_genes=mn_genes;
m_sample=ones(length(mn_id),1);%mn_id;
m_region=mn_region;
m_location=mn_location;
m_marker=mn_marker;
m_rank=string(mn_t3);
mmid=strcat("mm_",string(1:length(m_name)));

%% Romanov
disp('mr-new')

% read orignal xlsx
% fname='/data/Technion_analysis/zebrafish/sc_100410/comparative/hypoth_moldata_classification08-Mar-2017.xlsx';
% opts = detectImportOptions(fname);
% opts = setvartype(opts, 1,'string');
% hypo2=readtable(fname);
% opts = detectImportOptions(fname);
% opts.VariableTypes = "string";
% hypo2=loadCellFile_csv_turbo(fname);
% R = readcell(fname);
% hypo2=readtable('/data/Technion_analysis/zebrafish/sc_100410/comparative/hypoth_moldata_classification08-Mar-2017.xlsx');
% or load
load('/data/Technion_analysis/zebrafish/sc_100410/comparative/hyporoman.mat','hypog','hypom')
load('/data/Technion_analysis/zebrafish/sc_100410/comparative/hypof.mat','hypof')
hypon=(hypof(2,:)=="neurons");
% % m_nt=mn_nt;
mr_data=hypom(11:end,hypon(2:end));
mr_data=mr_data(:,2:end);
mr_name=hypof(3,hypon(2:end));
mr_name=mr_name(:,2:end)';
mr_genes=upper(hypog(11:end));
mr_sample=3*ones(1,size(mr_data,2))';%mn_id;
mr_region    = cell(size(mr_data,2),1);
mr_region(:) = {'Hypothalamus'};
mr_location=string(3*ones(1,size(mr_data,2)))';
mr_marker=string(3*ones(1,size(mr_data,2)))';
mr_rank=string(3*ones(1,size(mr_data,2)))';
mr_nt=string(3*ones(1,size(mr_data,2)))';
clear hypof hypom hypog
%% intersect m_data with m_ramon
disp('conc m&m')

[inBoth,xi,xf] = intersect(upper(m_genes), mr_genes);
m_data=m_data(xi,:);
mr_data=mr_data(xf,:);
m_genes=inBoth;

m_nt=[m_nt;mr_nt];
m_data=[m_data,mr_data];
m_name=[m_name;mr_name];
m_sample=[m_sample;mr_sample];%mn_id;
m_region=[m_region;mr_region];
m_location=[m_location;mr_location];
m_marker=[m_marker;mr_marker];
m_rank=[m_rank;mr_rank];
mmid=strcat("mm_",string(1:length(m_nt)));
save('/data/Technion_analysis/zebrafish/sc_100410/comparative/mm_ra.mat','m_genes','m_data','-v7.3')
%% random sampling mm
if sp==1
disp('sample mm')
n_m=unique(m_name,'stable');
y_all=[];
for rs=1:length(n_m)
    ig=find(m_name==n_m(rs));% for a cell cluster
    if  length(ig)>200 % if high no sampling
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
mx_data=m_data(:,y_all);
mmid=mmid(y_all);
end
if ~exist('mm_samples.mat')
    save('mm_samples.mat','y_all')
else
    disp('WARNIING USED THE PREVIOUS SAMPLE')
end
%%
if gads==0
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
    txtz=readtable('clusters names glut.xlsx');
    a6_name=table2array(txtz(:,end));% length is number of clusters
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
    for cl=1:length(u_c)
        ca6(n_celltype==u_c(cl))=  insertAfter(ca6(n_celltype==u_c(cl)),"zf_",[num2str(cl),'_']);
    end
    save('n_sorted.mat','n_celltype','n_sample','X','out_idx','T_cells_x','geneid_all','n_data','all_flags_sorted','ba6_genes','n_cellid','s_xy_sorted','u_c','-v7.3');
    clear data_orig_aull_sorted geneid_all c T_cells_tmp sample_sorted
    %% or load glut
    % load('/data/Technion_analysis/zebrafish/sc_100410/SLC17A6A/n_sorted.mat','n_celltype','n_sample','ba6_genes','X','out_idx','T_cells_x','n_data','all_flags_sorted','n_cellid','s_xy_sorted','u_c');
elseif gads==1 % gad 
    %% gad zf : 1
    disp('gad2')

    % nrmz=20000;

    cd('/data/Technion_analysis/zebrafish/sc_100410/SLC32A1/')
    load('Sorted.mat', 'data_orig_all_sorted','T_cells_tmp','sample_sorted','all_flags_sorted','s_xy_sorted','cellid_sorted')
    load('Orignal.mat', 'geneid_all')
    gaddata=data_orig_all_sorted;%round(data_orig_all_sorted./repmat(sum(data_orig_all_sorted),length(data_orig_all_sorted(:,1)),1)*nrmz);
    txtz=readtable('clusters names gaba.xlsx');
    gad_name=table2array(txtz(:,end));
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
    for cl=1:length(u_c)
        cgad(n_celltype==u_c(cl))=  insertAfter(cgad(n_celltype==u_c(cl)),"zf_",[num2str(cl),'_']);
    end
    n_sample=sgad;
    n_cellid=igad;
    save('n_sorted.mat','n_celltype','n_sample','X','out_idx','T_cells_x','bgad_genes','n_data','all_flags_sorted','n_cellid','geneid_all','s_xy_sorted','u_c','-v7.3');

    clear data_orig_all_sorted  c T_cells_tmp sample_sorted
else % all (gads+gluts)
       %% a6 zf : 3
    disp('zf')
    disp('a6')
    cd('/data/Technion_analysis/zebrafish/sc_100410/SLC17A6A/after')
    load('Sorted.mat', 'data_orig_all_sorted','T_cells_tmp','sample_sorted','all_flags_sorted','s_xy_sorted','cellid_sorted')
    load('Orignal.mat', 'geneid_all')


   
    all_flags_sorted=string(zeros(length(sample_sorted),1));
    a6data=data_orig_all_sorted;%round(data_orig_all_sorted./repmat(sum(data_orig_all_sorted),length(data_orig_all_sorted(:,1)),1)*nrmz);
    sa6=sample_sorted;
    ia6=cellid_sorted;
    ba6_genes=geneid_all;
    fa6=all_flags_sorted;
    u_c=extractBefore(ia6,'_10X');
    ca6=regexprep(u_c,'_','-');
    ca6=strcat(string(T_cells_tmp),'-',ca6);
    u6cell=unique(ca6);



 

    %% gad zf : 1
    disp('gad2')

    % nrmz=20000;

    cd('/data/Technion_analysis/zebrafish/sc_100410/SLC32A1/after')
    load('Sorted.mat', 'data_orig_all_sorted','T_cells_tmp','sample_sorted','all_flags_sorted','s_xy_sorted','cellid_sorted')
    load('Orignal.mat', 'geneid_all')
    gaddata=data_orig_all_sorted;%round(data_orig_all_sorted./repmat(sum(data_orig_all_sorted),length(data_orig_all_sorted(:,1)),1)*nrmz);


    sgad=sample_sorted;

    igad=cellid_sorted;

    % group simiilar clusters
    
  

    bgad_genes=geneid_all;
    fgad=all_flags_sorted;
    u_c=extractBefore(igad,'_10X');
    cgad=regexprep(u_c,'_','-');
    cgad=strcat(string(T_cells_tmp),'-',cgad);
    ugadcell=unique(cgad);
end
%% (just run for decoy names) concatante all nourons for bdata  g76 (optional)
disp('conc ng&7&6')

if gads==1 % gad
    cd('/data/Technion_analysis/zebrafish/sc_100410/comparative/GABA')
    b_genes=bgad_genes;
    ball_data=[gaddata];
    ball_datax=[gaddata];
    ball_c=[cgad];
    ball_s=[sgad];
    ball_f=[fgad];
    bu_c=[ugadcell];
elseif gads==0 % glut
    cd('/data/Technion_analysis/zebrafish/sc_100410/comparative/GLUT')
    b_genes=ba6_genes;
    ball_datax=[a6data];
    ball_data=[a6data];
    ball_c=[ca6];
    ball_s=[sa6];
    ball_f=[fa6];
    bu_c=[u6cell];
else % combined
    cd('/data/Technion_analysis/zebrafish/sc_100410/comparative/COMBINED')

    geneid67=ba6_genes;
    a76data=a6data;
    c76=ca6;
    s76=sa6;
    f76=fa6;
    u76cell=u6cell;
    [inBoth,xi,xf] = intersect(geneid67, bgad_genes) ;
    a76data=a76data(xi,:);
    gaddata=gaddata(xf,:);

    b_bar=[];
    b_genes=inBoth;
    ball_datax=[gaddata,a76data];
    ball_data=[gaddata,a76data];
    ball_c=[cgad;c76];
    ball_s=[sgad;s76];
    ball_f=[fgad;f76];
    bu_c=[ugadcell;u76cell];
end
zfid=strcat("zf_",string(1:length(ball_s)));
save('b_g76.mat','zfid','ball_f','ball_s','ball_c','ball_data','b_genes','bu_c','-v7.3')

%% update zf genes to account for a & b
disp('conc m&b')
% with a & b


% and combine similar genes (and matrix rows)
% with ensembel list :
zf2m=table2array(readtable('/data/Technion_analysis/zebrafish/sc_100410/comparative/zf2m.txt'));
[cc, aa,bb]=intersect(upper(b_genes),upper(zf2m(:,5)));
ball_data=ball_data(aa,:);
zf_genes=zf2m(bb,6);
[Cx, iax,iac] = unique(upper(zf_genes),'stable','row');
zf_genes=Cx;

dats=zeros(length(Cx),size(ball_data,2));
for gpx=1:size(ball_data,2)
    dats(:,gpx)=accumarray(iac,ball_data(:,gpx));
end
ball_data=dats;
%% random sampling zf
if sp==1
    disp('sample gf')
    if 1%~exist('mm_samples.mat')
    y_all=[];
    n_g=unique(ball_c,'stable');
    for rs=1:length(n_g)
        ig=find(ball_c==n_g(rs));% for a cell cluster
        if  length(ig)>200 % if high then no sampling
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
    zfid=zfid(y_all);
    end
else
    disp('WARNIING USED THE PREVIOUS SAMPLE')
end
%     g_name=ball_c(:);
%     g_ca=ball_c(:);
%     g_sample=ball_s(:);
%     g_data=ball_data(:,:);
%% concatante MM & B mutual genes - for mouse bat
disp('conc m&g')

[inBoth,xi,xf] = intersect(upper(m_genes), zf_genes) ;
m_data=mx_data(xi,:);
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
nt3=contains(m_nt,'3');

% ntz=contains(m_name,"INH");
% nty=contains(m_name,"GLU");
nnn=ntgd & ntv;
% nnn=~ntv & ~ ntgd;
ntf=zeros(length(m_nt),1);% 0=> NN
ntf(ntgd | ntach | ntd )=gads;% gad a& ach
ntf(ntv)=gluts; % glut
ntf(nt3)=1; % added to add ramon!!!!!!!
% nti(nnn)=3;
% ntf(nnn)=0;
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
mmid=mmid(ntx);
%% create gf flag
disp('gf g76')
ntgd=contains(g_ca,'GABA');
ntgl=contains(g_ca,'Glut');
% nnn=~ntv & ~ ntgd;
flag_gf=zeros(length(g_ca),1);% 0=> NN
flag_gf(ntgd)=gads;% gad a& ach
flag_gf(ntgl)=gluts;% gad a& ach
zfid=zfid(logical(flag_gf));

%% concatante
disp('conc all')
all_name=[m_name(logical(ntx));g_name(logical(flag_gf))];
all_data=[m_data(:,logical(ntx)),g_data(:,logical(flag_gf))];
sampleid=[string(m_sample(logical(ntx)));g_sample(logical(flag_gf))];
n_m = unique(m_name(logical(ntx)),'stable');
n_g = unique(g_name(logical(flag_gf)),'stable');
%% decoy ids & flags
disp('dec flg')
cellid=cell(length(all_name),1);
flag_mgf=2.*ones(length(all_name),1);% gf=3
% flag_mgf(1:(sum(ntx)-length(mr_name)))=1; %  mm=1
% flag_mgf(1+(sum(ntx)-length(mr_name)):sum(ntx))=2; %  mr=2
flag_mgf(1:sum(logical(ntx)))=1; %  mm=1
flag_g76=[flag_m(ntx);flag_gf(logical(flag_gf))]; % 0 1 2
flag_mgfg76=[flag_m(ntx);flag_gf(logical(flag_gf))+3];% +3 to make indexes 0 1 2  | 3 4
n_mg=[n_m;n_g];
flag_rgn=[m_region(ntx);strings(sum(logical(flag_gf)),1)];
flag_loc=[m_location(ntx);strings(sum(logical(flag_gf)),1)];
flag_gen=[m_marker(ntx);strings(sum(logical(flag_gf)),1)];

%% save
disp('save')
save("n2_mzf_10x.mat","zfid","mmid","all_name","all_data","geneid","flag_mgf","flag_g76","flag_mgfg76","cellid","sampleid",'n_mg','n_m','n_g','flag_rgn','flag_loc','flag_gen','gads','-v7.3')

%% co-scatter species -percent expression

id=find(strcmpi(geneid,'gad2'))
md=(mean(all_data(:,flag_mgf==1)>0,2));
gd=(mean(all_data(:,flag_mgf==2)>0,2));
figure('color','w');
scatter(md,gd)
%% co-scatter genes-species -percent expression
%% load
close all
clear all
clc

genedublex=["NPVF","HCRT"];
load('/data/Technion_analysis/zebrafish/sc_100410/comparative/mm_ra.mat','m_genes','m_data')
load('/data/Technion_analysis/zebrafish/sc_100410/comparative/COMBINED/b_g76.mat','ball_data','b_genes')
%% plot HCRT VS NPVF: 
% old 

% zfx=startsWith(all_name,'zf_');
% id=find(strcmpi(geneid,genedublex(1)));
% md=log2(all_data(id,zfx==0)+1);
% gd=log2(all_data(id,zfx==1)+1);
% id2=find(strcmpi(geneid,genedublex(2)));
% md2=log2(all_data(id2,zfx==0)+1);
% gd2=log2(all_data(id2,zfx==1)+1);

% new
id=find(strcmpi(m_genes,genedublex(1)));
idx=find(strcmpi(b_genes,genedublex(1)));

md=log2(m_data(id,:)+1);
gd=log2(ball_data(idx,:)+1);
id2=find(strcmpi(m_genes,genedublex(2)));
id2x=find(strcmpi(b_genes,genedublex(2)));
md2=log2(m_data(id2,:)+1);
gd2=log2(ball_data(id2x,:)+1);
figure('color','w');
% zf
scatter(gd,gd2,30,[104,188,134]/255,'filled')
hold on
% mm
scatter(md,md2,30,[244,88,226]/255,'filled')

xlabel(genedublex(1))
ylabel(genedublex(2))
legend({'Zebrafish','Mouse'})
% refline
title("Log2 expression")
axis tight;
%  xlim([0 10])
%  ylim([0 10])
% hline=refline;
% hline.Color = 'k';
%% plot Venn of geneduplex
figure('color','w');
subplot(1,2,1)
A = [sum(logical(gd)), sum(logical(gd2))]; I = sum(logical(gd) & logical(gd2));
[~,S]=venn(A,I,'FaceColor',{'y','c'},'EdgeColor','black');
axis equal; axis tight;axis off;
title('ZF(logical>0)')
text(S.ZoneCentroid(1,1), S.ZoneCentroid(1,2),genedublex(1));
text(S.ZoneCentroid(2,1), S.ZoneCentroid(2,2),genedublex(2));
% figure('color','w');
subplot(1,2,2)
A = [sum(logical(md)), sum(logical(md2))]; I = sum(logical(md) & logical(md2));
[~,S]=venn(A,I,'FaceColor',{'y','c'},'EdgeColor','black');
axis equal; axis tight;axis off;
text(S.ZoneCentroid(1,1), S.ZoneCentroid(1,2),genedublex(1));
text(S.ZoneCentroid(2,1), S.ZoneCentroid(2,2),genedublex(2));
title('MM (logical>0)')
