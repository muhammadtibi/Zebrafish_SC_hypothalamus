%% SCRATCH: code for scRNAseq and visseq
% for further clustring after csubcalss clustring !!!! READ THIIS LINE !!!!
% UPDATED: GO TO TO DO COMBINED AFTER SCRATCH TSNE 
%% construction
% fprintf(2,'Code under-construction, please try again later \nMuhammad \n')
% return;
%% bgn
clear all
warning('off','all')
uuu=cd;
cd('/data/runs/samples/')
set(0,'DefaultFigureWindowStyle','normal')

%% get parameters

try
    transcendentalDoc = '';
    clc
    S = {...
        struct('name','Cluster Type','type','enum','values',{{'CatClust','cClust','vClust'}},'doc',transcendentalDoc);...
        struct('name','Show(first=counter)','type','str','default','Snap25,Stmn2,Gad2,Slc32a1,Slc17a7,Slc17a6,Sst,Sim1,Foxj1,Pdgfra,Mog,C1qc,Flt1,Cldn5,Aqp4,Plp1,');...
        struct('name','Group Genes','type','str','default','C1qc,C1qa,C1qb,Gja1,Cx3cr1,Acta2,Ly6c1,Mfge8,Plp1,Aqp4,Vtn,Cldn5,Pdgfrb,Flt1,Slc1a3,Pdgfra,Foxj1,Olig1,Olig2,Sox10,Hbb-bs,Hbb-bt,Hba-a2,');...
        struct('name','exclude genes','type','str','default','Xist,Tsix,Eif2s3y,Ddx3y,Uty,Kdm5d,Btg2,Jun,Egr4,Fosb,Junb,Gadd45g,Fos,Arc,Nr4a1,Npas4,Coq10b,Tns1,Per2,Ptgs2,Rnd3,Tnfaip6,Srxn1,Tiparp,Ccnl1,Mcl1,Dnajb5,Nr4a3,Fosl2,Nptx2,Rasl11a,Mest,Sertad1,Egr2,Midn,Gadd45b,Dusp6,Irs2,Plat,Ier2,Rrad,Tpbg,Csrnp1,Peli1,Per1,Kdm6b,Inhba,Plk2,Ifrd1,Baz1a,Trib1,Pim3,Lrrk2,Dusp1,Cdkn1a,Pim1,Sik1,Frat2,Dusp5,VGF');...
        struct('name','Category Genes(1|2+3)','type','str','default','Gad2,Slc17a7,Slc17a6,');...
        struct('name','Normalization','type','int','default',30000);...
        struct('name','DBSCAN MinPts','type','int','default',30);...
        struct('name','DBSCAN Eps.|K','type','int','default',70);...
        struct('name','Valid  min.gen','type','int','default',1000);...
        struct('name','Valid  min.mol','type','int','default',1000);...
        struct('name','Valid  max.mol','type','int','default',5e4);...
        struct('name','Heatmap top genes','type','int','default',0);...
        struct('name','compare factor','type','int','default',1.5);...
        struct('name','Save fig and pdf','type','checkbox','default',1);...
        struct('name','Use Loaded','type','checkbox');...
        struct('name','Plot Flags','type','checkbox','default',1);...
        struct('name','Upper Genes','type','checkbox','default',1);...
        struct('name','Stops','type','checkbox','default',0);...
        struct('name','Load 10XMTX','type','checkbox','default',1);...
        struct('name','NAN|dbl|zero out','type','checkbox','default',1);...
        struct('name','Class non-neuro out','type','checkbox');...
        struct('name','Batch Correction','type','int','default',0);...
        struct('name','Fet. Sel.','type','enum','values',{{'Default','Alt','Spatial','InterSpecies'}});...
        struct('name','Cluster Method','type','enum','values',{{'DBSCAN','KNN-DBSCAN','K-means','Linkage','KNN-louvain'}});...
        struct('name','2D Embedding','type','enum','values',{{'t-SNE','UMAP','Phate','PCA'}});...
        struct('name','Filter(-out|+in)','type','int','default',0);...
        };
    
    %
    Param = Settings_GUI(S);
    clusty=string(Param(1));
    shower=string(Param(2));
    shower=strsplit(shower,',');
    shower=shower(1:end-1);
    list=convertStringsToChars(shower);
    if isempty(list)
        list={};
    end
    showerg=string(Param(3));
    showerg=strsplit(showerg,';');
    showerg=showerg(1:end-1);
    if ~isempty(showerg)
        groupy=1;
        list2={};
    else
        groupy=0;
        shower1=string(Param(3));
        shower1=strsplit(shower1,',');
        shower1=shower1(1:end-1);
        list2=convertStringsToChars(shower1);
        if isempty(list2)
            list2={};
        end
    end
    Excluder=string(Param(4));
    Excluder=strsplit(Excluder,',');
    Excluder=Excluder(1:end-1);
    ex_genes=convertStringsToChars(Excluder);
    if isempty(ex_genes)
        ex_genes={};
    end
    Exclust=string(Param(5));
    Exclust=strsplit(Exclust,',');
    Exclust=Exclust(1:end-1);
    ex_clust=string(convertStringsToChars(Exclust));
    if isempty(char(ex_clust))
        ex_clust={};
    end
    nrmz=cell2mat(Param(6));
    MinPts=cell2mat(Param(7));
    eps_prc=cell2mat(Param(8));
    valg=cell2mat(Param(9));
    valm=cell2mat(Param(10));
    valmx=cell2mat(Param(11));
    top_genes=cell2mat(Param(12));
    fct=cell2mat(Param(13));
    savefig_flag=cell2mat(Param(14));
    if savefig_flag
        savefig_pdf=1 ;
    end
    previous=cell2mat(Param(15));
    flags=cell2mat(Param(16));
    upgn=cell2mat(Param(17));
    if upgn==1
        list=upper(list);
        list2=upper(list2);
        ex_clust=upper(ex_clust);
        ex_genes=upper(ex_genes);
    end
    stops=cell2mat(Param(18));
    multicrops=cell2mat(Param(19));
    nndpzr=cell2mat(Param(20));
    classy=cell2mat(Param(21));
    batcorr=cell2mat(Param(22));
    fetsel=cell2mat(Param(23));
    clmethod=char(Param(24));
    if clmethod ~= "DBSCAN"
        minpts=0;
    end
    embd2d=string(Param(25));
    filtery=cell2mat(Param(26));
catch ME
    close all
    cd(uuu);
    return ;
end
direct=uigetdir('/data/Technion_analysis','Please choose MAIN analysis folder');% get save folder

close all
clc

%% loadings and swithes
if 1 % to run all sections
    if clusty=="CatClust"
        cd(direct)
        cd ..
        disp('Loading Sorted... please wait...')
        load('n_sorted.mat');
        MinPtsx=MinPts;
        eps_prcx=eps_prc;
        load('Clust_param.mat');
        MinPts=MinPtsx;
        eps_prc=eps_prcx;
        load('Orignal.mat', 'geneid_all')
        load('Orignal.mat', 'bar_ar')
        
        data_org=n_data;
        cellid_org=n_celltype;
        geneid_org=geneid_all;
        sample_org=n_sample;

        geneid=geneid_org;
        %     if grp~=0
        %         groups=1:grp
        %         grpn = input('Which Group?');
        %         cellid_gaba_glut = loadCellFile([num2str(grpn),'_exclust.txt']);
        %         gena = input('1=Genes/2=NoGenes');
        %         nclust;
        %         nclust;
        
   
        data = data_org;
        cellid = cellid_org;
        sample = sample_org;
        % data = round(data./repmat(sum(data),length(data(:,1)),1)*nrmz);
    end
    
    % if %vClust | cClust %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % chose flags
    if clusty=="vClust"
        are_samples= readcell('/bigdata/microscope_images/vis_samples.csv');
    else
        are_samples= readcell('/bigdata/microscope_images/10X_samples.csv');
        
    end
    if clusty=="CatClust" % create folder for catclust specific gene
    
        cd(direct)
        
        
        
        
        var_names=are_samples(1,:);
        are_samples(1,:)=[];
        are_sample=table(are_samples);
        flist=var_names(3:6);
    else
        var_names=are_samples(1,:);
        are_samples(1,:)=[];
        are_sample=table(are_samples);
        flist=var_names(3:6);
        [findx,~] = listdlg('ListString',flist,'SelectionMode','multiple','Name','Flags Selection','ListSize',[150,70]);
        if isempty(findx)
            return;
        end
        
        
        %%
        % loop through files
        if previous==1
            fprintf('Please Select '); fprintf(2, 'Sample Mother Folder/s \n');
            cd('/data/runs/samples')
            listyx=uipickfiles('num',[],'FilterSpec',[],'out','struct'); % loop images
            listy=listyx;
            clc
         
        else
            
        end
        fprintf('Please Select '); fprintf(2, 'Sample Mother Folder/s \n');
        cd('/data/runs/samples')
        listyx=uipickfiles('num',[],'FilterSpec',[],'out','struct'); % loop images
        listy=listyx;
        clc
        cd(direct)
        
        geneidx={};datax={};barcodesx={};sample=[];data=[];cellid=[]; bar_ar={};
        disp('Loading matrix files..')
        for lis=1:length(listy)
            fprintf([num2str(lis),' / ',num2str(length(listy))]);
            
            path=listy(lis).name;
            diry=strsplit(path,'/');
            cd(path)
            if  multicrops==0
                outs_files = dir([path filesep 'out*']); % get the contents of the image_folder
                % %                 if isempty(outs_files)
                % %                     disp('multicrop')
                % %                     cd ..
                % %                     cd ..
                % %                     path=cd;
                % %                 outs_files = dir([path filesep 'out*']); % get the contents of the image_folder
                % %                 outs_files = natsortfiles({outs_files.name});
                % %                 outs_files = char(outs_files(end))
                % %                 else
                outs_files = natsortfiles({outs_files.name});
                outs_files = char(outs_files(end))
                %                 end
                
                % % % % % %             try
                % % % % % %                 zzzzz
                % % % % % %                 load('10xMtxFile.mat') % check if there is already a file and read it
                % % % % % %                 disp('          read')
                % % % % % %             catch
                mtx_file = [path,'/',outs_files,'/outs/filtered_feature_bc_matrix/matrix.mtx'];
                bc_file = [path,'/',outs_files,'/outs/filtered_feature_bc_matrix/barcodes.tsv'];
                gene_file =[path,'/',outs_files,'/outs/filtered_feature_bc_matrix/features.tsv'];
                [datax,geneidx, barcodesx] = load10xMtxFile(mtx_file,bc_file,gene_file,[],[]);
                if nndpzr==1 % no duplication no nans
                    % grouping
                    disp('grouping similar genes...')
                    %                     if lis==1
%                     geneidx=geneidx(:,1);
%   load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_violin/violin/ort_geneid.mat','geneid')

                    geneidx=geneidx(:,1);
                    
                    [Cx, iax,iac] = unique(upper(geneidx),'stable','row');
                    
                    %                     end
                    geneidx=Cx;
                    
                    dats=zeros(length(Cx),size(datax,2));
                    for gpx=1:size(datax,2)
                        dats(:,gpx)=accumarray(iac,datax(:,gpx));
                    end
                    
                    disp('Sorting Out NaNs...')
                    
                    datax=dats;
                    % NANS OUT
                    nany=find(geneidx=="NaN");
                    geneidx(nany,:)=[];
                    datax(nany,:)=[];
                    % chck for locs
                    try
                        TT=readtable('/data/Aligment/Bat_genome/batLOC2Symbol.csv');
                        genelocs=table2array(TT(:,1));
                        disp('Converting LOCS...')
                        alloc=table2array(TT(:,2));
                        idloc=contains(string(geneidx),"LOC",'IgnoreCase',true);
                        sum_locs=sum(idloc)
                        gnloc=geneidx(idloc);
                        %                         fullid=1:length(idloc);
                        %                         fullid=fullid(idloc)';
                        [interg,ia,ib]=intersect(geneidx,alloc); % update genelocs acoording to TT bat annotated table
                        geneidx(ia)=genelocs(ib);
                        %                         for zz=1:length(ib)
                        %                         geneidx(ia(zz))={char(update_gene(zz))};
                        %                         end
                        %                         disp('Sorting Out LOCS...')
                    end
                    nany=contains(string(geneidx),"LOC",'IgnoreCase',true);
                    geneidx(nany,:)=[];
                    datax(nany,:)=[];
                    % grouping again if LOCS are repeated
                    disp('grouping similar genes...')
                    [Cx, iax,iac] = unique(geneidx(:,1),'stable','row');
                    
                    geneidx=Cx;
                    
                    dats=zeros(length(Cx),size(datax,2));
                    for gpx=1:size(datax,2)
                        dats(:,gpx)=accumarray(iac,datax(:,gpx));
                    end
                    
                    
                    datax=dats;
                    
                end
                save('10xMtxFile.mat','datax','geneidx','barcodesx','-v7.3');
            else
                load('10xMtxFile.mat') % check if there is already a file and read it
                disp('          read')
            end
            
            
            
            
            if upgn==1
                geneidx=upper(geneidx(:,1));
            end
            
            
            if batcorr==1 && lis>1 % intersect spceies
                [inBoth,xi,xf] = intersect(geneidx, geneid) ;
                data=data(xf,:);
                datax=datax(xi,:);
                geneid=inBoth;
            else
                geneid=geneidx;
            end
     
             %%%
            sample = [sample;repmat({char(diry(end))},length(barcodesx),1)];
            cellid = [cellid;barcodesx];
            data = [data,datax];
            if lis==1
                if isempty(list)
                    disp('Shower Genes is empty')
                    sainity1=1;
                else
                    sainity1= ismember(list,geneid);
                end
                if isempty(list2)
                    disp('exCategory Genes is empty')
                    sainity2=1;
                else
                    sainity2= ismember(list2,geneid);
                end
                if isempty(ex_clust)
                    disp('Category Genes is empty')
                    sainity3=1;
                else
                    sainity3= ismember(ex_clust,geneid);
                end
                if isempty(ex_genes)
                    disp('excluster Genes is empty')
                    sainity4=1;
                else
                    sainity4= ismember(ex_genes,geneid);
                end
                if ~all(sainity1(:))
                    disp('Shower Genes')
                    list(sainity1==0)
                    disp('Warning: Genes were not found in this file')
                    disp('Make sure the gene name is writen correctly with Capitliazation and comma afterwards.')
                    list(sainity1==0)=[];
                end
                if ~all(sainity2(:))
                    disp('exCategory Genes')
                    list2(sainity2==0)
                    disp('Warning: Genes were not found in this file')
                    disp('Make sure the gene name is writen correctly with Capitliazation and comma afterwards.')
                    list2(sainity2==0)=[];
                end
                if ~all(sainity3(:))
                    disp('Category Genes')
                    ex_clust(sainity3==0)
                    disp('Error: Genes were not found in this file')
                    disp('Make sure the gene name is writen correctly with Capitliazation and comma afterwards.')
                    %             ex_clust(sainty4==0)=[];
                    return;
                end
                if~all(sainity4(:))
                    disp('excluster Genes')
                    ex_genes(sainity4==0)
                    disp('Warning: Genes were not found in this file')
                    disp('Make sure the gene name is writen correctly with Capitliazation and comma afterwards.')
                    ex_genes(sainity4==0)=[];
                end
            end
            disp('XY & Spatial barcodes')
            ln=length(barcodesx);
            str_barcodes=string(barcodesx);
            if clusty=="vClust"
                bar_xy=cell(ln,3);
                if multicrops==0
                    tissue_positions_list=readtable([path,'/',outs_files,'/outs/spatial/tissue_positions_list.csv']);
                    xst=tissue_positions_list.Var5;yst=tissue_positions_list.Var6;
                    bar_org=string(tissue_positions_list.Var1);
                    
                    for ib=1:ln
                        curr_barcode=str_barcodes(ib);
                        bar_ind=find(curr_barcode==bar_org);
                        bar_xy(ib,:)={[char(curr_barcode),'_',char(diry(end))],xst(bar_ind),yst(bar_ind)};
                    end  %barcode
                    mapped_xy=cell2mat(bar_xy(:,2:3));
                    save([path,'/mapped_xy.mat'],'mapped_xy')
                else
                    grp=0;
                    load([path,'/10xMtxFile.mat'], 'barcodesx')
                    load([path,'/mapped_xy.mat'],'mapped_xy')
                    for ib=1:ln
                        bar_xy(ib,:)={[char(barcodesx(ib,1)),'_',char(diry(end))],mapped_xy(ib,1),mapped_xy(ib,2)};
                    end
                end
                bar_ar=[bar_ar;bar_xy];
            end


        end
        %  spatialid=cellid;
        % addpath(genpath(direct));

        %%
           %% exttra for bat and mouse
        disp('extra bg76')
        direct='/data/Technion_analysis/zebrafish/sc_100410/comparative/COMBINED';
        cd(direct);
        load('scratch_starter.mat');
        load('b_g76.mat','b_genes','ball_c','ball_data','ball_f','ball_s','bu_c')
        sampleid=ball_s;
        all_data=ball_data;
        geneid=b_genes;
        cellid=ball_c;
        all_flags(1).fc_time=ball_f;
        %%

        sample=sampleid;
        data=all_data;
        listy={'1'}
        geneid_org=geneid;
        data_org=data;
        cellid_org=cellid;
        sample_org=sample;
        cd(direct)
%         disp('saving loads mat..')
        % % %         save('Loads','data_org','geneid','cellid_org','bar_ar','sample_org','geneid_org','cellid_org','-v7.3');
        
    end
    listy=unique(sample);
    sample_uni=listy;
    
    %% take out NANS and duplications grouping and zeros row
%     if nndpzr==1
%         
%         % zeros out
%         disp('Sorting Out zero rows of Data...')
%         zerow=~any(data,2); % index with zero rows == no gene expresson at all in every cell
%         geneid(zerow, : ) = [];  %zero rows out first genes and then data
%         data(zerow, : ) = [];  %zero rows out
%         %         geneid_org=geneid;
%         %         data_org=data;
%     end
    % catclust or cvclust
    
    %% flags
    %
    %     listyx=listy;
    
    
    % flags to store data structers
    try
    %     if multicrops==0
    if batcorr==1
        project=string(are_samples(string(are_samples(:,2))==string(listy(1)),1));
        samples_id=find(string(are_samples(:,1))==project);
        project=string(are_samples(string(are_samples(:,2))==string(listy(end)),1)); % last project should bethe second list!
        samples_id2=find(string(are_samples(:,1))==project2);
        samples_idx=[samples_id ;samples_id2];
    else
        project=string(are_samples(string(are_samples(:,2))==string(listy(1)),1))
        samples_idx=find(string(are_samples(:,1))==project);
    end
    %     else
    %     listyx=string(listyx(1).name);
    %     diryx=strsplit(listyx,'/');
    %     diryx=char(diryx(end-2));
    %      project=string(are_samples.Project(string(are_samples.SampleID)==string(diryx)))
    %     end
    
    findx=findx+2;
    x_samples=are_samples(samples_idx,:);
    
    for f_i=1:length(findx) % how many conditions
        fc_time = zeros(size(cellid));
        colflag=convertCharsToStrings(x_samples(:,findx(f_i)));
        %         colflag = rmmissing(string(colflag));
        fc_time_uni=unique(colflag);
        for uf_i=1:size(fc_time_uni,1)% how many subconditions
            
            s_idx=find(string(colflag)==string(fc_time_uni(uf_i)));
            
            cond_samples=x_samples(s_idx,2);
            for csid=1:length(cond_samples) %find chosen samples and condition
                f_sample = find(string(sample)==string(cond_samples(csid)));
                fc_time(f_sample)=uf_i;%update fc_times
                all_flags(f_i).fc_time=fc_time;
            end
        end
    end
    
    
    catch 
    
    %%% try catch  temp for mouse&bat
all_flags(1).fc_time=flag_mb; 
all_flags(2).fc_time=all_c; 
all_flags(3).fc_time=flag_g67;
all_flags(4).fc_time=flag_mbg67;
    end
    
    %% distribution
    set(0,'DefaultFigureWindowStyle','docked');
    figure('Name','distribution','NumberTitle','off');
    while 1
        tStart = tic;
        disp("Normalizing...")
        tot_mol = sum(data);
        maximal_total_molecule=max(tot_mol)
        perct_90th_molecule=round(prctile(tot_mol,90));
%         nrmz=perct_90th_molecule % dynamic normalize 
%         tot_mol(tot_mol>nrmz) = nrmz;%%what value?
        tot_genes = sum(data>0);
        set(gcf,'color','w')
        % molecules
        subplot(2,2,1);
        [f,xi] = ksdensity(tot_mol);
        plot(xi,f);
        axis tight;
        set(gca,'xlim',[0,prctile(tot_mol,98)]);
        xlabel('total mol')
        ylabel('kdensity val')
        title('Molecules')
        %
        subplot(2,2,3);
        [f,xi] = ecdf(tot_mol);
        % % % % xisave=xi;fsave=f;
        plot(xi,f);
        axis tight;
        set(gca,'xlim',[0,prctile(tot_mol,98)]);
        xlabel('total mol')
        ylabel('ecdf %')
        % genes
        subplot(2,2,2);
        [f,xi] = ksdensity(tot_genes);
        plot(xi,f);
        axis tight;
        set(gca,'xlim',[0,prctile(tot_genes,98)]);
        xlabel('total genes')
        ylabel('kdensity val')
        title('Genes')
        %
        subplot(2,2,4);
        [f,xi] = ecdf(tot_genes);
        plot(xi,f);
        axis tight;
        set(gca,'xlim',[0,prctile(tot_genes,98)]);
        xlabel('total genes')
        ylabel('ecdf %')
        %
        sgtitle(['pre-valid: #cells:',num2str(length(cellid_org)),' | #genes:',num2str(length(geneid))])
        %% normalization
        
        if stops==1
            shg;
            disp('PRESS any button to continue')
            pause;
            commandwindow;
            redo = input('redo(y)?','s');
            if redo=="y"
                nrmz
                nrmz = input('Normalization:');
                valm
                valm = input('valid Molecules:');
                valg
                valg = input('valid Genes:');
            else
                break;
            end
        else
            break;
        end
    end
    if savefig_flag==1
        savefig(gcf,'distribution.fig') % take property name
    end
    
%     if clusty ~="CatClust"
%     data = normalize(data,'norm');
%     data = round(data./repmat(sum(data),length(data(:,1)),1)*nrmz);
%     end
    %% clear for memory
    % save('Orignal','data_org','cellid_org','sample_org','bar_ar','-v7.3');
    clear  barcodesx datax geneidx
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %%
    % stmn2 = data(strcmpi(geneid,'Stmn2'),:);
    % snap25 = data(strcmpi(geneid,'Snap25'),:);
    % exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
    %     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
    % % % exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plxnb3','Cldn11'....
    % % %     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc25a18','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2','Ttr'};
    % % % [~,loc] = ismember(exclude_markers,geneid);
%     tot_mol=sum(data,2);
    % validcells = (tot_mol>4000 & tot_mol<5e4 & tot_genes>2000 & sum(data(loc,:)>0)<=2 & (stmn2>0 | snap25>0) & amy_flag');
    disp('filtring non-valid cells...')
%     validcells = (tot_mol>valm & tot_mol<valmx & tot_genes>valg); % all cells%5 value??
    validcells=1: size(data,2);
    for i=1:length(sample_uni)
        fprintf(['valid cells in ',sample_uni{i},' = ', num2str(sum(validcells(strcmpi(sample,sample_uni{i})))),'\n']);
    end
    sum(validcells)
    data = data(:,validcells);
    cellid = cellid(validcells);
    sample = sample(validcells);
    sample_org=sample;
    for flaged=1:size(all_flags,2)
        all_flags(flaged).fc_time = all_flags(flaged).fc_time(validcells);
    end
    all_flags_org= [all_flags.fc_time];
    cellid = cellfun(@(x,y) [x,'_',y], cellid, sample,'UniformOutput',0);
    
    % stmn2 = data(strcmpi(geneid,'Stmn2'),:);
    snap25 = data(strcmpi(geneid(:,1),shower(1)),:);
    
    
    
    for i=1:length(sample_uni)
        fprintf([[char(shower(1)),' in before '],sample_uni{i},' = ', num2str(sum(snap25(strcmpi(sample,sample_uni{i}))>0)),'\n']);
    end
    
    %%or load it
    % load afterloading_FC_0_2_8_24_02-Aug-2020
    tot_mol = sum(data);
%     tot_mol(tot_mol>nrmz) = nrmz; % ??why twice
    tot_genes = sum(data>0);
    %
    
    data_orig_all = data;
    geneid_all = geneid;
    
    %% feature selection
    % % % %     disp("Excluding genes...")
    % % % %     ex_genes
    cd(direct)
    
    if classy==0
        fetsel="Default"
                ex_genes=upper({'jun','fos','atf3','egr','agr2','hsp70.1','plk3','btg','ubb','gadd45','MT-','rrad', 'npas4', 'dnajb', 'nr4a1', 'hsps1', 'dusp1', 'dusp5', 'jdp2b','zgc:','si:'});

        disp("Feature selection...")
        in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),ex_genes));
%         in=1:size(data,1);
        if fetsel=="Default"
            corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);  % org
        elseif fetsel=="InterSpecies"
            batchid = all_flags(1).fc_time;%  1 flag for the species !!
            in_m = find(sum(data(:,batchid==1)>0,2)>5 & sum(data(:,batchid==1)>0,2)<length(data(1,batchid==1))*0.5 & ~ismember(geneid(:,1),ex_genes));
            in_b = find(sum(data(:,batchid==2)>0,2)>5 & sum(data(:,batchid==2)>0,2)<length(data(1,batchid==2))*0.5 & ~ismember(geneid(:,1),ex_genes));
            corr_filt1 = cv_vs_m_selection(data(in_m,batchid==1),geneid(in_m),[],1,0);
            corr_filt2 = cv_vs_m_selection(data(in_b,batchid==2),geneid(in_b),[],1,0);
            % %             corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);  % org
            corr_filt=[in_m(corr_filt1);in_b(corr_filt2)];
            corr_filt=unique(corr_filt);
            
        elseif fetsel=="Alt" % new filter
            corr_filt = m_vs_f_selection(data(in,:),geneid(in),[],1,0);
        else % "Spatial"
            pdist_vis
        end
        
        data = data(in(corr_filt),:);
        geneid = geneid(in(corr_filt));
        
        %% gene autocorelation
        
        % batch_genes = {'Fau','Rps27','Rps12','Rpl10','Rpl27','Rpl23a','Pam16','Amd1','Uba52','Rpl35','Rpl7a','Gm10076','Lars2',....
        %     'Cxx1b','Erdr1','Entpd6','AY036118','Eno1','Rplp0','Pigt'};
        % sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};
        disp('auto-genecorlation...')
        maxcl=50;
        geneid_org=geneid;
        data_org=data;
        while 1
            z = linkage(data,'ward','correlation'); % linkage on genes NOT cells== data rows  not col
            
            idx = cluster(z,'maxclust',maxcl);
            
            d = corr_mat(data');
            dsum = sum(d>0.2,2);
            leaforder = optimalleaforder(z,squareform(1-d,'tovector'));
            
            % [~,xi] = sort(idx);
            figure('Name','gcorrelation','NumberTitle','off');
            set(gcf,'color','w')
            h1 = axes('position',[0.2,0.1,0.7,0.85]);
            imagesc(d(leaforder,leaforder),[0,0.3]);
            set(gca,'ytick',[1:length(leaforder)],'YTickLabel',geneid(leaforder))
            freezeColors(h1)
            % % text([1:length(idx)],[1:length(idx)],num2str(idx(leaforder)))
            h2 = axes('position',[0.91,0.1,0.04,0.85]);
            imagesc(idx(leaforder));
            colormap('prism');
            h3 = axes('position',[0.95,0.1,0.05,0.85]);
            plot(sum(d(leaforder,leaforder)>0.2,2),[1:length(leaforder)],'.');
            set(h3,'ydir','reverse');
            linkaxes([h1,h2,h3],'y')
            axis tight;%axis equal;
            axis off;
            
            if savefig_flag==1
                savefig(gcf,[direct,'/gcorr.fig']) % take property name
            end
            hold off
            if stops==1
                shg;
                disp('PRESS any button to continue')
                pause;
                commandwindow;
                redo=input('redo(y)? ','s');
                if redo=='y'
                    maxcl
                    maxcl=input('max clust for linkage:');
                else
                    break;
                end
            else
                break;
            end
        end
        %     strx=lower(geneid(leaforder));
        strx=geneid(leaforder);
        % for making the
        %     sidx=regexp([' ' strx],'(?<=\s+)\S','start')-1
        %     strx(sidx)=upper(strx(sidx))
        %     propermsg = strcat(upper(strx(1)),lower(impropermsg(2:end)));
        table1 = [strx,m2c(idx(leaforder))];
        saveCellFile(table1,[direct,'/genes_correlation.csv'])
    else  %%  class exclusion classy==1
        
        gcorr_table=readcell('genes_correlation.csv');
        disp('excluding non-nouron genes')
        classid=find(gcorr_table(:,3)~="neuron");
        classexclude=gcorr_table(classid,1)';
        if upgn==1
            classexclude=upper(classexclude)
        end
        exy=[ex_genes,classexclude]
        disp("Feature selection...")
        batchid=all_flags(1).fc_time;
% %         data_orig_all = data;
% %         geneid_all = geneid;
        in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),exy));
        if fetsel=="Default"
            corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);  % org
        elseif fetsel=="InterSpecies"
            batchid = all_flags(1).fc_time;%  1 flag for the species !!
            in_m = find(sum(data(:,batchid==1)>0,2)>5 & sum(data(:,batchid==1)>0,2)<length(data(1,batchid==1))*0.5 & ~ismember(geneid(:,1),ex_genes));
            in_b = find(sum(data(:,batchid==2)>0,2)>5 & sum(data(:,batchid==2)>0,2)<length(data(1,batchid==2))*0.5 & ~ismember(geneid(:,1),ex_genes));
            corr_filt1 = cv_vs_m_selection(data(in_m,batchid==1),geneid(in_m),[],1,0);
            corr_filt2 = cv_vs_m_selection(data(in_b,batchid==2),geneid(in_b),[],1,0);
            % %             corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);  % org
            corr_filt=[in_m(corr_filt1);in_b(corr_filt2)];
            corr_filt=unique(corr_filt);
            data = data(in(corr_filt),:);
            geneid = geneid(in(corr_filt));
        elseif fetsel=="Alt" % new filter
            corr_filt = m_vs_f_selection(data(in,:),geneid(in),[],1,0);
            data = data(in(corr_filt),:);
            geneid = geneid(in(corr_filt));
        else % "Spatial"
            pdist_vis%(string(sample),data,200,5,5,string(geneid),bar_ar);
            corr_filt=corr_filt(corr_filt>0);
            data = data(corr_filt,:);
            geneid = geneid(corr_filt);
            in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),exy));
            data = data(in,:);
            corr_filt=corr_filt(in);
            geneid = geneid(in);
            save('corr_filt.mat','corr_filt'); 
            
        end
       
    end
   
    
    cd(direct)
    %%
% %    data = normalize(data,'norm'); % not always for bat & mouse

    disp('PCA...')
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    moldata = data;
    datalog_tmp = cent_norm([(log2(moldata+1))]);
    data_tsne = datalog_tmp;
    % data_tsne = cent_norm(log2(datamarkers+1));
    initial_dims = length(corr_filt);
    [prj,m,D,V,Q] = pca_wis(data_tsne',initial_dims);
    D = diag(D);
    initial_dims = findknee(D);
    figure('Name','PCA','NumberTitle','off');
    set(gcf,'color','w')
    subplot(2,1,1)
    plot(cumsum(D)); hold on;
    plot(initial_dims,sum(D(1:initial_dims)),'sk');
    subplot(2,1,2)
    plot(D); hold on;
    plot(initial_dims,D(initial_dims),'sk');
    title(['opt PC = ',num2str(initial_dims)]);
    % initial_dims = 30;
    prj = prj(:,1:initial_dims);
    
    %% HARMONY python addon
    if batcorr>0
        pyenv("ExecutionMode","OutOfProcess");
        disp('Batch correction with Harmony py...')
        hprj = prj;
        batchid = all_flags(1).fc_time;%  1 flag for the species !!
        usepylib = 1;
        [sout]=harmonypy(hprj,batchid,usepylib);
        prj = sout;
        terminate(pyenv)
    end
    %%
    cd(direct)
    disp('projecting...')
    init = prj(:,1:2)/std(prj(:,1))*1e-4;
    init = init-repmat(mean(init),length(init),1);
    if size(prj,1)>100000 % prj size (ncells, ndimensions) for big data 
        D = squareform(pdist(prj(2:2:end,:),'correlation'),'tomatrix');
    else % for smaller size data
        D = squareform(pdist(prj,'correlation'),'tomatrix');
    end
    [Dsort,~] = sort(D,'ascend');
    per_range = 500;%length(D);
    x = 1:per_range;%length(D);
    optk = zeros(length(D),1);
    for i=1:length(D)
        %     i;
        y = Dsort(1:per_range,i);
        x = x(:);
        y = y(:);
        a = atan((y(end)-y(1))/(x(end)-x(1)));
        xn = x*cos(a) + y*sin(a);
        yn = -x*sin(a) + y*cos(a);
        [~,imax] = max(yn);
        optk(i) = round((x(imax)));
    end
    disp('t-SNE...')
    perplexity = median(optk);
    options = statset('MaxIter',1000);
    mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
        'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',size(prj,1)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
    
    
    if savefig_flag==1
        savefig(gcf,'PCA.fig') % take property name
    end
    %
    % this is just the initial tsne, can be commented later
    clear Dsort
    %% intial tsne
    figure('Name','intial ts','NumberTitle','off');
    set(gcf,'color','w')
    dscatter(mapped_xy(:,1),mapped_xy(:,2));axis tight; axis off; axis equal;colorbar
    title(['#cells: ',num2str(length(mapped_xy)),'|%cells/all:',num2str(round(length(mapped_xy)/length(cellid_org))),' |max UMI: ',num2str(maximal_total_molecule),' |90th% UMI:',num2str(perct_90th_molecule)])
    
    %%% plot by sample
    % %     % sample_uni = {'8-1','10-1','18-1','19-1','23-1','23-3','41-1','41-2','42-1','45-1','45-2','45-3','46-1','46-2','46-3'};
    % %     colors = distinguishable_colors(length(sample_uni)+1);
    % %     figure('Name','all samples','NumberTitle','off');
    % %     set(gcf,'color','w','position',[1000,20,900,900]);
    % %     % [zha, pos] = tight_subplot(1, 3, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
    % %     % axes(ha(1))
    % %     for i=1:length(sample_uni)
    % %         s = scatter(mapped_xy(strcmpi(sample,sample_uni{i}),1),mapped_xy(strcmpi(sample,sample_uni{i}),2),10,colors(i,:),'filled'); hold on;
    % %     end
    % % %         axis tight
    % % %         axis equal
    % % %         axis off
    % %     legend(sample_uni)
    if savefig_flag==1
        savefig(gcf,'tsne_all.fig')
    end
    %% per sample
%     disp('samples tsne')
%     figure('Name','per sample ts','NumberTitle','off');
%     set(gcf,'color','w');
%     [px,~]=numSubplots(length(sample_uni));
%     [ha, pos] = tight_subplot(px(1),px(2), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
%     colors = distinguishable_colors(length(sample_uni)+1);
%     
%     for i=1:length(sample_uni)
%         axes(ha(i))
%         s=dscatter(mapped_xy(strcmpi(sample,sample_uni{i}),1),mapped_xy(strcmpi(sample,sample_uni{i}),2));%,10,colors(i,:),'filled');
%         %     alpha(s,0.4);
%         axis tight
%         axis equal
%         axis off
%         legend(sample_uni(i));
% % %         title([x_samples(i,findx),['tot:',num2str(length(mapped_xy(strcmpi(sample,sample_uni{i}),1))),'|',char(shower(1)),':', num2str(sum(snap25(strcmpi(sample,sample_uni{i}))>0))]],'FontSize',8)
%     end
%     if savefig_flag==1
%         savefig(gcf,'tsne_per_sample.fig')
%     end
%     
    %% all samples
    figure('Name','all samples ts','NumberTitle','off');
    set(gcf,'color','w');
    colors = distinguishable_colors(length(sample_uni));
    
    for i=1:length(sample_uni)
        s = scatter(mapped_xy(strcmpi(sample,sample_uni{i}),1),mapped_xy(strcmpi(sample,sample_uni{i}),2),10,colors(i,:),'filled'); hold on;
        %     alpha(s,0.4);
    end
    legend(sample_uni);
    axis tight
    axis equal
    axis off
    if savefig_flag==1
        savefig(gcf,'tsne_all_samples.fig')
    end
     %% all cellid (here is the cell types flags!!)
    figure('Name','all samples ts','NumberTitle','off');
    set(gcf,'color','w');
%     [uc,un,ui]=unique(cellid,'stable');
%     ucx=regexprep(uc,'_','-');

    colors = distinguishable_colors(length(uc));
    
    for i=1:length(uc)
        s = scatter(mapped_xy(strcmpi(cellid,uc{i}),1),mapped_xy(strcmpi(cellid,uc{i}),2),10,colors(i,:),'filled'); hold on;
        %     alpha(s,0.4);
    end
    legend(ucx);
    axis tight
    axis equal
    axis off
    if savefig_flag==1
        savefig(gcf,'tsne_all_samples.fig')
    end
    % legend({'8-1,Amy ct','10-1,Amy FC'...
    %     ,'18-1,Amy FC','19-1,Amy FC','23-1,Amy ct','23-3,Amy ct','41-1,Amy,24h','41-2,Amy,24h'....
    %     ,'42-1,Amy,8h','42-2,Amy,8h','45-1,Amy,24h','45-2,Amy,24h','45-3,Amy,ct','46-1,Amy,8h'....
    %     ,'46-2,Amy,8h','46-3,Amy,ct'})
    %% plot  per flags
    % plot by CONDITION
    disp('plot flags')
    if flags==1
        for ii=1:size(all_flags,2)
% % %             legflg=string(unique(are_samples(samples_idx,findx(ii))));
            fc_time_uni=unique(all_flags(ii).fc_time,'stable');
            colors = distinguishable_colors(length(fc_time_uni));
            [px,~]=numSubplots(length(fc_time_uni));
            figure('Name',['per flag',num2str(ii)],'NumberTitle','off');
            set(gcf,'color','w')
            
            [ha, pos] = tight_subplot(px(1),px(2), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
            for i=1:length(fc_time_uni)
                axes(ha(i))
                s = dscatter(mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),1),mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),2));%,10,colors(i,:),'filled');
                axis tight
                axis equal
                axis off
%                 legend('m');
                title(num2str(length(mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),1))))
            end
            if savefig_flag==1
% %                 savefig(gcf,['tsne_perflag_',char(flist(findx(ii)-2)),'.fig']) % take property name
                savefig(gcf,['tsne_per_flag_',num2str(ii),'.fig'])

            end
            %  % plot  all flags
            figure('Name',['flag ts',num2str(ii)],'NumberTitle','off');
            set(gcf,'color','w')
            
            for i=1:length(fc_time_uni)
                hold on
                s = scatter(mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),1),mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),2),10,colors(i,:),'filled'); hold on;
            end
            axis tight
            axis equal
            axis off
% % % %             title(flist(findx(ii)-2))
% legend('b','m') % take legend name
            legend(string(fc_time_uni),'Interpreter', 'none') % take legend name
% % % %             legend(legflg) % take legend name
            if savefig_flag==1
% %                 savefig(gcf,['tsne_all_flags_',char(flist(findx(ii)-2)),'.fig']) % take property name
savefig(gcf,['tsne_all_flags_',num2str(ii),'.fig'])
            end % silute to here
        end
        
    end
    
    % % % % %
    % Z = linkage(mapped_xy,'ward','euclidean');
    % c = cluster(Z,'Maxclust',120);
    % idx = c;
    % clustering with dbscan
    %% NO ! clustring !
    %% clustring method
    if 1 
    disp('clustring...')
    clmethod
    if clmethod=="K-means"
        markfig=figure('Name','preview','NumberTitle','off');
        set(gcf,'color','w')
        
        while 1
            
            idx = kmeans(prj,eps_prc);
            
            gscatter(mapped_xy(:,1),mapped_xy(:,2),idx);
            axis off
            axis tight
            axis equal
            hold off
            if stops==1
                redo=input('redo(y)? ','s');
                if redo=='y'
                    K=eps_prc
                    eps_prc=input('K: ');
                else
                    break;
                end
            else
                break;
            end
        end
        
    elseif clmethod=="KNN-louvain"
        markfig=figure('Name','preview','NumberTitle','off');
        set(gcf,'color','w')
        
        
        adj = knn_mat(prj,perplexity,2);
        addpath('/data/matlab_functions/louvain_scgeatoolbox/')
        [com, Q, comty] = louvain(adj);
        
        partition=length(comty.COM)
        while 1
            try
                idx = comty.COM{eps_prc}';
                gscatter(mapped_xy(:,1),mapped_xy(:,2),idx);
                axis off
                axis tight
                axis equal
                hold off
            catch
                disp('K chosen is bigger than partition')
            end
            if stops==1
                redo=input('redo(y)? ','s');
                if redo=='y'
                    K= eps_prc
                    eps_prc=input('K: ');
                else
                    break;
                end
            else
                break;
            end
            
        end
        
    elseif clmethod=="Linkage"
        markfig=figure('Name','preview','NumberTitle','off');
        set(gcf,'color','w')
        
        
        zx = linkage(prj,'ward','correlation');
        Dx = pdist(prj,'correlation');
        leaforder = optimalleaforder(zx,Dx);
        while 1
            
            idx = cluster(zx,'MaxClust',eps_prc);
            
            gscatter(mapped_xy(:,1),mapped_xy(:,2),idx);
            axis off
            axis tight
            axis equal
            hold off
            redo=input('redo(y)? ','s');
            if stops==1
                if redo=='y'
                    K=eps_prc
                    eps_prc=input('K: ');
                else
                    break;
                end
            else
                break;
            end
        end
        
    else %dbscan and knn-dbscan
        %% markers 1 + dbscan clustring
% %         data_orig_allx=data_orig_all;% non 
% % mapped_xyz=mapped_xy;
%         data_orig_all=data_orig_allx(:,batchid~=0);
%         mapped_xy=mapped_xyz(batchid~=0,:);
        disp("Scatter genes on 2D embedding")
        %     markfig=figure('Name','preview','NumberTitle','off');
        %     set(gcf,'color','w')
        
        %     vClust and cclsut and catclust
        while 1
            markfig=figure('Name','preview','NumberTitle','off');
            set(gcf,'color','w')
            [p2,~]=numSubplots(5);
            compumi=zeros(size(mapped_xy,1),(length(ex_clust)+2));
            %     set(0,'CurrentFigure',markfig)
            [ha, pos] = tight_subplot(p2(1), p2(2), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
            
            ax_link=[];
            for i=1:length(ex_clust)
                genePlot = ex_clust{i};
                markergene = (data_orig_all(strcmpi(geneid_all(:,1),genePlot),:));
                %         markergenex=markergene;
                %         markergenex(markergene<5)=0;% filter for compare
                try
                compumi(:,i)=markergene;%markergenex';
                % catch
                %         compumi(:,i)=markergene(1,:)';%markergenex';
                
                inpos = markergene>0;
                tmpthlow = prctile(markergene(markergene>0),1);
                tmpthhigh = prctile(markergene(markergene>0),99);
                if tmpthlow==tmpthhigh
                    tmpthlow = 0;
                end
                markergene(markergene>tmpthhigh) = tmpthhigh;
                markergene(markergene<tmpthlow) = tmpthlow;
                c_rgb = [1,0,0];rand([1,3]);
                %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
                %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
                markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
                axes(ha(i));
                scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
                scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
                set(gca,'xlim',[-150,150],'ylim',[-150,150])
                ax_link=[ax_link,ha(i)];% I linked the axeses
                
                
                title([genePlot,':',num2str(length(mapped_xy(inpos,1)))]);
                axis tight
                axis equal
                axis off
                end
            end
            if groupy==0
                grp=0;
                [~,loc] = ismember(list2,geneid_all);
                list2(loc==0)=[];
                [~,loc] = ismember(list2,geneid_all);
                nonneuro = sum(data_orig_all(loc,:));
                axes(ha(4));
                
                markergene = nonneuro;
                %        markergenex=markergene;
                %         markergenex(markergene<5)=0;% filter  value to compare <5 zeroing
                compumi(:,4)=markergene;% markergenex';
                inpos = markergene>0;
                tmpthlow = prctile(markergene(markergene>0),1);
                tmpthhigh = prctile(markergene(markergene>0),90);
                markergene(markergene>tmpthhigh) = tmpthhigh;
                markergene(markergene<tmpthlow) = tmpthlow;
                c_rgb = [1,0,0];rand([1,3]);
                markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
                scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
                scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
                axis off; axis tight; axis equal
                
                ax_link=[ax_link,ha(4)];
                linkaxes(ax_link,'xy');
                title(['Group:',num2str(length(mapped_xy(inpos,1)))]);
                
            end
            
            axes(ha(5));
            for icomp=1:length(compumi)
                if (compumi(icomp,4)> fct*compumi(icomp,1)) && (compumi(icomp,4)> fct*compumi(icomp,2)) && (compumi(icomp,4)> fct*compumi(icomp,3))
                    compumi(icomp,5)=4;% group
                elseif (compumi(icomp,1)> fct*compumi(icomp,2)) && (compumi(icomp,1)> fct*compumi(icomp,3))
                    compumi(icomp,5)=1; % g1
                elseif (compumi(icomp,2)> fct*compumi(icomp,1)) && (compumi(icomp,2)> fct*compumi(icomp,3))
                    compumi(icomp,5)=2; % g2
                elseif (compumi(icomp,3)> fct*compumi(icomp,2)) && (compumi(icomp,3)> fct*compumi(icomp,1))
                    compumi(icomp,5)=3; % g3
                    
                    %          else
                    % zero
                end
            end
            gscatter(mapped_xy(:,1),mapped_xy(:,2),compumi(:,5));
            title(['NA:',num2str(sum(compumi(:,5)==0)),'|',char(ex_clust(1)),':',num2str(sum(compumi(:,5)==1)),'|',char(ex_clust(2)),':',num2str(sum(compumi(:,5)==2)),'|',char(ex_clust(3)),':',num2str(sum(compumi(:,5)==3)),'|Group:',num2str(sum(compumi(:,5)==4))])
            legx=legend(["NA",ex_clust,"Group"],'AutoUpdate','off');% delete legend by not adding leg(end)=[];
            ax_link=[ax_link,ha(5)];
            axis off; axis tight; axis equal
            linkaxes(ax_link,'xy');
            
            % % % dbscan
            set(gca,'FontSize',8) % Creates an axes and sets its FontSize to 8
            disp("DBSCAN Clustring...")
            tStartx=tic;
            % MinPts = 15
            % eps_prc = 70
            axes(ha(6))
            axis off; axis tight; axis equal
            text(0,0,'Clustring & Filtering please wait...')
            drawnow;
            % here is the actual clustring :
            try
                if eps_prc==eps_prcx && MinPts==MinPtsx
                    idx=idxx;
                else
                    %           [idx, isnoise] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);
                    epsilon= dbscan_epsprc_mipts_ed(mapped_xy,eps_prc,MinPts);
                    disp('...')
                    tic;
                    [idx, ~] = dbscan(Dx,epsilon,MinPts,'Distance','precomputed');
                    idx(idx==-1)=0;
                    toc
                    idxx=idx;eps_prcx=eps_prc ; MinPtsx=MinPts;
                end
            catch
                %             [idx, isnoise] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);
                epsilon= dbscan_epsprc_mipts_ed(mapped_xy,eps_prc,MinPts);
                disp('...')
                tic;
                Dx = pdist(mapped_xy);
                [idx, ~] = dbscan(Dx,epsilon,MinPts,'Distance','precomputed');
                idx(idx==-1)=0;
                toc
                idxx=idx;eps_prcx=eps_prc ; MinPtsx=MinPts;
            end
            if filtery~=0
                disp('Filtering...')
                
                % filtery is range of KNN
                adj=knn_mat(mapped_xy,abs(filtery),1);
                itemp=zeros(size(idx));
                for i=1:max(idx)
                    ii=find(idx==i);
                    if length(ii)<=abs(filtery)
                        i
                        for ix=1:length(adj)
                            
                            if filtery>0 % if positive filtery make the points as the closest point;
                                in =adj(ix,:)>0;
                                maj=mode(idx(in));
                                itemp(ix)=maj;
                            else% if negative filtery make the points as 0;
                                itemp(ix)=0;
                            end
                        end
                        idx(idx==i)=itemp(idx==i);% do it only for idx==0 if general filtring  do it idx==itemp
                    end
                end
                % check for gaps iin the clusters after filtering e.g. 13 15 21
                % 22 25....
                unidx=unique(idx);
                gapy=diff(unidx);gapy=[1;gapy];
                for iz=1:length(unidx)
                    if gapy(iz)>1
                        idx(idx==unidx(iz))=iz;
                    end
                end
            end
            
            
            
            colors = distinguishable_colors(length(unique(idx)));
            
            %             zcolor=0;
            for i=unique(idx)'
                if i>0
                    %                 txt = [num2str(i),':',num2str(sum(idx==i))]; % legend in a loop
                    %                zcolor=zcolor+1;
                    ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i,:),'markersize',5); hold on; %% 'DisplayName',txt)
                    xc=mapped_xy(ii,1);
                    yc=mapped_xy(ii,2);
                    scbound = boundary(xc,yc,0.1);
                    plot(xc(scbound),yc(scbound),'--','LineWidth',2,'Color','k');
                    hold on
                    ht = text(median(mapped_xy(ii,1)),median(mapped_xy(ii,2)),num2str(sum(idx==i)));
                    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',5);
                elseif i==0
                    ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'+','color',[0.5 0.5 0.5],'markersize',2); hold on;
                    
                end
            end
            ax_link=[ax_link,ha(6)];
            linkaxes(ax_link,'xy');
            axes(ha(6))
            axis off; axis tight; axis equal
            set(gca,'FontSize',8) % Creates an axes and sets its FontSize to 18
            title(['prp:',num2str(perplexity),'|min:',num2str(MinPts),'|eps:',num2str(eps_prc),'|C:',num2str(length(unique((idx)))),'|in/out:',num2str(sum(idx>=0)),'/',num2str(sum(idx==0)), ' ,%out:',num2str(round(100*sum(idx==0)/length(compumi))),'%',' ,Filter:',num2str(filtery)]);
            % % %
            axes(ha(5))
            for i=1:max(idx)
                hold on
                ii=find(idx==i);
                xc=mapped_xy(ii,1);
                yc=mapped_xy(ii,2);
                scbound = boundary(xc,yc,0.1);
                plot(xc(scbound),yc(scbound),'--','LineWidth',2,'Color','k');% show border line
                inx = idx==i;% show textbox
                ht = text(median(mapped_xy(inx,1)),median(mapped_xy(inx,2)),num2str(mode(compumi((idx==i),5))));
                set(ht,'BackgroundColor',0.2*[1,1,1],'fontsize',8);
                set(ht,'Color','w')
            end
            % end timer
            tEndx = toc(tStartx);
            fprintf('%d minutes and %f seconds\n', floor(tEndx/60), rem(tEndx,60));
            
            
            
            
            if stops==1
                redo = input('redo(y)?','s');
                if redo=="y"
                    eps_prc
                    eps_prc = input('epsilon:');
                    MinPts
                    MinPts = input('MinPts:');
                    fct
                    fct=input('compare factor:');
                    filtery
                    filtery=input('Filter:');
                else
                    break;
                end
            else
                break;
            end
            
        end
        
        % %     else % vclust
        % %         disp("DBSCAN Clustring...please wait...")
        % %         while 1
        % %
        % %             tStartx=tic;
        % %             [idx, isnoise] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);
        % %             idxuni = unique(idx);
        % %             colors = distinguishable_colors(length(unique(idx)));
        % %
        % %
        % %
        % %             for i=unique(idx)'
        % %                 if i>0
        % %                     %                 txt = [num2str(i),':',num2str(sum(idx==i))]; % legend in a loop
        % %                     ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i,:),'markersize',5); hold on; %% 'DisplayName',txt)
        % %
        % %                 elseif i==0
        % %                     ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'+','color',[0.5 0.5 0.5],'markersize',2); hold on;
        % %                 end
        % %             end
        % %             %         ax_link=[ax_link,ha(6)];
        % %             axis off; axis tight; axis equal;
        % %             %                 linkaxes(ax_link,'xy');
        % %             title(['prp:',num2str(perplexity),'|min:',num2str(MinPts),'|eps:',num2str(eps_prc),'|C:',num2str(max(idx)),'|in/out:',num2str(sum(idx>=0)),'/',num2str(sum(idx==0)), ' ,%out:',num2str(round(100*sum(idx==0)/length(compumi))),'%']);
        % %             % % %
        % %             tEndx = toc(tStartx);
        % %             fprintf('%d minutes and %f seconds\n', floor(tEndx/60), rem(tEndx,60));
        % %             if stops==1
        % %                 redo=input('redo(y) ?','s');
        % %                 if redo=='y'
        % %                     eps_prc = input('epsilon:');
        % %                     MinPts = input('MinPts:');
        % %                 else
        % %                     break;
        % %                 end
        % %             else
        % %                 break;
        % %             end
        % %
        % %         end
        % %     end
    end
end
    
    idx=ui;
    idx_zero=sum(idx==0)


    %% final  cluster  figure
    colors=distinguishable_colors(max(idx));
    T_cells_all=idx;
    fdbs=figure('Name',clmethod,'NumberTitle','off');
    set(gcf,'color','w')
    for i=unique(idx)'
        if i>0
            %                 txt = [num2str(i),':',num2str(sum(idx==i))]; % legend in a loop
            ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i,:),'markersize',5); hold on; %% 'DisplayName',txt)
            
        elseif i==0
            ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'+','color',[0.5 0.5 0.5],'markersize',2); hold on;
        end
    end
    %  % %           legend show
    for i=1:max(idx)
        in = idx==i;
        ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),num2str(sum(idx==i)));
        set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',8);
    end
    % for i=idxuni'
    %     if i>=0
    %         in = idx==i;
    %         ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),num2str(i));
    %         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
    %     end
    % end
    axis tight;
    axis equal
    axis off
    title(['peplexity=',num2str(perplexity),', MinPts=',num2str(MinPts),', epsprc/K=',num2str(eps_prc),',#C=',num2str(max(idx)),',#in=',num2str(sum(idx>=0)),',#out=',num2str(sum(idx==0)),...
        ' ,%out=',num2str(round(100*sum(idx==0)/length(idx))),'% ,Filter:',num2str(filtery),' , Clutser counts:'],'fontsize',8);
    
    
    if savefig_flag==1
        savefig(gcf,[char(clmethod),'.fig'])
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% sort the data by the clusters and remove outliers
    disp("Sort the data by the clusters and remove outliers")
    
    
    [idx,xi] = sort(idx);
    xi(idx==0) = [];
    idx(idx==0) = [];
    idxuni = unique(idx);
    data_sorted_all = data(:,xi);
    data_orig_all_sorted = data_orig_all(:,xi);
    cellid_sorted = cellid((xi));
    sample_sorted = sample((xi));
    s_xy_org=mapped_xy;% orignal t-sne xy
    mapped_xy = mapped_xy(xi,:);% filtered t-sne xy
    for ii=1:size(all_flags,2)
        all_flags(ii).fc_time_sorted = all_flags(ii).fc_time(xi);
    end
    % % amy_flag_sorted = amy_flag(xi);
    % % piri_flag_sorted = piri_flag(xi);
    prj_sorted = prj(xi,:);
    
    %%
    disp('inner t-SNE')
    no_dims = 1;
    initial_dims = 10;
    perplexity = 5;
    % epsilon = 100;
    dist_flag = 2;
    theta = 0.5;
    rand_seed = 13;
    data_tsne = cent_norm(log2(data_sorted_all+1));
    xi = [1:length(idx)];
    for i=1:length(idxuni)
        i
        ind = find(idx==i);
        %     if length(ind)>20000
        try
            tmp1d = tsne((data_tsne(:,ind))','Algorithm','barneshut','Distance','correlation','NumDimensions',no_dims,'NumPCAComponents',initial_dims,.....
                'Perplexity',perplexity,'Standardize',true,'LearnRate',100,'Theta',theta,'Verbose',1,'Options',options,'Exaggeration',20);
            [~,xitmp] = sort(tmp1d);
            xi(ind) = xi(ind((xitmp)));
            %     elseif length(ind)<20000%% && length(ind)>100
            % %         tmp1d = fast_tsne((data_tsne(:,ind))', no_dims, initial_dims, perplexity,theta, rand_seed);
            % %         [~,xitmp] = sort(tmp1d);
            % %         xi(ind) = xi(ind((xitmp)));
            %     else % below %
            %     end
        catch
        end
    end
    
    data_sorted_all = data_sorted_all(:,xi);
    data_orig_all_sorted = data_orig_all_sorted(:,xi);
    cellid_sorted = cellid_sorted((xi));
    sample_sorted = sample_sorted((xi));
    mapped_xy = mapped_xy(xi,:);
    for ii=1:size(all_flags,2)
        all_flags(ii).fc_time_sorted = all_flags(ii).fc_time_sorted(xi);
    end
    % % % fc_time_sorted = fc_time_sorted(xi);
    
    prj_sorted = prj_sorted(xi,:);
    
    meangr_mat = zeros(length(moldata(:,1)),length(idxuni));
    clust_cent = zeros(length(idxuni),2);
    for jjj=1:length(idxuni)
        jjj
        meangr_mat(:,jjj) = mean(log2(data_sorted_all(:,idx==idxuni(jjj))+1),2);
        clust_cent(jjj,:) = [median(mapped_xy(idx==idxuni(jjj),1)),median(mapped_xy(idx==idxuni(jjj),2))];
    end
    meangr_mat1 = meangr_mat;
    % meangr_mat1(loc,:) = [];
    %% dendogram
    
    %     meangr_mat1 = cent_norm(meangr_mat1(:,leaforder));
    if max(idxuni)>9
        [prj,m,D,V,Q] = pca_wis(meangr_mat1',initial_dims);
    else
        [prj,m,D,V,Q] = pca_wis(meangr_mat1',max(idxuni));
    end
    % [prj,m,D,V,Q] = pca_wis(meangr_mat1);
    % prj=prj(1:max(idxuni),:);
    Zpca = linkage(prj,'ward','correlation');
    Dpca = pdist(prj,'correlation');
    leaforder_pca = optimalleaforder(Zpca,Dpca);
    disp("Dendogram tree")
    figure('Name','dendogram','NumberTitle','off');
    set(gcf,'color','w')
    axes('position',[0.03,0.03,0.3,0.93])
    % hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
    hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
    axis off
    set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
    axes('position',[0.35,0.03,0.63,0.93])
    x=squareform(Dpca); imagesc(x(leaforder_pca,leaforder_pca));
    colormap('summer')
    set(gca,'ytick',[1:length(leaforder_pca)],'xtick',[],'fontsize',8,'ydir','normal')
    
    if savefig_flag==1
        savefig(gcf,'tree.fig')
    end
    leaforder = leaforder_pca;
    leaforder=1:length(ucx);
    Zpca_post = linkage(prj(leaforder_pca,:),'ward','correlation');
    T_cells_tmp_new = zeros(size(idx));
    for i=1:length(leaforder)
        T_cells_tmp_new(idx==idxuni(leaforder(i))) = i;
    end
    idxuni_new = unique(T_cells_tmp_new);
    [~,xi] = sort(T_cells_tmp_new);
    T_cells_tmp_new = T_cells_tmp_new(xi);
    data_sorted_all = data_sorted_all(:,xi);
    data_orig_all_sorted = data_orig_all_sorted(:,xi);
    cellid_sorted = cellid_sorted(xi);
    % % tissue = tissue(xi);
    mapped_xy = mapped_xy(xi,:);
    for ii=1:size(all_flags,2)
        all_flags(ii).fc_time_sorted = all_flags(ii).fc_time_sorted(xi);
    end
    % % % fc_time_sorted = fc_time_sorted(xi);
    cells_bor_2 = find(diff(T_cells_tmp_new)>0)+1;
    sample_sorted = sample_sorted(xi);
    prj_sorted = prj_sorted(xi,:);
    T_cells_tmp = T_cells_tmp_new;
    T_cells_tmp_uni = unique(T_cells_tmp);
    % % % %
    idx = T_cells_tmp_new;
    idxuni = idxuni_new;
    % save junctions
    try  %????
        [table1, table2] = dendrogram_split_markers(cellfun(@(x) num2str(x), m2c(T_cells_tmp_uni),'UniformOutput',0).....
            ,cellfun(@(x) num2str(x), m2c(T_cells_tmp),'UniformOutput',0),Zpca_post,data_sorted_all,geneid);
        saveCellFile(table1,'junctions.txt');
        saveCellFile(table2,'junctions.txt');
    end
    %% tsne with annoated clusters
    %      if clmethod=='DBSCAN'
%      u_c=regexprep(u_c,'_','-');
%     u_c=regexprep(u_c,'b-','');
%     u_c=u_c(leaforder);
    %     colors = distinguishable_colors(length(unique(idx))+1);
    disp('filtered')
    figure('Name','filtered','NumberTitle','off');
    set(gcf,'color','w')
%     subplot(1,2,1)
    
    for i=unique(idx)'
        if i==0
            ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',[0.5,0.5,0.5],'markersize',3); hold on;
        else
            ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i,:),'markersize',5); hold on;
            groupval(i)=sum(ii);
        end
    end
    for i=idxuni'
        in = idx==i;
        ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),num2str(i));%u_c(i) or num2str(i)
%         ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),u_c(i));%u_c(i) or num2str(i)
        set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',8)
    end
    axis tight;
    axis equal
    axis off
    title(['MinPts=',num2str(MinPts),', epsprc=',num2str(eps_prc),',#C=',num2str(max(idx)),' ,#tot=',num2str(length(mapped_xy)),' , Clutser numbers:']);
    
%     subplot(1,2,2)
%     axis tight;
%     axis equal
%     axis off
%     donut(groupval,{},colors);
%     legend(string(round(100*groupval/sum(groupval))));
%     text(0,0,"%")
    if savefig_flag==1
        savefig(gcf,'filtered.fig')
    end
    %     end
    
    %% heatmap
    uci=unique(cellid_sorted,'stable');
    [ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp,data_sorted_all,1);% top_genes by me
    % % % % % % % % %
    datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
    
    datamarkers_cn = cent_norm(log2(datamarkers+1));
    % gr_tmp_mark = gr_tmp_mark(xi);
    gr_tmp_mark = geneid(ind_gr_tmp_mark);
    
    figure('Name','HeatMap','NumberTitle','off');
    set(gcf,'color','w')
    ax1 = axes('position',[0.1,0.1,0.88,0.73]);
    imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
    hold on;
    linewid =0.5;
    bor_color = 'grey11';%'green1';%
    for jj=1:length(cells_bor)
        plot(cells_bor(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
    end

%     set(gca,'xtick',gr_center,'xticklabel',[1:length(gr_center)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 8)
    set(gca,'xtick',gr_center,'xticklabel',uci,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 8)
    
    colormap('summer');
    freezeColors(gca);
    
    % sample_uni = {'8-1','10-1','18-1','19-1','23-1','23-3','41-1','41-2','42-1','45-1','45-2','45-3','46-1','46-2','46-3'};
    samples_num = false(length(sample_uni),length(sample_sorted));
    for i=1:length(sample_uni)
        samples_num(i, strcmpi(sample_sorted,sample_uni{i})) = true;
    end
    %sample axes
    ax5 = axes('position',[0.1,0.83,0.88,0.135]);
    imagesc(~samples_num); hold on;
    grid on
    colormap('gray');
    freezeColors(gca);
    set(gca,'xtick',[],'ytick',[1:length(sample_uni)],'yticklabel',sample_uni,'fontsize',6);
    %faxes
    f_link=[];
    if batcorr==2 || batcorr==-1
        starter=2;
        posx=2;
    else
        starter=1;
        posx=1;
    end
    for ii=starter:size(all_flags,2)
        % fc_uni=unique(all_flags(ii).fc_time_sorted)
        fl_name=string(unique(are_samples(samples_idx,findx(ii))));
        fl_name(cellfun('isempty',fl_name)) = []; % remove empty cells from mm
        for iii=1:length(fl_name)
            ax2 = axes('position',[0.1,0.75+((iii-1)*0.012)+((ii-posx)*length(fl_name)*0.012),0.88,0.012]);
            imagesc(~(all_flags(ii).fc_time_sorted'==iii)); hold on;
            colormap('gray');
            freezeColors(gca);
            set(gca,'xtick',[],'ytick',[1],'yticklabel',fl_name(iii),'fontsize',7);
            ax2.TickDirMode = 'auto';
            f_link=[f_link,ax2];
        end
    end
    
    gads=[];
    g_link=[];
    %gaxes
    for gena=1:length(ex_clust)
        ax6 = axes('position',[0.1,0.965+((gena-1)*0.01),0.88,0.01]);
        gad2 = data_orig_all_sorted(strcmpi(geneid_all(:,1),ex_clust(gena)),:);%cluster genes
        imagesc(~gad2); hold on;
        axes(ax6)
        colormap('gray');
        freezeColors(ax6);
        set(ax6,'xtick',[],'ytick',[1],'yticklabel',ex_clust(gena),'fontsize',5);
        g_link=[g_link,ax6];% I linked the axeses
        gads=[gads;gad2];
    end
    linkaxes([ax1,f_link,ax5,g_link],'x')
    if savefig_flag==1
        savefig(gcf,'markertable.fig')
    end
    % % if savefig_pdf==1
    % %     eval('export_fig markertable.pdf');
    % % end
    % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % %
    % list = {'Snap25','Stmn2','Gad2','Slc32a1','Slc17a7','Slc17a6','Sst','Sim1','Foxj1','Pdgfra','Mog','C1qc','Flt1','Cldn5','Aqp4','Plp1'};
    %% groupmap
    
    
 
    datamarkers = meangr_mat(ind_gr_tmp_mark,leaforder_pca);
    datamarkers_cn = datamarkers;%cent_norm(datamarkers);(lo,:)
    
    disp("GroupMap")
    figure('Name','GroupMap','NumberTitle','off');
    set(gcf,'color','w')
    ax1 = axes('position',[0.1,0.14,0.88,0.65]);
    imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
    hold on;
    % % clusteruni={1:size(meangr_mat,2)};
    set(gca,'xtick',[1:length(gr_center)],'xticklabel',[1:length(gr_center)],'XTickLabelRotation',45,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)%,'xticklabel',regexprep(clusteruni(leaforder_pca),'_','-')....
    %     ,'XTickLabelRotation',45,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)
    colormap('summer');
    freezeColors(gca);
    
    h2 = axes('position',[0.1,0.81,0.88,0.15]);
    % hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
    hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','top');
    axis off
    set(gca,'xlim',[0.5,length(leaforder_pca)+0.5])
    linkaxes([ax1,h2],'x')
    if savefig_flag==1
        savefig(gcf,'grouptable.fig')
    end
    %% search gene summary
    t=cd;
    cd('/bigdata');
    gene_summary=readtable('bigdata/gene_summary_all.csv');
    cd(t);
    gene_name=table2array(gene_summary(:,1));
    gsummary={};
    for gs=1:length(gr_tmp_mark)
        igs=find(gene_name==string(gr_tmp_mark(gs)))
        if ~isempty(igs)
        gsummary(gs,[1 2 3])=table2array(gene_summary(igs(1),:));
        end
    end
    gsummary=table(gsummary)
    writetable(gsummary,'gsummary.csv')
    %% search cell type names 
    suggested_ID=run.alona(data_orig_all_sorted,string(geneid_all),T_cells_tmp,'species','mouse');
    sid= table2array(suggested_ID);
    sidx=sid(1,1:2:end);
    sidv=double(sid(1,2:2:end));
    NeuroID=(sidv>40) & (sidx=="Neurons") |(sidx=="Interneurons");
    TCI=1:length(sidx);
    sum(NeuroID)
    num_cells=sum(sum(T_cells_tmp==TCI(NeuroID),2))
    sidx(NeuroID)="EXCLUDE";
    sidx=sidx';

    %% markers 2
%     disp("Scatter genes on t-SNE")
%     figure('Name','genes & ts','NumberTitle','off');
%     set(gcf,'color','w')
%     [p2,~]=numSubplots(length(list));
%     [ha, pos] = tight_subplot(p2(1), p2(2), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
%     link_axes=[];
%     for i=1:length(list)
%         genePlot = list{i};
%         markergene = (data_orig_all_sorted(strcmpi(geneid_all(:,1),genePlot),:));
%         if ~isempty(markergene)
%             %             markergene = mean(markergene);
%             inpos = markergene>0;
%             tmpthlow = prctile(markergene(markergene>0),1);
%             tmpthhigh = prctile(markergene(markergene>0),99);
%             if tmpthlow==tmpthhigh
%                 tmpthlow = 0;
%             end
%             markergene(markergene>tmpthhigh) = tmpthhigh;
%             markergene(markergene<tmpthlow) = tmpthlow;
%             c_rgb = [1,0,0];rand([1,3]);
%             %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%             %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
%             markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
%                 interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
%                 ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
%             axes(ha(i));
%             scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
%             scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
%             set(gca,'xlim',[-150,150],'ylim',[-150,150])
%             title([genePlot,':',num2str(length(mapped_xy(inpos,1)))]);
%             axis tight
%             axis equal
%             axis off
%         end
%     end
%     
    
    %% groups
    if clusty~="vClust"
        % % % % show gad2 and slc17a7 slc17a6
        % %  list2 = ({'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
        % %     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'});
        if groupy==0
            
            grp=0;
            bar_ar_sorted=[];
            [~,loc] = ismember(list2,geneid_all);
            list2(loc==0)=[];
            [~,loc] = ismember(list2,geneid_all);
            nonneuro = sum(data_orig_all_sorted(loc,:));
            gabaglut = zeros(size(T_cells_tmp_uni));
            gabaglut_sc = zeros(size(T_cells_tmp));
            
            for jjj=1:length(T_cells_tmp_uni)
                jjj   % cluster number
                tmp=[];
                for ggg=1:size(gads,1)
                    gad2=gads(ggg,:);
                    tmp=[tmp;mean(gad2(:,T_cells_tmp==jjj)>0)];
                end
                tmp = [tmp;mean(nonneuro(:,T_cells_tmp==jjj)>0)];
                
                tmpsort = [tmp(1),max(tmp(2:3)),tmp(4)];
                tmpsort= sort(tmpsort,'descend');
                if tmpsort(1)>fct*tmpsort(2)% compare factor
                    [~,gabaglut(jjj)] = max(tmp);
                else
                    gabaglut(jjj) = 4;
                end
                gabaglut_sc(T_cells_tmp==jjj) = gabaglut(jjj);
            end
            figure('Name','Groups & ts','NumberTitle','off');
            set(gcf,'color','w')
            subplot(1,2,1);
            colorx=['b','g','r','k','m']; %distinguishable_colors(5);%length(unique(gabaglut_sc)));
            groupvalx=zeros(1,5);
            %% count gluts and gads together
            if clusty== "cClust"
                gabaglut_sc(gabaglut_sc==2)=1;
                gabaglut_sc(gabaglut_sc==3)=1;
                %         gabaglut_sc(gabaglut_sc==5)=1;
            end
            %%
            for idx=1:5
                ii=find(gabaglut_sc==idx);
                h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colorx(idx),'markersize',3); hold on;
                groupvalx(idx)=length(ii);
            end
            title([' #:-',char(ex_clust(1)),':',num2str(groupvalx(1)),'|',char(ex_clust(2)),':',num2str(groupvalx(2)),'|',char(ex_clust(3)),':',num2str(groupvalx(3)),'|Groups:',num2str(groupvalx(4)),'|Doublets:',num2str(groupvalx(5))])
            %             title(['1:',num2str(groupval(1)),'|2:',num2str(groupval(2)),'|3:',num2str(groupval(3)),'|4:',num2str(groupval(4)),'|5:',num2str(groupval(5))])
            axis tight
            axis equal
            axis off
            % title('GABA/Glut1/Glut2')
            leg = {num2str(ex_clust(1)),num2str(ex_clust(2)),num2str(ex_clust(3)),'Group','doublets'};
            legend(leg(unique(gabaglut_sc)))
            subplot(1,2,2);
            donut(groupvalx,{},{[0 0 1],[0 1 0],[1 0 0],[0 0 0],[1 0 1]});
            axis tight
            axis equal
            axis off
            legend(string(round(100*groupvalx/sum(groupvalx))));
            text(0,0,"%")
            %         title([' %:-',char(ex_clust(1)),':',num2str(100*groupvalx(1)/sum(groupvalx)),'|',char(ex_clust(2)),':',num2str(100*groupvalx(2)/sum(groupvalx)),'|',char(ex_clust(3)),':',num2str(100*groupvalx(3)/sum(groupvalx)),'|Groups:',num2str(100*groupvalx(4)/sum(groupvalx)),'|Doublets:',num2str(100*groupvalx(5)/sum(groupvalx))])
            if savefig_flag==1
                savefig(gcf,'category.fig')
            end
            %             table1 = [cellid_sorted,m2c(gabaglut_sc)];
            %             saveCellFile(table1,'exclust.txt');
            % % legend('GABA','Glut1','Glut2','non-neurons','doublets')
            
        else
            %% for group analysis
            for grp=1:length(showerg)
                shower1=strsplit(showerg(grp),',');
                list2=convertStringsToChars(shower1);
                if upgn==1
                    list2=upper(list2)
                end
                [~,loc] = ismember(list2,geneid_all);
                list2(loc==0)=[];
                [~,loc] = ismember(list2,geneid_all);
                nonneuro = sum(data_orig_all_sorted(loc,:));
                gene_count=[gene_count,markergene];
                gabaglut = zeros(size(T_cells_tmp_uni));
                gabaglut_sc = zeros(size(T_cells_tmp));
                for jjj=1:length(T_cells_tmp_uni)
                    jjj;   % cluster number
                    % %             for ggg=1:size(gads,1)
                    % %                 gad2=gads(ggg,:);
                    % %                 tmp=[tmp;mean(gad2(:,T_cells_tmp==jjj)>0)];
                    % %             end
                    tmp = mean(nonneuro(:,T_cells_tmp==jjj)>0);
                    %                     tmpsort = [tmp(1),max(tmp(2:3)),tmp(4)];
                    %                     tmpsort= sort(tmpsort,'descend');
                    %                     if tmpsort>0
                    [~,gabaglut(jjj)] = max(tmp);
                    %                     else
                    %                         gabaglut(jjj) = 2;
                    %                     end
                    gabaglut_sc(T_cells_tmp==jjj) = gabaglut(jjj);
                end
                % %                 figure('Name','Groups & ts','NumberTitle','off');
                % %                 set(gcf,'color','w')
                % %                 colors = distinguishable_colors(5);%length(unique(gabaglut_sc)));
                %                 for idx=unique(gabaglut_sc)'
                %                     ii=find(gabaglut_sc==idx);
                %                     h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(idx,:),'markersize',3); hold on;
                %                     groupval(idx)=length(ii);
                %                 end
                cid = cellfun(@(x,y) [x,'_',y], cellid_sorted, sample_sorted,'UniformOutput',0);
                h=figure('Name',['Grp & ts ',num2str(grp)],'NumberTitle','off');
                set(gcf,'color','w')
                markergene = nonneuro;
                inpos = markergene>0;
                tmpthlow = prctile(markergene(markergene>0),1);
                tmpthhigh = prctile(markergene(markergene>0),90);
                markergene(markergene>tmpthhigh) = tmpthhigh;
                markergene(markergene<tmpthlow) = tmpthlow;
                c_rgb = [1,0,0];rand([1,3]);
                markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
                scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
                scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
                axis tight
                axis equal
                axis off
                newStr = join(list2,"-");
                xlabel((newStr))
                title([num2str(grp),'G: ',num2str(length(mapped_xy(inpos)))])
                %                 axis tight
                %                 axis equal
                %                 axis off
                % title('GABA/Glut1/Glut2')
                
                h = gca;                  %fixed relative to jonas's suggestion
                h.XAxis.Label.Visible='on';
                
                
                if savefig_flag==1
                    savefig(gcf,[num2str(grp),'_group.fig'])
                end
                % % % % % % % % % % % %         cid = cellfun(@(x,y) [x,'_',y], cellid_sorted, sample_sorted,'UniformOutput',0);
                % % % % % % % % % % % %         figure('Name','Groups & ts','NumberTitle','off');
                % % % % % % % % % % % %         set(gcf,'color','w')
                % % % % % % % % % % % %
                % % % % % % % % % % % %         markergene = nonneuro;
                % % % % % % % % % % % %         inpos = markergene>0;
                % % % % % % % % % % % %         tmpthlow = prctile(markergene(markergene>0),1);
                % % % % % % % % % % % %         tmpthhigh = prctile(markergene(markergene>0),90);
                % % % % % % % % % % % %         markergene(markergene>tmpthhigh) = tmpthhigh;
                % % % % % % % % % % % %         markergene(markergene<tmpthlow) = tmpthlow;
                % % % % % % % % % % % %         c_rgb = [1,0,0];rand([1,3]);
                % % % % % % % % % % % %         markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                % % % % % % % % % % % %             interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                % % % % % % % % % % % %             ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
                % % % % % % % % % % % %         scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
                % % % % % % % % % % % %         scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
                % % % % % % % % % % % %         title('Grouped genes');
                % % % % % % % % % % % %         axis tight
                % % % % % % % % % % % %         axis equal
                % % % % % % % % % % % %         axis off
                
                table1 = [cellid_sorted,m2c(double(inpos'))];
                saveCellFile(table1,[num2str(grp),'_exclust.txt']);
            end % move over groups
            % % %
            
        end % if groupy==1
        
    else  %vclust
        disp('reordering spatial barcodes with XY')
        xst=cell2mat(bar_ar(:,2));    yst=cell2mat(bar_ar(:,3));
        bar_ar_sorted=cell(length(cellid_sorted),3);
        for ib=1:length(cellid_sorted)
            curr_barcode=string(cellid_sorted(ib));
            bar_ind=find(curr_barcode==bar_ar(:,1));
            bar_ar_sorted(ib,:)={cellid_sorted(ib),xst(bar_ind),yst(bar_ind)};
        end
        gabaglut_sc=[];
    end
    %%
    
    if flags==1000
        %% Samples pie
        disp("Sample pie")
        figure('Name','P-Samples','NumberTitle','off');
        [p2,~]=numSubplots(size(gads,1));
        
        set(gcf,'color','w')
        
        
        [Gx,Gz] = groupcounts(sample_sorted(:,1));
        %         Gx=(Gx,1);
        %         subplot(p2(1),p2(2),ff)
        pie(Gx);
        title('All samples')
        axis off; axis tight; axis equal;
        
        legend(Gz)
        if savefig_flag==1
            savefig(gcf,'P-sample.fig')
        end
        %% samples bar
        disp("Bar Samples")
        num_clust=length(unique(T_cells_tmp));
        all_flags_sorted= [all_flags.fc_time_sorted];
        
        % samples
        % CLUTERS
        figure('Name','B-sample','NumberTitle','off');
        set(gcf,'color','w')
        sampleg=unique(sample);
        Ga=zeros(num_clust,length(sampleg));
        
        for acl=1:num_clust
            %     subplot(p2(1),p2(2),acl)
            [Gx,Gy] = groupcounts(sample_sorted(T_cells_tmp==acl));
            %     pie(Gx');
            %     title(num2str(acl))
            Gnorm=zeros(1,length(unique(sample)));
            Gnorm(ismember(sampleg,Gy))=100*Gx'./sum(Gx);
            %          if Gnorm==100
            %              Gnorm=[100,0];
            %
            %          end
            Ga(acl,:)=Gnorm;
        end
        cmapx=distinguishable_colors(length(sampleg));
        hc=bar(Ga,'stacked');
        for cxc=1:length(sampleg) % change bar colors distinguishly
            hc(cxc).FaceColor=cmapx(cxc,:);
        end
        title('Cluster by Sample')
        xt = get(gca, 'XTick');                                                         % Get Y-Tick Values
        xtn=1:num_clust;                                                       % New Y-Tick Labels
        set(gca, 'XTick',1:max(xtn), 'XTickLabel',xtn)
        ylim([0,100])
        xtickangle(45)
        legend(Gz)
        if savefig_flag==1
            savefig(gcf,'B_Sample.fig')
        end
        
        %% flags bar
        disp("Bar Flags")
        
        % samples
        % CLUTERS
        
        for ff=starter: size(all_flags,2)
            figure('Name',['B-Flags_',num2str(ff)],'NumberTitle','off');
            set(gcf,'color','w')
            legflg=string(unique(are_samples(samples_idx,findx(ff))));
            Ga=zeros(num_clust,2);
            for acl=1:num_clust
                %     subplot(p2(1),p2(2),acl)
                [Gx,Gy]= groupcounts(all_flags_sorted((T_cells_tmp==acl),ff));
                
                %     pie(Gx');
                %     title(num2str(acl))
                %          Gnorm=flip(Gnorm,2)
                Gnorm=zeros(1,2);
                Gnorm(ismember(1:length(legflg),Gy))=100*Gx'./sum(Gx);
                
                Ga(acl,:)=Gnorm;
            end
            hc=bar(Ga,'stacked');
            
            title(['Cluster by Flag: ',num2str(ff)])
            xt = get(gca, 'XTick');                                                         % Get Y-Tick Values
            xtn=1:num_clust;                                                       % New Y-Tick Labels
            set(gca, 'XTick',1:max(xtn), 'XTickLabel',xtn)
            ylim([0,100])
            xtickangle(45)
            legend(legflg)
            if savefig_flag==1
                savefig(gcf,['B-Flag_',num2str(ff),'.fig'])
            end
        end
        
        
        %% genes bar
        disp("Bar Genes")
        
        for ff=1:size(gads,1)
            figure('Name',['B-Genes ',num2str(ff)],'NumberTitle','off');
            set(gcf,'color','w')
            Ga=zeros(num_clust,2);
            sel_genes=gads'>0;
            for acl=1:num_clust
                [Gx,Gy] = groupcounts(sel_genes((T_cells_tmp==acl),ff));
                Gnorm=zeros(1,2);
                Gnorm(ismember([0 1],Gy))=100*Gx'./sum(Gx);
                Ga(acl,:)=flip(Gnorm,2);
            end
            hc=bar(Ga,'stacked');
            hc(1).FaceColor='b';
            hc(2).FaceColor='w';
            title('Cluster by Genes')
            xt = get(gca, 'XTick');                                                         % Get Y-Tick Values
            xtn=1:num_clust;                                                       % New Y-Tick Labels
            set(gca, 'XTick',1:max(xtn), 'XTickLabel',xtn)
            ylim([0,100])
            xtickangle(45)
            legend(ex_clust(ff))
            if savefig_flag==1
                savefig(gcf,['B-Gene_',num2str(ff),'.fig'])
            end
        end
        
        
        
        %% flags pie
        disp('Flags pie')
        all_legflag={};
        [p2,~]=numSubplots(size(all_flags,2));
        figure('Name','P-Flags','NumberTitle','off');
        set(gcf,'color','w')
        for ff=starter: size(all_flags,2)
            legflg=string(unique(are_samples(samples_idx,findx(ff))));
            subplot(p2(1),p2(2),ff)
            [Gx,Gy] = groupcounts(all_flags_sorted(:,ff));
            Gx=flip(Gx,1);
            pie(Gx');
            title(['Flag',num2str(ff)])
            axis off; axis tight; axis equal;
            
            legend(legflg)
            all_legflag=[all_legflag,{legflg}];
            
            
        end
        if savefig_flag==1
            savefig(gcf,'P-Flags.fig')
        end
        %% Genes pie
        disp("Genes pie")
        figure('Name','P-Genes','NumberTitle','off');
        [p2,~]=numSubplots(size(gads,1));
        
        set(gcf,'color','w')
        %      for ff=1:size(gads,1)
        
        for ff=1: size(gads,1)
            subplot(p2(1),p2(2),ff)
            [Gx,Gy] = groupcounts(sel_genes(:,ff));
            Gx=flip(Gx,1);
            pie(Gx);
            title(['Gene:',char(ex_clust(ff))])
            axis off; axis tight; axis equal;
            
            legend(ex_clust(ff))
            
            
        end
        if savefig_flag==1
            savefig(gcf,'P-Genes.fig')
        end
        %     end
    else
        all_flags_org=[];
        all_flags_sorted=[];
        num_clust=length(unique(T_cells_tmp));
        all_legflag=[];
    end

    %% save
    set(0,'DefaultFigureWindowStyle','normal')
    warning('on','all')
    %
    bar_ar_sorted=[];
    gabaglut_sc=[];
    s_xy_sorted=mapped_xy;
    %     all_flags_sorted= [all_flags.fc_time_sorted];
    datex=date;
    geneid_small=geneid;
    findx=findx-2;
    grp=0; % temporary
    cd(direct)
    disp('Saving mats... this may take a while... please wait...')
    %
    disp('Saving Sorted...')
    save('Sorted','data_orig_all_sorted','data_sorted_all','geneid_small','cellid_sorted','sample_sorted','bar_ar_sorted','s_xy_sorted','T_cells_tmp','all_flags_sorted','gabaglut_sc','-v7.3');
    %
    disp('Saving Cluster Parameters...')
    save('Clust_param','clmethod','num_clust','nrmz','MinPts','eps_prc','ex_genes','ex_clust','list2','maximal_total_molecule','initial_dims','perplexity','valm','valg','valmx','perct_90th_molecule','listy','datex','colors','fct','findx','sample_uni','all_legflag','filtery','fetsel','batcorr','-v7.3');
    %
    disp('Saving Orignal...')
    save('Orignal','data_orig_all','data_org','geneid_all','cellid_org','sample_org','bar_ar','s_xy_org','T_cells_all','all_flags_org','-v7.3');
    % sidx
    save('sidx.mat','sidx')

    
    tEnd = toc(tStart);
    fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    
    % end
    disp('Finished')
    
    if stops==0
        disp('Figures & variables will be cleared in 1 min - to cancel: ctrl+c')
        pause(60)% 1 min and delete
    else
        disp('Figures & variables will be cleared in 10 min - to cancel: ctrl+c')
        pause(60*10)% 10 min and delete
    end
    clear all
    close all
    clc
    disp('Finished')
end