%% SCRATCH: code for scRNAseq and visseq
%% construction
% fprintf(2,'Code under-construction, please try again later \nMuhammad \n')
% return;
%% bgn
clear all
close all
warning('off','all')
uuu=cd;
% cd('/data/runs/samples/')
set(0,'DefaultFigureWindowStyle','normal')

%% get parameters

try
    transcendentalDoc = '';
    clc
    S = {...
        struct('name','Cluster Type','type','enum','values',{{'cClust','vClust','CatClust'}},'doc',transcendentalDoc);...
        struct('name','Show(first=counter)','type','str','default','Snap25,Stmn2,Gad2,Slc32a1,Slc17a7,Slc17a6,Sst,Sim1,Foxj1,Pdgfra,Mog,C1qc,Flt1,Cldn5,Aqp4,Plp1,');...
        struct('name','Group Genes','type','str','default','C1qc,C1qa,C1qb,Gja1,Cx3cr1,Acta2,Ly6c1,Mfge8,Plp1,Aqp4,Vtn,Cldn5,Pdgfrb,Flt1,Slc1a3,Pdgfra,Foxj1,Olig1,Olig2,Sox10,Hbb-bs,Hbb-bt,Hba-a2,');...
        struct('name','exclude genes','type','str','default','Xist,Tsix,Eif2s3y,Ddx3y,Uty,Kdm5d,Btg2,Jun,Egr4,Fosb,Junb,Gadd45g,Fos,Arc,Nr4a1,Npas4,Coq10b,Tns1,Per2,Ptgs2,Rnd3,Tnfaip6,Srxn1,Tiparp,Ccnl1,Mcl1,Dnajb5,Nr4a3,Fosl2,Nptx2,Rasl11a,Mest,Sertad1,Egr2,Midn,Gadd45b,Dusp6,Irs2,Plat,Ier2,Rrad,Tpbg,Csrnp1,Peli1,Per1,Kdm6b,Inhba,Plk2,Ifrd1,Baz1a,Trib1,Pim3,Lrrk2,Dusp1,Cdkn1a,Pim1,Sik1,Frat2,Dusp5,VGF');...
        struct('name','Category Genes(1|2+3)','type','str','default','Gad2,Slc17a7,Slc17a6,');...
        struct('name','Normalization','type','int','default',5000);...
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
        struct('name','NAN|dbl|zero out','type','checkbox','default',0);...
        struct('name','Class non-neuro out','type','checkbox');...
        struct('name','Batch Correction','type','int','default',1);...
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

close all
clc

%% loadings and swithes
if 1 % to run all sections
    %% CHANGE HERE
    % exttra for bat and mouse
    disp('extra m&b')
    direct='/data/Technion_analysis/zebrafish/sc_100410/comparative/COMBINED';% get save folder
    cd(direct)
    %     load('/data/Technion_analysis/goldfish/scRNAseq_gf/n_mg_10x.mat')
    load([direct,'/n2_mzf_10x.mat'])

    %         load('/data/Technion_analysis/bat/NN_mouse_bat_harmony/nn_mb_10x.mat')

    sample=sampleid;
    data=all_data;
    listy={'1'}
    %
    geneid_org=geneid;
    data_org=data;
    cellid_org=cellid;
    sample_org=sample;
    cd(direct)
    %         disp('saving loads mat..')
    % % %         save('Loads','data_org','geneid','cellid_org','bar_ar','sample_org','geneid_org','cellid_org','-v7.3');


    listy=unique(sample);
    sample_uni=listy;



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
        all_flags(1).fc_time=flag_mgf;
        all_flags(2).fc_time=all_name;
        all_flags(3).fc_time=flag_g76;
        all_flags(4).fc_time=flag_mgfg76;
        all_flags(5).fc_time=flag_rgn;
        all_flags(6).fc_time=flag_loc;
        all_flags(7).fc_time=flag_gen;
    end

    %% distribution
%     set(0,'DefaultFigureWindowStyle','docked');
    figure('Name','distribution','NumberTitle','off');
    while 1
        tStart = tic;
        disp("Normalizing...")
        tot_mol = sum(data);
        maximal_total_molecule=max(tot_mol)
        perct_90th_molecule=round(prctile(tot_mol,90));
        nrmz=perct_90th_molecule % dynamic normalize
        tot_mol(tot_mol>nrmz) = nrmz;%%what value?
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

    if clusty ~="CatClust"
        %     data = normalize(data,'norm');
        %     data = round(data./repmat(sum(data),length(data(:,1)),1)*nrmz);
    end

    % clear for memory
    % save('Orignal','data_org','cellid_org','sample_org','bar_ar','-v7.3');
    clear  barcodesx datax geneidx
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% VALIDCELLS
    % stmn2 = data(strcmpi(geneid,'Stmn2'),:);
    % snap25 = data(strcmpi(geneid,'Snap25'),:);
    % exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
    %     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
    % % % exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plxnb3','Cldn11'....
    % % %     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc25a18','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2','Ttr'};
    % % % [~,loc] = ismember(exclude_markers,geneid);
    % validcells = (tot_mol>4000 & tot_mol<5e4 & tot_genes>2000 & sum(data(loc,:)>0)<=2 & (stmn2>0 | snap25>0) & amy_flag');
    disp('filtring non-valid cells...')
    %     validcells = (tot_mol>valm & tot_mol<valmx & tot_genes>valg); % all cells%5 value??
    validcells=1:size(data,2); % !!!!!!!!!!!!!!!!!!!!!!!! NO VALID CELLS
    for i=1:length(sample_uni)
        fprintf(['valid cells in ',sample_uni{i},' = ', num2str(sum(validcells(strcmpi(sample,sample_uni{i})))),'\n']);
    end
    sum(validcells>0)
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
    tot_mol(tot_mol>nrmz) = nrmz; % ??why twice
    tot_genes = sum(data>0);
    %

    data_orig_all = data;
    geneid_all = geneid;

    %% feature selection
    % % % %     disp("Excluding genes...")
    % % % %     ex_genes
    cd(direct)
    geneid=geneid_all;
    data=data_orig_all;
    if classy==0

        disp("Feature selection...")
        in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),ex_genes));
        fetsel="InterSpecies";
        if fetsel=="Default"
            corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);  % org
            data = data(in(corr_filt),:);
            geneid = geneid(in(corr_filt));
        elseif fetsel=="InterSpecies"
            batchid = all_flags(1).fc_time;%  1 flag for the species !!
            in_m = find(sum(data(:,batchid==1)>0,2)>5 & sum(data(:,batchid==1)>0,2)<length(data(1,batchid==1))*0.5 & ~ismember(geneid(:,1),ex_genes));
            in_b = find(sum(data(:,batchid==2)>0,2)>5 & sum(data(:,batchid==2)>0,2)<length(data(1,batchid==2))*0.5 & ~ismember(geneid(:,1),ex_genes));
            corr_filt1 = cv_vs_m_selection(data(in_m,batchid==1),geneid(in_m),[2000],1,0);
            corr_filt2 = cv_vs_m_selection(data(in_b,batchid==2),geneid(in_b),[2000],1,0);
            %             corr_filt=[in_m(corr_filt1);in_b(corr_filt2)];
            corr_filt=intersect(in_m(corr_filt1),in_b(corr_filt2));
            %             corr_filt=unique(corr_filt);
            data = data(corr_filt,:);
            geneid = geneid(corr_filt);
        elseif fetsel=="mt_tmp"
            batchid = all_flags(1).fc_time;%  1 flag for the species !!
            in_m = find(sum(data(:,batchid==1)>0,2)>5 & sum(data(:,batchid==1)>0,2)<length(data(1,batchid==1))*0.5 & ~ismember(geneid(:,1),ex_genes));
            in_b = find(sum(data(:,batchid==2)>0,2)>5 & sum(data(:,batchid==2)>0,2)<length(data(1,batchid==2))*0.5 & ~ismember(geneid(:,1),ex_genes));
            [TC,TI,TX]=unique(string(all_flags(2).fc_time),'stable');
            in_mg=unique([in_m;in_b]);
            [ind_gr_tmp_mark_fs,vi0] = markertablefeatures_tmp(TX,data(in_mg,:),10,1);
            corr_filt=ind_gr_tmp_mark_fs(:);
            corr_filt=unique(in_mg(corr_filt));
            data = data(corr_filt,:);
            geneid = geneid(corr_filt);
            geneidfs=geneid;
        elseif fetsel=="Alt" % new filter
            corr_filt = m_vs_f_selection(data(in,:),geneid(in),[],1,0);
        else % "Spatial"
            pdist_vis
        end


        fs_gen=geneid;

        %% gene autocorelation

        % batch_genes = {'Fau','Rps27','Rps12','Rpl10','Rpl27','Rpl23a','Pam16','Amd1','Uba52','Rpl35','Rpl7a','Gm10076','Lars2',....
        %     'Cxx1b','Erdr1','Entpd6','AY036118','Eno1','Rplp0','Pigt'};
        % sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};
        disp('auto-genecorlation...')
        maxcl=50;
        geneid_org=geneid;
        data_org=data;
        r=1;
        while r==0
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

            %     strx=lower(geneid(leaforder));
            strx=geneid(leaforder);
            % for making the
            %     sidx=regexp([' ' strx],'(?<=\s+)\S','start')-1
            %     strx(sidx)=upper(strx(sidx))
            %     propermsg = strcat(upper(strx(1)),lower(impropermsg(2:end)));
            table1 = [strx,m2c(idx(leaforder))];
            saveCellFile(table1,[direct,'/genes_correlation.csv'])
        end
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
    %     initial_dims = 100;
    prj = prj(:,1:initial_dims);
    prjx=prj;
    prj_all{1}=prjx;

    %% tsne 0 (optinal)
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
   mapped_xy_start=mapped_xy;

    % this is just the initial tsne, can be commented later
    clear Dsort
    %% HARMONY python addon
    if batcorr>0
        num_har=5;
        terminate(pyenv)
        batchid = all_flags(1).fc_time;%  1 flag for the species !!
        sd(1)=1000000000000000;
        prj=prjx;
        for ph=1:num_har
            ph
            pyenv("ExecutionMode","OutOfProcess");
            disp('Batch correction with Harmony py...')
            hprj = prj;
            usepylib = 1;
            [sout]=harmonypy(hprj,batchid,usepylib);
            prj = sout;
            terminate(pyenv)
            %             if ph>1
            %                 if  sd(ph)<=sd(ph-1)
            %
            prj_all{ph+1}=prj;
            %                     px1= pdist(prj_all{ph}); % local delta for x
            %                     px2= pdist(prj_all{ph+1}); % local delta for x+1
            %                     sd(ph+1)=sqrt(sum((px1 - px2) .^ 2));% global delta
            %                 else
            %                     disp('!break!')
            %                     ph
            %                     break;
            %
            %                 end
            %             else
            %                 sd(ph+1)=100000000000;
            %                 prj_all{ph+1}=prj;
            %             end
        end
        num_har=ph;
    
        %% check & plot delta between the prjs
        cd(direct)

        for ph=1:(num_har-1)
            ph
            px1= pdist(prj_all{ph}); % local delta for x
            px2= pdist(prj_all{ph+1}); % local delta for x+1
            sd(ph)=sqrt(sum((px1 - px2) .^ 2));% global delta

            figure('color','w');
            plot(1:length(sd),sd,'k');
            xlabel('H-H iteration'); ylabel('delta');

            if savefig_flag==1
                savefig(gcf,[direct,'/hrmny_itr.fig']) % take property name
            end
        end
    end  % batcorr
            %% choose projection and tsne

            prj=prj_all{end};% {itr+1} % CHANGE HERE!! % uncomment when not in a
            %     loop
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
%             mapped_xy_all{ph+1}=mapped_xy;

            % this is just the initial tsne, can be commented later
            clear Dsort
        
%     end % batcorr
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
    colors = distinguishable_colors(length(sample_uni)+1);

    %% plot  per flags
    % plot by CONDITION
    disp('plot flags')
    if flags==1
        for ii=[1:7]%:size(all_flags,2) % 2 is all cell types not needed
            % % %             legflg=string(unique(are_samples(samples_idx,findx(ii))));
            fc_time_uni=unique(all_flags(ii).fc_time,'stable');
            colors =  distinguishable_colors(length(fc_time_uni));
            colors(end,:)=[0.5 0.5 0.5];
            [px,~]=numSubplots(length(fc_time_uni));
            %             figure('Name',['per flag',num2str(ii)],'NumberTitle','off');
            %             set(gcf,'color','w')
            %
            %             [ha, pos] = tight_subplot(px(1),px(2), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
            %             for i=1:length(fc_time_uni)
            %                 axes(ha(i))
            %                 s = dscatter(mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),1),mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),2));%,10,colors(i,:),'filled');
            %                 axis tight
            %                 axis equal
            %                 axis off
            % %                 legend('m');
            %                 title(num2str(length(mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),1))))
            %             end
            %             if savefig_flag==1
            % % %                 savefig(gcf,['tsne_perflag_',char(flist(findx(ii)-2)),'.fig']) % take property name
            %                 savefig(gcf,['tsne_per_flag_',num2str(ii),'.fig'])
            %
            %             end
            %  % plot  all flags
            figure('Name',['flag ts',num2str(ii)],'NumberTitle','off');
            set(gcf,'color','w')

            for i=1:length(fc_time_uni)
                hold on
                s = scatter(mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),1),mapped_xy(all_flags(ii).fc_time==fc_time_uni(i),2),10,colors(i,:),'filled'); hold on;
            end
            legend
            axis tight
            axis equal
            axis off
            % % % %             title(flist(findx(ii)-2))
            % legend('b','m') % take legend name
            if ii==1
                legend({'Mouse1','Mouse2','Zebrafish'},'Location','northeastoutside')
            elseif ii==3
                legend({'GABA','GLUT'},'Location','northeastoutside')
            elseif ii==4
                legend({'m-GABA','m-GLUT','g-GABA','g-GLUT'},'Location','northeastoutside')
            else
                legend(string(fc_time_uni),'Interpreter', 'none','Location','northeastoutside') % take legend name
            end
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


    %% save
    disp('save')
    save([direct,'/m_vs_gf_h.mat'],"geneid","prj_all","perplexity",'mapped_xy','mapped_xy_start','all_flags','num_har','geneid_all','data_orig_all','data','-v7.3')
    set(0,'DefaultFigureWindowStyle','normal')
end