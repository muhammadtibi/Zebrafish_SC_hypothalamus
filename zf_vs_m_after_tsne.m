%% %% %% %% after scratch (before dbscan (no clustring))
% cd('/data/Technion_analysis/goldfish')
% save("/data/Technion_analysis/goldfish/m_vs_gf_h.mat","prj",'mapped_xy','all_flags','data_org','geneid_all','data','-v7.3'); % prj&xy=> H0

%% load
% cd('/data/Technion_analysis/goldfish')
clear all
close all
clc
%%
direct='/data/Technion_analysis/zebrafish/sc_100410/comparative/COMBINED';
cd(direct);
load('/bigdata/all_mouse_dataset/mn_markeragg.mat','mn_markeragg','mn_clagg')
load([direct,'/m_vs_gf_h.mat'],"prj_all","perplexity",'geneid',"num_har","all_flags","data_orig_all","geneid_all","data",'mapped_xy','mapped_xy_start');% prj is harmony 0
% load([direct,'/lvidx.mat'],'idx')
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureVisible','on');% off / on
% load('clr.mat','clr','usort')
% choosen_run=5;
% mapped_xy=mapped_xy_start; %  uncomment start = no harmony !!!!!!!!!!!!!

prj=prj_all{end};
% perplexity=47;% !!
% lv2ct=table2array(readtable('Tlv2ct.csv'));
%%

% num_har=5;
nms=all_flags(2).fc_time;
mgf=all_flags(1).fc_time;
loc=all_flags(5).fc_time;
rgn=all_flags(6).fc_time;
agcls=mn_clagg;
agnms=mn_markeragg;
ugf=unique(nms,'stable');
n_gi=contains(ugf,'zf-');
n_mx=ugf(~n_gi);
n_m= regexprep(n_mx,'_','-');
n_gx=ugf(n_gi);
n_g= regexprep(n_gx,'_','-');
n_g= regexprep(n_g,'zf-','');
loci=strings(length(n_m),1);
rgni=strings(length(n_m),1);
geni=strings(length(n_m),1);

for i=1:length(n_m)
    ni= find(nms==n_m(i));
    loci(i)=unique(loc(ni));
    rgni(i)=unique(rgn(ni));
    agi= find(agcls==n_m(i));
%     geni(i)=agnms(agi);
end

% create color vector for gf cell types


r1=n_g(contains(n_g,'GABA'));
% make colormap
gabacol=winter(length(r1));


r2=n_g(contains(n_g,'Glut-'));

if isempty(r2)
    prefix='GABA-';
else
    prefix='Glut-';
end
glutcol=spring(length(r2));

gabaglut=[glutcol;gabacol];
%% optional : graph-dendro method
grphix=1;
if grphix==1

    %% clustering (lovian)- old 
    disp('Clustering...')
    markfig=figure('Name','preview','NumberTitle','off');
    set(gcf,'color','w')

    eps_prc=1;
    adj = knn_mat(prj,perplexity,2);
    addpath('/data/matlab_functions/louvain_scgeatoolbox/')
    [com, Q, comty] = louvain(adj);

    partition=length(comty.COM)
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
        save('lv_cluter_param',"idx","MinPts","eps_prc")

    %% dbscan new 
    eps_prc=45;
    MinPts=60;
    dbscanit
    save('db_cluter_param',"idx","MinPts","eps_prc")
     %% or KNNmat instead (tranfer label)
        load('CustomColormap.mat')
%         mapped_xy=prj; % c/uc to tsne/pca calc. 
        drang=50;
        mm=sum(mgf==1);
        m_xy=prj(1:mm,:);
        g_xy=prj((mm+1):end,:);
        disty=pdist(prj);
        dd=squareform(disty);
        subdd=dd(1:mm,(mm+1):end);
        [sdd,idd]=sort(subdd,1,'ascend');
        iddtop=idd(1:drang,:);
        sddtop=sdd(1:drang,:);
        [ga,gb,gc]=unique(nms,'stable');
        inum=zeros(size(subdd,2),1);
        istr=strings(size(subdd,2),1);
        for ii=1:length(g_xy)
         subi=iddtop(:,ii);
         maj=mode(gc(subi));
         inum(ii)=maj;
         istr(ii)=ga(maj);
        end
    %  
    % create a ratio matrix (for KNN) MM VS GF (NO AUTO-compare)

    lv2ct=zeros(length(n_m),length(n_g));
    gfnms=nms((mm+1):end);
    for j=1:length(n_g)
        jj=find(gfnms==n_gx(j));
        inn=inum(jj); % take only cyurrent gf ct maj mm ct -idx
        iss=istr(jj);  % take only cyurrent gf ct maj mm ct -num
        [iy,~,ix] = unique(inn,'stable');
        C = accumarray(ix,1);
        lv2ct(iy,j)=C/sum(C);
    end
     all(lv2ct,2) 

      
    Tlv2ct=table([["cell type",n_g'];[n_m,lv2ct]]);
    writetable(Tlv2ct,[direct,'/',prefix,'Tknn2ct.csv']);
    %
    % create hm for knnmat
    
    figure('color','w');
    if drang==50
    [~,mxi]=max(lv2ct');
    [~,is]=sort(mxi);
    end
    imagesc(lv2ct(is,:),[0 0.5]);
    colormap(CustomColormap)
    set(gca,'ytick',1:length(n_m),'yticklabel',strcat(rgni(is),'-',n_m(is)), 'fontsize', 10)
    set(gca,'xtick',1:length(n_g),'xticklabel',n_g, 'fontsize', 10)
    colorbar
    title(drang)
  %% check cluster ratios! (clusters:dbscan | louvian)

    for i=1:length(n_m)
        ni= find(nms==n_m(i));
        loci(i)=unique(loc(ni));
        rgni(i)=unique(rgn(ni));
        agi= find(agcls==n_m(i));
%         geni(i)=agnms(agi);
    end
    locia=[loci;strings(length(n_g),1)]';
    lv2ct=zeros(max(idx),length(ugf));
    for j=1:length(ugf)
        jj=find(nms==ugf(j));
        [iv,~,ix] = unique(idx(jj));

        C = accumarray(ix,1);
        lv2ct(iv,j)=C/sum(C);
    end

    
    %% create dendrogram
    n_g= regexprep(n_g,'Glut-','');
    n_g= regexprep(n_g,'GABA-','');
    appended=append(prefix,string(1:length(n_g)));
    nugf=[n_m;appended'];
    % y= clusters
    Zpcax = linkage(lv2ct,'ward','correlation');
    Dx = pdist(lv2ct,'correlation');
    leafordery = optimalleaforder(Zpcax,Dx);
    cmap=turbo(length(leafordery));
    % x= cell types
figure('color','w')

%     axes('position',[0.03,0.1,0.3,0.8])
        subplot(1,2,1)
    hden = dendrogram(Zpcax,length(leafordery),'Reorder',leafordery,'Orientation','left');
     axis off; 
    set(gca,'ylim',[0.5,length(leafordery)+0.5])
    Zpca = linkage(lv2ct','ward','correlation');
    Dx = pdist(lv2ct','correlation');
    leaforderx = optimalleaforder(Zpca,Dx);
    for ii=1:length(leafordery)

    subplot(1,2,2)
     axis off; 
     hold on
    scatter(1,ii-0.5,100,cmap(ii,:),'filled')
    text(1.1,ii-0.5,num2str((length(leafordery)+1)-ii));
    end 
    Tlv2ct=table([["cell type","location",string(1:max(idx))];[ugf,locia',lv2ct(leafordery,:)']]);
    writetable(Tlv2ct,[direct,'/',prefix,'Tlv2ct.csv']);

    %% hm 
figure('color','w')
subplot(1,2,1)

    hden = dendrogram(Zpcax,length(leafordery),'Reorder',leafordery,'Orientation','left');
        set(gca,'ylim',[0.5,length(leafordery)+0.5])

axis off
subplot(1,2,2)

    imagesc(lv2ct(leafordery,leaforderx));
    colormap('summer')
    set(gca,'xtick',[1:length(nugf)],'xticklabel',nugf(leaforderx), 'fontsize', 7);
    set(gca,'ytick',[1:length(leafordery)],'yticklabel',[1:length(leafordery)], 'fontsize', 7);%leafordery
%         set(gca,'ytick',[1:length(leafordery)],'yticklabel',leafordery, 'fontsize', 7);%org

    
    % figure loc legend

%%  clusterts sorted figure
figure('color','w');
axis tight;
axis equal;
axis off;
hold on
t_cells_tmp=zeros(size(idx));
  for i=1:max(idx)

        ii=find(idx==leafordery(i));
        t_cells_tmp(ii)=i; 
        xc=mapped_xy(ii,1);
        yc=mapped_xy(ii,2);
        scatter(xc,yc,30,cmap(length(leafordery)+1-i,:),'.');
        scbound = boundary(xc,yc,0.1);
%         plot(xc(scbound),yc(scbound),'--','LineWidth',2,'Color','k');% show border line
%         inx = idx==i;% show textbox
        ht = text(median(xc),median(yc),num2str(i));
        set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',8);
        set(ht,'Color','k')
  end

 %% heatmap with top genes 
    data_sorted_all=data_orig_all;
    [sidx,ssidx]=sort(t_cells_tmp,'ascend');
    t_cells_tmpx=t_cells_tmp(ssidx);
    [ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(t_cells_tmpx,data_sorted_all(:,ssidx),3);% top_genes by me
    % % % % % % % % %
    datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
      Zpcay = linkage(datamarkers,'ward','correlation');
    Dx = pdist(datamarkers,'correlation');
    leafordery = optimalleaforder(Zpcay,Dx);
    % gr_tmp_mark = gr_tmp_mark(xi);
    gr_tmp_mark = geneid_all(ind_gr_tmp_mark);
     datamarkers= datamarkers(leafordery,:);
    datamarkers_cn = cent_norm(log2(datamarkers+1));

    figure('Name','HeatMap','NumberTitle','off');
    set(gcf,'color','w')
%     ax1 = axes('position',[0.1,0.02,0.88,0.73]);
    imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
    hold on;

hold on
alltmp=1:size(data_sorted_all,2);
idtmp=alltmp(logical(diff(t_cells_tmpx)));
idtmp=[1,idtmp];
% idtmp(end)=[];
xline(idtmp)

set(gca,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark(leafordery), 'fontsize', 10)
set(gca,'xtick',gr_center,'xticklabel',1:max(t_cells_tmpx), 'fontsize', 10)
colormap('summer');

 

    %%  shatttered hm
    figure('Color','w')
    colormap(summer)
    maxima=max(max(lv2ct));
    ha = @(m,n,p) subtightplot (m, n, p,[.001 .02],[.001 .001],[.02 .01]);
    for hi=1:length(leafordery)
        ax(hi)=ha(1,length(leafordery),hi);
        [vstrip,istrip]=sort(lv2ct(leafordery(hi),:),'descend');
        imagesc(vstrip(1:10)');
        [R, C] = ndgrid(1:10, 1);
        R = R(:); C = C(:) - 1/4;
        text(C, R, string(round(vstrip(1:10),2)), 'color', 'k','FontSize',7)
        axis tight;axis equal;
        % title(leafordery(hi))
%       set(gca,'xtick',[1],'xticklabel',leafordery(hi), 'fontsize', 12);
        set(gca,'xtick',[1],'xticklabel',hi, 'fontsize', 12);
        set(gca,'ytick',[1:10],'yticklabel',nugf(istrip(1:10)), 'fontsize', 7);
        caxis manual
        caxis([0 1]);
    end
    % This sets the limits of the colorbar to manual for the second plot

    %% heatmap- subset !
    d_mat=pdist(lv2ct');
    d_mat=squareform(d_mat);
    hm=d_mat(1:length(n_m),length(n_m)+1:end);
    % y= clusters
    Zpcay = linkage(hm,'ward','correlation');
    Dx = pdist(hm,'correlation');
%     leafordery = optimalleaforder(Zpcay,Dx);
    % x= clusters
    Zpcax = linkage(hm','ward','correlation');
    Dx = pdist(hm','correlation');
    %     leaforderx = optimalleaforder(Zpcax,Dx);
    leaforderx=1:length(Zpcax);
    hm=hm(leafordery,leaforderx);
    figure('color','w');
    ax1=subplot(3, 4, [2:4,6:8,10:12]); % top row of the 3x3 grid
    imagesc(hm);

    L=10; %  number of data points
    indexValue = 0;     % value for which to set a particular color
    topColor = [1 1 1];         % color for maximum data value (red = [1 0 0])
    indexColor = [0.5 0 0];       % color for indexed data value (white = [1 1 1])
    bottomcolor = [1 0 0];      % color for minimum data value (blue = [0 0 1])
    % Calculate where proportionally indexValue lies between minimum and
    % maximum values
    largest = max(max(hm));
    smallest = min(min(hm));
    index = L*abs(indexValue-smallest)/(largest-smallest);
    % Create color map ranging from bottom color to index color
    % Multipling number of points by 100 adds more resolution
    customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
        linspace(bottomcolor(2),indexColor(2),100*index)',...
        linspace(bottomcolor(3),indexColor(3),100*index)'];
    % Create color map ranging from index color to top color
    % Multipling number of points by 100 adds more resolution
    customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
        linspace(indexColor(2),topColor(2),100*(L-index))',...
        linspace(indexColor(3),topColor(3),100*(L-index))'];
    customCMap = [customCMap1;customCMap2];  % Combine colormaps
    % colormap(customCMap)
    % psudo = pcolor(hm);
    colormap(ax1,customCMap)
    colorbar('location','eastoutside')
    set(gca, 'XTick',1:length(n_g), 'XTickLabel',extractAfter(n_g(leaforderx),'-'))
    set(gca ,'YTick',1:length(n_m), 'YTickLabel',n_m(leafordery))
    % set(gca,'xtick',[])
    xtickangle(45)
    grid on
    yx=1:length(n_m);
    hold on
    locia(locia=="")=[];
    [u,uu,uuu]=unique(locia(leafordery),'stable');

    ux=find(diff(uuu)~=0);

    for i=1:length(ux)

        % yline(yx(ux(i))+0.5,'k',u(uuu(ux(i))),'LineWidth',2,'LabelHorizontalAlignment','right')
        yline(yx(ux(i))+0.5,'k','LineWidth',2,'LabelHorizontalAlignment','right')

    end
    % yline(length(n_m)+0.5,'k',u(uuu(end)),'LineWidth',2,'LabelHorizontalAlignment','right')% for last one
    yline(length(n_m)+0.5,'k','LineWidth',2,'LabelHorizontalAlignment','right')% for last one

    ax2=subplot(3,4,[1 5 9]);
    imagesc(uuu)
    [ia,ib,ic]=intersect(u,usort);
    clr=(double(usort(2:4,:))/255)';
    colormap(ax2,clr(ic,:))
    for i=1:length(ux)

        yline(yx(ux(i))+0.5,'k',u(uuu(ux(i))),'LineWidth',2,'LabelHorizontalAlignment','right')
        % yline(yx(ux(i))+0.5,'k','LineWidth',2,'LabelHorizontalAlignment','right')

    end
    yline(length(n_m)+0.5,'k',u(uuu(end)),'LineWidth',2,'LabelHorizontalAlignment','right')% for last one
    % yline(length(n_m)+0.5,'k','LineWidth',2,'LabelHorizontalAlignment','right')% for last one

    axis off




end
%% scatter both species in parrallel (goldfish)
figure;
set(gcf,'Color','w')
% % % goldfish

% subplot(1,2,1)
scatter(mapped_xy(:,1),mapped_xy(:,2),30,[0.5, 0.5, 0.5],'.'); % mouse
hold on
% scatter(mapped_xy(mgf==2,1),mapped_xy(mgf==2,2),20,'b','.'); % mouse
all_medianxy=zeros(length(n_gx),2);

for i=1:length(n_gx)
    xgf=contains(nms,n_gx(i));
    scatter(mapped_xy(xgf,1),mapped_xy( xgf,2),30,gabaglut(i,:),'.'); % mouse
    all_medianxy(i,:)=median([mapped_xy(xgf,:)]);


end





textPosition=calcTextPosInScatterPlot(all_medianxy);
for x=1:length(n_gx)
    plot([textPosition(x,1) all_medianxy(x,1)], ...
        [textPosition(x,2) all_medianxy(x,2)],'k-','LineWidth',2);
    text(textPosition(x,1), textPosition(x,2), ...
        n_gx(x), 'HorizontalAlignment','center','FontSize',5,'FontWeight','bold');%, ...        % 'BackgroundColor',[0.8 0.8 0.8]

end
axis off;  axis equal;

% title('Goldfish')



%% % % mouse (supra-regions)
% load('legcol.mat', 'usort')
% usort(2,2) = "1";
% usort(3,2) = "159";
% usort(4,2) = "171";
clr=jet();%(double(usort(2:4,:))/255)';
% t=readtable([direct,'/',num2str(num_har),'_Tm2gf_rgn.csv']);
% mt=table2array(t(:,1));
% gt=table2array(t(:,5));
figure('color','w')
% subplot(1,2,2)
scatter(mapped_xy(:,1),mapped_xy(:,2),30,[0.5, 0.5, 0.5],'.'); % mouse
% scatter(mapped_xy(mgf==1,1),mapped_xy(mgf==1,2),20,'r','.'); % mouse


hold on
usort_all=[];
for i=1:length(usort)

    xm=strcmpi(loc,usort(1,i));
    % xt=strcmpi(mt,n_mx(i));
    % yt=contains(n_gx,regexprep(gt(xt),'-','_'));
    % scatter(mapped_xy(xm,1),mapped_xy(xm,2),30,gabaglut(yt,:),'.'); % mouse
    scatter(mapped_xy(xm,1),mapped_xy(xm,2),30,clr(i,:),'.'); % mouse
    % all_medianxy(i,:)=median([mapped_xy(xm,:)]);
    if sum(xm)>0
        usort_all=[usort_all;usort(1,i)];
    end

end

for ii=1:length(n_m)
    xm=strcmpi(nms,n_mx(ii));

    all_medianxy(ii,:)=median([mapped_xy(xm,:)]);

end

textPosition=calcTextPosInScatterPlot(all_medianxy);
for x=1:length(n_mx)
    plot([textPosition(x,1) all_medianxy(x,1)], ...
        [textPosition(x,2) all_medianxy(x,2)],'k-','LineWidth',2);
    text(textPosition(x,1), textPosition(x,2), ...
        n_m(x,:), 'HorizontalAlignment','center','FontSize',5,'FontWeight','bold');%, ...        % 'BackgroundColor',[0.8 0.8 0.8]
end

axis off;  axis equal;
% title('Mouse')
% legend(['';usort_all])

  for i=1:max(idx)

        ii=find(idx==leafordery(i));
        xc=mapped_xy(ii,1);
        yc=mapped_xy(ii,2);
%         scatter(xc,yc,30,cmap(length(leafordery)+1-i,:),'.');
        scbound = boundary(xc,yc,0.1);
        plot(xc(scbound),yc(scbound),'--','LineWidth',2,'Color','k');% show border line
%         inx = idx==i;% show textbox
%         ht = text(median(xc),median(yc),num2str(i));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',8);
%         set(ht,'Color','k')
  end
%%
% figure loc legend
figure('Color','w');
hold on
axis off;
for ii=1:length(usort)
    scatter(1,0+ii,100,clr(ii,:),'filled')
    text(1.1,0+ii,num2str(usort(1,ii)));
end
%% tsne visual (3)y
% plot(mapped_xy(:,1),mapped_xy(:,2),'.')
% bmcol=distinguishable_colors(2);
figure('Color','w');
scatter(mapped_xy(all_flags(1).fc_time==1,1),mapped_xy(all_flags(1).fc_time==1,2),30,[244,88,226]/255,'.')
hold on
scatter(mapped_xy(all_flags(1).fc_time==2,1),mapped_xy(all_flags(1).fc_time==2,2),30,[0.5,0.5,0.5]/255,'.')
% scatter(mapped_xy(all_flags(1).fc_time==3,1),mapped_xy(all_flags(1).fc_time==3,2),30,[104,188,134]/255,'.')
axis off; axis equal; axis tight;
legend('Mouse', 'Zebrafish')
%% tsne visual (2)
% plot(mapped_xy(:,1),mapped_xy(:,2),'.')
% bmcol=distinguishable_colors(2);
figure('Color','w');
scatter(mapped_xy(all_flags(1).fc_time==1,1),mapped_xy(all_flags(1).fc_time==1,2),30,[244,88,226]/255,'.')
hold on
scatter(mapped_xy(all_flags(1).fc_time==2,1),mapped_xy(all_flags(1).fc_time==2,2),30,[104,188,134]/255,'.')
axis off; axis equal; axis tight;
legend('Mouse', 'Zebrafish')
   %% create graph
    %     zlv2ct=lv2ct;
    %     zlv2ct(lv2ct<(1/max(idx)))=0;
    %     upmat=(lv2ct>(2/max(idx)));
    %     ids=1:length(ugf);
    %     adjmat=zeros(length(ugf),length(ugf));
    %     %     adjmat(upmat)=lv2ct(upmat);
    %     for i=1:length(ugf)
    %         i
    % %         [vpmat,ipmat]=max(lv2ct(:,i));
    %         adjmat(i,lv2ct(ipmat,:)>=0.5)=lv2ct(ipmat,lv2ct(ipmat,:)>=0.5);
    % %         adjmat(lv2ct(ipmat,:)>=0.5,i)=lv2ct(ipmat,lv2ct(ipmat,:)>=0.5)';
    %
    % %         [vss,iss]=sort(lv2ct(ipmat,:),'descend');
    %         %         [vss,iss]=max(lv2ct(ipmat,idss));
    % %         if iss(1)==i % if the row is also the max
    % %             ix=iss(2);
    % %             vx=vss(2);
    % %         else
    % %             ix=iss(1);
    % %             vx=vss(1);
    % %         end
    % %         adjmat(i,ix)=vx;
    % %         adjmat(ix,i)=vx;
    %
    %     end
    %     %     adjmat(logical(eye(length(adjmat))))=0;
    %     colors=zeros(length(ugf),3);
    %     colors([length(n_m)+1:end],:)=gabacol;
    %     G = graph(adjmat,'omitselfloops');
    %     figure('color','w');
    %     p=plot(G,'Layout','force','NodeColor',colors,'Nodelabel',[n_m;n_g]);
    %     axis off;
    %     % p.Marker = 's';
    %     highlight(p,1:length(ugf),'MarkerSize',15)
    %     p.NodeFontSize =10;

%% (sainity) plot gene of mouse and goldfish to view proximity
% NPT % c/uc
% genet=string(table2array(readtable('/data/Technion_analysis/zebrafish/sc_100410/np_tf_list.xlsx')));
% genex=upper([genet(:,1);genet(:,2);genet(:,3)]);
% or TF % c/uc
% genet=readtable('/data/Technion_analysis/zebrafish/sc_100410/TFs.xlsx');
% genex=upper(string(table2array(genet(:,1))));
% genex(cellfun('isempty',genex)) = []; % remove empty cells from mm
% genex=["GAD2","SLC32A1","GAD1","SLC17A7","SLC17A6"]; % or comment or umcoomment
genex=upper(["lhx6"
"crhbp"
"gal"
"galn"
"kiss1"
"sst"
"otpa"
"otpb"
"nts"
"trh"
"avp"
"nr5a2"
"dusp2"
"oxt"
"vip"
"lhx9"
"rfx4"
"hcrt"
"npvf"
"hmx3a"]);

for i=1:length(genex)
    i

% scatter(mapped_xy(:,1),(mapped_xy(:,2)),10,'.');
% ax1=subplot(1,2,2)

fg=genex(i);
try
gx=find(string(geneid_all)==fg);
markergene=data_orig_all(gx,:);
tmpthlow = 0;%prctile(markergene(markergene>0),30);
tmpthhigh = 1;%prctile(markergene(markergene>0),80);
markergene(markergene>tmpthhigh) = tmpthhigh;
markergene(markergene<tmpthlow) = tmpthlow;
c_rgb =[244,88,226]/255;% pink

markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
c_rgb = [104,188,134]/255;% green
markergene_color2 = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
% scatter(mapped_xy(mgf==1,1),mapped_xy(mgf==1,2),'.','r'); % mouse
% % hold on
% scatter(mapped_xy(mgf==2,1),mapped_xy(mgf==2,2),'.','b'); % goldfish
figure('color','w');
scatter(mapped_xy(mgf==1,1),mapped_xy(mgf==1,2),30,markergene_color((mgf==1),:),'.'); % mouse
axis off; axis tight; axis equal;
hold on
scatter(mapped_xy(mgf==2,1),mapped_xy(mgf==2,2),30,markergene_color2((mgf==2),:),'.'); % mouse
axis off; axis tight; axis equal;
title(fg)
namepdf=char(genex(i));
eval(['export_fig ',namepdf,'.pdf  -nocrop -r 300']);
end
close all
end % genes on tsne
%% XXXXX DO it with mxy : with prj( or t-sne) mean of mxy
% [fc_time_uni,~,fc_time_all]=unique(all_flags(2).fc_time,'stable');


loc(all_flags(5).fc_time=="")=nms(all_flags(5).fc_time=="");
fc_time_uni=unique(loc,'stable');
%  to go back put nms instead of loc
prjx=zeros(length(fc_time_uni),size(prj,2));
for xx=1:length(fc_time_uni)
    prjx(xx,:)=mean(prj((loc==fc_time_uni(xx)),:),1);
end

u_c = squareform(pdist(prjx,'euclidean'),'tomatrix');
% Dx=u_c(1:length(n_m),(length(n_m)+1):end); % uncomment to go back
Dx=u_c(1:9,10:end);


u_cx=Dx;
% u_cx=u_cx(1:length(n_mx),length(n_mx)+1:end);

%% (OPTINAL) RUN in ALL PRJ_ALL in  a LOOP
grp=strings(length(n_g),32);
grp(:,1)=n_g(:);
for p=1:length(prj_all)
    p
    prj=prj_all{p};
    fc_time_uni=unique(all_flags(2).fc_time,'stable');
    prjx=zeros(length(fc_time_uni),size(prj,2));
    for xx=1:length(fc_time_uni)
        prjx(xx,:)=mean(prj((all_flags(2).fc_time==fc_time_uni(xx)),:),1);
    end

    u_c = squareform(pdist(prjx,'euclidean'),'tomatrix');
    Dx=u_c(1:length(n_m),(length(n_m)+1):end);

    leafordery=1:length(n_m);

    Dy=Dx(leafordery,:);
    n_mx=n_m(leafordery);
    n_gx=n_g;
    u_cx=u_c;
    u_cx=u_cx(1:length(n_mx),length(n_mx)+1:end);
    [s_u,s_i]=sort(u_cx);
    for x=1:length(n_g)
        % grp(x,6)=s_u(1,x);
        grp(x,p+1)=n_m(s_i(1,x));
    end

end
% grp(:,2)=mode(categorical(grp(:,2:end)),2);
tgrp=table(grp);
writetable(tgrp,'all_Tgf2m_without_midbrain.csv')
%% heatmap the grp
grp_index=zeros(size(grp(:,3:end)));
for i=1:length(n_g)
    i
    [ia,ib,ic]=unique(grp(i,2:end));
    grp_index(i,:)=logical(diff(ic'));
end
figure;imagesc(grp_index)
colormap("summer")
xlabel('Harmony Itreration')
ylabel('Goldfish Cell Type')

numg=1:numel(n_g);
set(gca, 'YTick',numg, 'YTickLabel',n_g)
f_sum=find(sum(grp_index)==0);
num_har=f_sum(1)+1
%% plot matrix (full u_c , sub: u_cx)
% set(0,'DefaultFigureWindowStyle','normal')
figure('Color','w');
cmap=redblue(256);
colormap(flip(cmap,2))
u_cx=u_c;
% u_cx(55,:)=[]; % derbalek ana kemtehn I remove it because its problamatic
% u_cx(:,55)=[];
u_cx=u_cx(1:length(n_mx),length(n_mx)+1:end);
% Yup = prctile(u_cx(:),80);% thresholding imagesc
% Ydn = prctile(u_cx(:),50);
Z1 = linkage(u_cx,'ward','euclidean');
% leaforder1 = optimalleaforder(Z1,pdist(u_cx));
leaforder1 = 1:length(n_m);
Z2 = linkage(u_cx','ward','euclidean');
% leaforder2 = optimalleaforder(Z2,pdist(u_cx'));
leaforder2 = 1:length(n_g);
imagesc(u_cx(leaforder1,leaforder2),[0 20]);
xlabel('Goldfish')
ylabel('Mouse')
% ylabel('Bat')
% xlabel('Mouse')

numm=1:numel(n_m);
numb=1:numel(n_g);
set(gca, 'YTick',numm, 'YTickLabel',n_m(leaforder1))
set(gca ,'XTick',numb, 'XTickLabel',n_g(leaforder2))
xtickangle(45)
% xline(1,'-k',{'Mouse'},'LineWidth',2)
% xline(length(n_mx),'-k',{'Goldfish','GABA'},'LineWidth',2)
% xline(123,'-k',{'Goldfish','vGlut'},'LineWidth',2)
% yline(1,'-k',{'Mouse'},'LineWidth',2)
% yline(length(n_mx),'-k',{'Goldfish','GABA'},'LineWidth',2)
% yline(123,'-k',{'Goldfish','vGlut'},'LineWidth',2)
colormap;
%
% n_mby=n_mbx;
% n_mby(55)=[];
Tx=[n_g(leaforder2)';u_cx];
zz=([['Mouse\Goldfish';n_m(leaforder1)],['Marker';geni],['Location';loci],['Region';rgni]]);
Ty=[zz,Tx];
Tt=table(Ty);
writetable(Tt,'Tm2gf.csv')

%% for nms: group table m2gf
[s_u,s_i]=sort(u_cx,2);
grp=strings(length(n_m),6);
for x=1:length(n_m)
    grp(x,6)=s_u(x,1);
    grp(x,5)=n_g(s_i(x,1));
    grp(x,4)=geni(x);
    grp(x,3)=loci(x);
    grp(x,2)=rgni(x);
    grp(x,1)=n_m(x);

end
tgrp=table(grp);
writetable(tgrp,[num2str(num_har),'_Tm2gf_rgn.csv'])
%% for nms: group table gf2m
[s_u,s_i]=sort(u_cx);
grp=strings(length(n_g),16);
for x=1:length(n_g)
    grp(x,6)=s_u(1,x);
    grp(x,5)=n_m(s_i(1,x));
    grp(x,4)=geni(s_i(1,x));
    grp(x,3)=loci(s_i(1,x));
    grp(x,2)=rgni(s_i(1,x));
    grp(x,1)=n_g(x);
    grp(x,11)=s_u(2,x);
    % grp(x,10)=n_m(s_i(2,x));
    % grp(x,9)=geni(s_i(2,x));
    grp(x,8)=loci(s_i(2,x));
    grp(x,7)=rgni(s_i(2,x));
    grp(x,16)=s_u(3,x);
    grp(x,15)=n_m(s_i(3,x));
    grp(x,14)=geni(s_i(3,x));
    grp(x,13)=loci(s_i(3,x));
    grp(x,12)=rgni(s_i(3,x));
end
tgrp=table(grp);
writetable(tgrp,[num2str(num_har),'_loc_Tgf2m_rgn.csv'])

%% for locs: group table gf2m
[s_u,s_i]=sort(u_cx);
grp=strings(length(n_g),5);
for  x=1:length(n_g)
    grp(x,1)=n_g(x);
    grp(x,2)=fc_time_uni((s_i(1,x)));
    grp(x,3)=s_u(1,x);
    grp(x,4)=fc_time_uni(s_i(2,x));
    grp(x,5)=s_u(2,x);
end
tgrp=table(grp);
writetable(tgrp,[num2str(num_har),'_locs_Tgf2m_rgn.csv'])


%% categorize data
catl=["Amygdala","Cortex","Cortex","Dentate gyrus","Hippocampus","Hippocampus","Hypothalamus","Midbrain","Midbrain","Striatum"];% contains
catr=[["Basolateral amygdala";"."],["Piriform";"olfactory"],["Piriform";"pyramidal"],["";"."],["CA1";"CA3"],[".";"."],["";"."],[".";"."],[".";"."],["";"."]];% contains
catc=[".",".",".",".",".","TEINH",".","GLU","INH","."];% contains

cat2=['Cortex'];
cat2b=['pirifiorm'',olfactory'];
cat3=['Cortex'];
cat3b=['~pirifiorm'',Pyramidal']; % ~~~~~~~ for 3
cat4=['Dentate gyrus'];
cat4b=[''];
cat5=['Hippocampus'];
cat5b=['CA'];
cat6=['Hippocampus'];
cat6a=['TEINH'];
cat7=['Hypothalamus'];
cat7b=[''];
cat8=['Midbrain'];
cat8a=['GLU'];
cat9=['Midbrain'];
cat9a=['INH'];
cat10=['Striatum'];
cat10b=[''];
catcol=distinguishable_colors(length(catl));
%%  my subparallel plot

gf_names=grp(:,1);
rgn_names=grp(:,2);
loc_names=grp(:,3);
cor_val=str2double(grp(:,6));
m_names=grp(:,5);

rangx=length((cor_val));

for x=1:length(catl)

    subr=contains(rgn_names,catr(1,x),'IgnoreCase',true);
    subr2=contains(rgn_names,catr(2,x),'IgnoreCase',true);
    subl=contains(loc_names,catl(x),'IgnoreCase',true);
    subc=contains(m_names,catc(x),'IgnoreCase',true);
    subx=subl & (subr2 | subr | subc);

    if x==3 % for ~piriform
        subr=~contains(rgn_names,catr(1,x),'IgnoreCase',true);
        subr3=contains(rgn_names,'cingulate','IgnoreCase',true);

        subx=subl &  subr & (subr2 | subc | subr3);

    end

    gfx=gf_names(subx);
    [my,mi,mii]=unique(m_names(subx),'stable');
    mx=m_names(subx);
    rx=rgn_names(subx);
    lx=cor_val(subx);

    % set(h,'Rotation',45);
    dsize=100;
    fsize=10;
    % if length(gfx)>10
    % rangx=length(gfx);
    % else
    % rangx=10;
    % end
    spectrum=round(rangx.*(lx./max(cor_val)));
    cmapx=flip(turbo(rangx));%blue-red
    cmap=cmapx(spectrum,:);
    lwidth=flip(linspace(1,15,rangx));
    lwmap=lwidth(spectrum);
    % pos=length(gfx);
    [~,ci]=sort(lx);% sort color values

    % if length(gfx)>1
    % addstep=(length(gfx)-length(mi))/2 -1;%(length(n_gx)-sum(isnotzro))/2;
    % else
    addstep=0;
    % end
    % spacing=linspace(1,length(gfx),length(mi));
    spacefold=length(gfx)/length(mi);
    % subplot(length(catl),1,x)
    figure('color','w','units','normalized','outerposition',[0 0 1 1]);
    ax1 = axes('Position',[0.1300    0.1100    0.7750    0.8150]);

    axis off;axis tight;
    hold on
    currcol=gabaglut(subx,:);
    for i=1:length(lx)
        clw=lwmap(ci(i)) ;% current line width
        cc=cmap(ci(i),:);% current color
        id=mii(i);
        y=plot([-pos pos], [i id*(spacefold)+addstep],'MarkerEdgeColor','k','LineWidth',5,'color',cc);
        scatter(-pos,i,dsize,currcol(i,:),'filled','MarkerEdgeColor','k');%
        scatter(pos,id*(spacefold)+addstep,dsize,catcol(x,:),'s','filled','MarkerEdgeColor','k');
        text(-pos*1.5,i,gfx(i),'FontSize',fsize);
        text(pos*1.1,id*(spacefold)+addstep,[char(mx(i)),'-',char(rx(i)),'-',num2str(sum(mx(i)==mx))],'FontSize',10);
        %     dtRows = [dataTipTextRow("B-M:",[char(n_bx(i)),' = ',char(n_my(i))])];
        %     y.DataTipTemplate.DataTipRows(end+1) = dtRows;
    end
    scatter(-pos*1.5,1,'w');% just to force extend figure window
    scatter(pos*1.7,1,'w');% just to force extend figure window
    % s=scatter(-pos*ones(length(lx),1),(1:1:length(lx)),dsize,'r','filled');% add values to upper +20*ones(1,length(n_bx))
    % dtRows = [dataTipTextRow("GF:",gfx)];
    % s.DataTipTemplate.DataTipRows(end+1) = dtRows;
    % sx=scatter(pos*ones(length(my),1),(1:spacefold:length(lx))+addstep*ones(1,length(my)),dsize,'r','filled');
    % dtRows = [dataTipTextRow("M:",mx)];
    % sx.DataTipTemplate.DataTipRows(end+1) = dtRows;
    % n_by= regexprep(lx,'g-','');
    % text(-pos*1.5,length(gfx)+3,'Goldfish','FontSize',14)
    % text(pos*1.1,length(gfx)+3,'Mouse','FontSize',14)
    % text(0,length(gfx)+3,[char(catl(1,x)),'|',char(catr(1,x)),'|',char(catr(2,x)),'|',char(catc(1,x))],'FontSize',15)
    title([char(catl(1,x)),'|',char(catr(1,x)),'|',char(catr(2,x)),'|',char(catc(1,x))])

    % lxc=lx(ci);
    % [uclmp,~,iclmp]=unique(clmp,'rows','stable');

    % save2png(['/data/Technion_analysis/goldfish/scRNAseq_gf/comparative/n_mouse/',[num2str(x),'_gf2m']],gcf,300)
    if x==1
        colormap(flip(cmapx));
        colorbar('location','southoutside','Ticks', [linspace(0,1,length(cmapx))], 'TickLabels',[flip(sort(cor_val))]);
    end
    set(gcf, 'PaperPosition', [0.1300    0.1100    0.7750    0.8150]); %
    eval(['export_fig ',[num2str(x),'_gf2m'],'.pdf -r 600']);
    close all
end

%% colormap

% START of CUSTOM COLORBAR
L=max(idx); %  number of data points
indexValue = 0;     % value for which to set a particular color
topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and
% maximum values
largest = max(max(hm));
smallest = min(min(hm));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
    linspace(bottomcolor(2),indexColor(2),100*index)',...
    linspace(bottomcolor(3),indexColor(3),100*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
    linspace(indexColor(2),topColor(2),100*(L-index))',...
    linspace(indexColor(3),topColor(3),100*(L-index))'];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
customCMap=redblue(256);
cmpx=customCMap(round(length(customCMap)*spratio),:);
% colormap(customCMap)
%% plot centroid cell types
spratio=zeros(max(idx),1);% ratio between species
for u=1:max(idx)
    spratio(u)=sum(mgf(idx==u)==1)/length(mgf(idx==u));
end
customCMap=redblue(256);
cmpx=customCMap(round(length(customCMap)*spratio),:);
mxy=zeros(max(idx),2);
sxy=zeros(max(idx),1);
for xx=1:max(idx)
    mxy(xx,:)=median(mapped_xy(idx==xx,:));   % calculate the centorids by median | mean
    sxy(xx)=sum(idx==xx);% GET SUM OF SIZE OF EACH CLUSTER
end
% cmap1=distinguishable_colors(length(n_mx));
% cmap2=distinguishable_colors(length(n_gx));
cmap=clr;%distinguishable_colors(length(n_m));
nxy=5*rescale(sxy,5,1000);
figure('Color','w');
% ax1=subplot(1,2,1);
scatter(mapped_xy(:,1),mapped_xy(:,2),30,[0.5, 0.5, 0.5],'.'); % mouse

for ct=1:max(idx)
    hold on

    scatter(mxy(ct,1),mxy(ct,2),nxy(ct),cmpx(ct,:),'.')
end
axis off; axis equal; axis tight;
colormap(customCMap)
colorbar
% gscatter(mxy(1:96,1),mxy(1:96,2),n_mbx(1:96),cmap1(1:96,:),'.') % mouse
% hold on
% gscatter(mxy(97:end,1),mxy(97:end,2),n_mbx(97:end),cmap2(97:end,:),'+')   % goldfish
% legend off;
% axis off; axis equal; axis tight;
%%
figure('Color','w');
% ax2=subplot(1,2,2);
% scatter(mapped_xy(:,1),mapped_xy(:,2),1,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])

for ct=1:length(n_m)
    hold on

    scatter(mxy(ct,1),mxy(ct,2),nxy(ct),'.r')
end
for ct=(length(n_m)+1):length(n_mbx)
    hold on

    scatter(mxy(ct,1),mxy(ct,2),nxy(ct),'.b')
end
legend('Mouse')
axis off; axis equal; axis tight;
% linkaxes([ax1,ax2],'xy')
%% plot centroids species
% fc_time_uni=unique(all_flags(2).fc_time,'stable');
% mxy=zeros(length(fc_time_uni),2);
% sxy=zeros(length(fc_time_uni),1);
% for xx=1:130
%  mxy(xx,:)=median(mapped_xy((all_flags(1).fc_time==1),:));   % calculate the centorids by median | mean
%  sxy(xx)=sum((all_flags(1).fc_time==1));% GET SUM OF SIZE OF EACH CLUSTER
% end
% % bat cells
% for xx=131:length(fc_time_uni)
%  mxy(xx,:)=median(mapped_xy((all_flags(1).fc_time==2),:));   % calculate the centorids by median | mean
%  sxy(xx)=sum((all_flags(1).fc_time==2));% GET SUM OF SIZE OF EACH CLUSTER
% end
% cmap=distinguishable_colors(187);
% nxy=rescale(sxy,5,1000);
% figure('Color','w');
% axis off; axis equal; axis tight;
% for ct=1:187
%     hold on
%
% scatter(mxy(ct,1),mxy(ct,2),nxy(ct),cmap(ct,:),'.')
% end

%%
% figure('Color','w');
% scatter(mapped_xy(all_flags(1).fc_time==1,1),mapped_xy(all_flags(1).fc_time==1,2),'.r')
% hold on
% scatter(mapped_xy(all_flags(1).fc_time==2,1),mapped_xy(all_flags(1).fc_time==2,2),'.b')
% axis off; axis equal; axis tight;
% legend('Mouse', 'Bat')

% hold on
% scatter(mapped_xy(all_flags(1).fc_time==2,1),mapped_xy(all_flags(1).fc_time==2,2),'.b')
%
% hold on
% scatter(mxy(1:130,1),mxy(1:130,2),'*g')% mouse
% scatter(mxy(131:end,1),mxy(131:end,2),'+r')% bat
% for ix=1:187
% yc=mapped_xy(fc_time_all==ix,2);
% xc=mapped_xy(fc_time_all==ix,1);
% scbound = boundary(xc,yc,0.1);
% plot(xc(scbound),yc(scbound),'--','LineWidth',2,'Color','k');
% end
%% bat - mouse cell type
[M,I]=sort(Dx,1);
% M(2)
M_norm=1-M/max(max(M));
gmc=strings(length(n_gx),2);
for ix=1:length(n_gx)
    gmc(ix,1)=n_gx(ix);
    gmc(ix,2)=n_mx(I(1,ix));
    gmc(ix,3)=M(1,ix);
end


catbin = discretize(double(gmc(:,3)),10); % 10 is number of bins edges to create
[catbin,icatbin]=sort(catbin);
gmc=gmc(icatbin,:);
writematrix(gmc,"gmc.txt");
%% mouse - bat cell type
[M,I]=sort(Dx,2);
% M(2)
mgc=strings(length(n_mx),3);
for ix=1:length(n_mx)
    mgc(ix,1)=n_mx(ix);
    mgc(ix,2)=n_gx(I(ix,1));
    mgc(ix,3)=M(ix,1);
end
writematrix(mgc,"mgc.txt");
%%
[M,I]=sort(Dx,2);

D_knn=zeros(size(Dx));
for ix=1:length(n_mx)
    D_knn(ix,I(ix,[1]))=1;

end
% sd_knn=sort(D_knn,1,'rows');
[~,index] = sortrows(D_knn');
index=flipud(index(:));% in ecludian smaller means closer
figure; imagesc(D_knn(:,index))
hold on
set(gca, 'YTick',numm, 'YTickLabel',n_mx)
set(gca, 'XTick',numb, 'XTickLabel',n_gx(index))
xtickangle(45)
% xline(1,'-r',{'Mouse','GABA'},'LineWidth',2)
% xline(57,'-r',{'Mouse','Glut2'},'LineWidth',2)
% xline(89,'-r',{'Mouse','Glut1'},'LineWidth',2)
% xline(131,'-r',{'Bat','GABA'},'LineWidth',2)
% xline(161,'-r',{'Bat','Glut2'},'LineWidth',2)
% xline(179,'-r',{'Bat','Glut1'},'LineWidth',2)
% yline(1,'-r',{'Mouse','GABA'},'LineWidth',2)
% yline(57,'-r',{'Mouse','Glut2'},'LineWidth',2)
% yline(89,'-r',{'Mouse','Glut1'},'LineWidth',2)
% yline(131,'-r',{'Bat','GABA'},'LineWidth',2)
% yline(161,'-r',{'Bat','Glut2'},'LineWidth',2)
%% redblue correlation bat mouse heatmap
figure('Color','w');
cmap=redblue(256);
colormap(flip(cmap,2))
Dsort=Dz(:,index);
Yup = prctile(Dsort(:),95);% thresholding imagesc
Ydn = prctile(Dsort(:),1);
Dsortx=Dsort';
Dsortx(:,55)=[];
% imagesc(Dsortx,[Ydn Yup])
imagesc(Dz',[Ydn Yup])

ylabel('Goldfish')
xlabel('Mouse')
% ylabel('Bat')
% xlabel('Mouse')

numm=1:numel(n_mx);
numb=1:numel(n_gx);
set(gca, 'XTick',numm,'fontsize',5,'XTickLabel',n_mx)
set(gca, 'YTick',numb, 'YTickLabel',n_gx)
xtickangle(45)
% xline(1,'-k',{'Mouse','GABA'},'LineWidth',2)
% xline(57,'-k',{'Mouse','Glut2'},'LineWidth',2)
% xline(89,'-k',{'Mouse','Glut1'},'LineWidth',2)
%
% yline(1,'-k',{'Bat','GABA'},'LineWidth',2)
% yline(31,'-k',{'Bat','Glut2'},'LineWidth',2)
% yline(49,'-k',{'Bat','Glut1'},'LineWidth',2)

%%
[M,I]=sort(Dx,1);
D_knn=zeros(size(Dx));
for ix=1:length(n_gx)
    D_knn(I(ix,[1]),ix)=1;

end
% sd_knn=sort(D_knn,1,'rows');
[~,index] = sortrows(D_knn);
index=flipud(index(:));
figure; imagesc(D_knn(index,:))
hold on
set(gca, 'YTick',numm, 'YTickLabel',n_mx(index))
set(gca, 'XTick',numb, 'XTickLabel',n_gx)
xtickangle(45)
% xline(1,'-r',{'Mouse','GABA'},'LineWidth',2)
% xline(57,'-r',{'Mouse','Glut2'},'LineWidth',2)
% xline(89,'-r',{'Mouse','Glut1'},'LineWidth',2)
% xline(131,'-r',{'Bat','GABA'},'LineWidth',2)
% xline(161,'-r',{'Bat','Glut2'},'LineWidth',2)
% xline(179,'-r',{'Bat','Glut1'},'LineWidth',2)
% yline(1,'-r',{'Mouse','GABA'},'LineWidth',2)
% yline(57,'-r',{'Mouse','Glut2'},'LineWidth',2)
% yline(89,'-r',{'Mouse','Glut1'},'LineWidth',2)
% yline(131,'-r',{'Bat','GABA'},'LineWidth',2)
% yline(161,'-r',{'Bat','Glut2'},'LineWidth',2)
% yline(179,'-r',{'Bat','Glut1'},'LineWidth',2)
%% bat 2 mouse
[M,I]=sort(Dx,1);
D_knn=zeros(size(Dx));
for ix=1:length(n_g)
    D_knn(I(1,ix),ix)=1;
end
% sd_knn=sort(D_knn,1,'rows');
[~,index] = sortrows(D_knn);
index=flipud(index(:));
figure; imagesc(D_knn(index,:))
hold on
set(gca, 'YTick',numm, 'YTickLabel',n_mx(index))
set(gca, 'XTick',numb, 'XTickLabel',n_gx)
xtickangle(45)

%% all polarplot mbc
% figure;
% set(gcf,'color','w');
% num_tiks=length(n_bx);
% deg=360/num_tiks; % find deg ;
% rdeg=2*pi/num_tiks; % find deg in radians
% tick_v=deg:deg:360;
% theta = rdeg:rdeg:2*pi;
% % Isix=I(:,1:6);
% Isix=1-(Dx/max(max(Dx)));
% [p,~]=numSubplots(num_tiks);
% for isp=1:length(n_mx)
%     subplot(p(1),p(2), isp);
%     % rho=mbc_v(isp,[2:end]);
%     rho=Isix(isp,:);
%     polarplot(theta,rho,'-o')
%     thetaticks([tick_v])
%     thetaticklabels(n_bx);
%     subtitle(n_mx(isp))
% end
% % title('mouse to bat celll type comparison')
% %% sub- polarplot mbc : according to numticks
% figure;
% set(gcf,'color','w');
% num_tiks=6; %length(n_bx);
% deg=360/num_tiks; % find deg ;
% rdeg=2*pi/num_tiks; % find deg in radians
% tick_v=deg:deg:360;
% theta = rdeg:rdeg:2*pi;
% M_norm=1-M/max(max(M));
% Isix=M_norm(:,1:num_tiks);
%
% % Isix=1-(Dx/max(max(Dx)));
% % [p,~]=numSubplots(length(n_mx));
% [p,~]=numSubplots(6);
% for isp=1:6%length(n_mx)
%     subplot(p(1),p(2), isp);
%     % rho=mbc_v(isp,[2:end]);
%     rho=Isix(isp,:);
%     polarplot(theta,rho,'-o')
%     thetaticks([tick_v])
%     thetaticklabels(n_bx(I(isp,1:6)));
%     subtitle(n_mx(isp))
%
% end
% % title('mouse to bat celll type comparison')
%% parallel plot  sort(Isix,'descend')
% create  matrix that shows the two high cell types
% t_bm=array2table(Isix');

t_bm= table(gmc(:,1),gmc(:,2),catbin,'VariableNames',{'Bat','Mouse','strength'});
figure('Units','normalized','Position',[0.3 0.3 0.45 0.4],'Color','w');
coordvars=[1,2];
p = parallelplot(t_bm,'CoordinateVariables',coordvars,'GroupVariable','strength');
p.CoordinateTickLabels = {'Bat','Mouse'};
p.Color = hot(10);
%%  my parallel plot
isnotzro = any(D_knn(index,:),2); % check zero rows in D_knn
D_my=D_knn(index,:);
n_my=n_mx(index);
n_my=n_my(isnotzro);
figure('color','w');
%  set(h,'Rotation',-45);
hold on
% set(h,'Rotation',45);

text(-20,-5,'Goldfish')
text(12,-5,'Mouse')
cmap=flip(bone(length(n_g)));%blue-red
lwmap=linspace(0.1,3,length(n_g));
axis off; axis tight; axis equal;

[cv,ci]=sort(M(1,:));% sort color values
addstep=0;%(length(n_gx)-sum(isnotzro))/2;
spacing=linspace(1,length(n_gx),sum(isnotzro));
spacefold=length(n_gx)/sum(isnotzro);
for i=1:length(n_gx)
    clw=lwmap(ci(i)) ;% current line width
    cc=cmap(ci(i),:);% current color
    id=find(D_my(:,i)==1); %what value in imagesc D_my is one
    y=plot([-10 10], [i id*(spacefold)-1+addstep],'k','LineWidth',clw,'color',cc);
    %     dtRows = [dataTipTextRow("B-M:",[char(n_bx(i)),' = ',char(n_my(i))])];
    %     y.DataTipTemplate.DataTipRows(end+1) = dtRows;
end
s=scatter(-10*ones(length(n_gx),1),(1:1:length(n_gx)),30,'r','filled');% add values to upper +20*ones(1,length(n_bx))
dtRows = [dataTipTextRow("GF:",n_gx)];
s.DataTipTemplate.DataTipRows(end+1) = dtRows;
sx=scatter(10*ones(length(n_my),1),(1:spacefold:length(n_gx))+addstep*ones(1,length(n_my)),30,'r','filled');
dtRows = [dataTipTextRow("M:",n_my)];
sx.DataTipTemplate.DataTipRows(end+1) = dtRows;
n_by= regexprep(n_gx,'g-','');

text(-20*ones(length(n_gx),1),(1:1:length(n_gx)),n_by,'FontSize',5);
text(12*ones(length(n_my),1),(spacing)+addstep*ones(1,length(n_my)),n_my,'FontSize',5);
colormap(cmap)
title('Goldfish to Mouse')
%%
Dtot=logical(eye(length(u_c)));
for xx=1:size(Dtot,1)
    %     [~,Ix]=sort(u_c,2);
    [mv,mi]=min(u_c(xx,1:130));
    if mv==0
        [mv,mi]=sort(u_c(xx,1:130));
        mi=mi(2);
    end
    [bv,bi]=min(u_c(xx,131:187));
    if bv==0
        [bv,bi]=sort(u_c(xx,131:187));
        bi=bi(2);
    end
    Dtot(xx,bi+130)=1;
    Dtot(xx,mi)=1;
end
Dtot = Dtot+Dtot';
Dtot(Dtot==2)=1;
figure;imagesc(Dtot)
xline(131,'-r',{'M|B'},'LineWidth',2)
yline(131,'-r',{'M|B'},'LineWidth',2)

% axis equal
%% tsne_distribution
labels=[];
no_dims=2;
perplexity=5;
mappedX = tsne_d(Dtot,[], no_dims, perplexity);
sizes=1:187;


colors=turbo(length(Dtot));
figure('color','w');axis equal; axis off;
scatter(mappedX(:,1),mappedX(:,2),sizes,colors,'filled');
%% gaba vglut1 vglut2 graph plot
G = graph(Dtot,'omitselfloops');

colors=zeros(length(Dtot),3);
% bat
bgaba=contains(n_mbx,'g-GABA');
colors(bgaba,:)=ones(sum(bgaba),3).*([0,0,1]);
bglut1=contains(n_mbx,'g-Glut1');
colors(bglut1,:)=ones(sum(bglut1),3).*([1,0,0]);
bglut2=contains(n_mbx,'g-Glut2');
colors(bglut2,:)=ones(sum(bglut2),3).*([0,1,0]);
% mouse
bgaba=contains(n_mbx(1:130),'GABA');
colors(bgaba,:)=ones(sum(bgaba),3).*([0,0,1]);
bglut1=contains(n_mbx(1:130),'Glut1');
colors(bglut1,:)=ones(sum(bglut1),3).*([1,0,0]);
bglut2=contains(n_mbx(1:130),'Glut2');
colors(bglut2,:)=ones(sum(bglut2),3).*([0,1,0]);
% plot graph
figure('color','w');
% p=plot(G,'NodeLabel',n_mbx,'Layout','force','NodeColor','k');

p=plot(G,'Layout','force','NodeColor',colors);
axis equal; axis off;
% p.Marker = 's';
highlight(p,1:187,'MarkerSize',5)
p.NodeFontSize = 5;
p.NodeFontAngle='italic';
h = zeros(3, 1);
hold on
h(1) = plot(NaN,NaN,'.b');
h(2) = plot(NaN,NaN,'.r');
h(3) = plot(NaN,NaN,'.g');
hleg=legend('','GABA','Glut1','Glut2');
% hleg.String(1) = []; % delete the last legend entry of the very last plot
%% bat-mouse: gaba vglut1 vglut2 graph plot

colors=zeros(length(Dtot),3);
% mouse
bgaba=contains(n_mbx(1:130),'GABA');
colors(bgaba,:)=ones(sum(bgaba),3).*([0,0,1]);
bglut1=contains(n_mbx(1:130),'Glut1');
colors(bglut1,:)=ones(sum(bglut1),3).*([1,0,0]);
bglut2=contains(n_mbx(1:130),'Glut2');
colors(bglut2,:)=ones(sum(bglut2),3).*([0,1,0]);
% bat
bgaba=contains(n_mbx,'b-GABA');
colors(bgaba,:)=ones(sum(bgaba),3).*([0,0,0.1724]);
bglut1=contains(n_mbx,'b-Glut1');
colors(bglut1,:)=ones(sum(bglut1),3).*([1.0000   , 0.1034  ,  0.7241]);
bglut2=contains(n_mbx,'b-Glut2');
colors(bglut2,:)=ones(sum(bglut2),3).*([1.0000   , 0.8276     ,    0]);

% plot graph
figure('color','w');
p=plot(G,'Layout','force','NodeColor',colors);
axis equal; axis off;
% p.Marker = 's';
highlight(p,1:187,'MarkerSize',7)
p.NodeFontSize = 5;
h = zeros(3, 1);
hold on
h(1) = scatter(NaN,NaN,'MarkerFaceColor',[0,0,0.1724],'MarkerEdgeColor',[1,1,1]);
h(2) = scatter(NaN,NaN,'MarkerFaceColor',[1.0000   , 0.1034  ,  0.7241],'MarkerEdgeColor',[1,1,1]);
h(3) = scatter(NaN,NaN,'MarkerFaceColor',[1.0000   , 0.8276     ,    0],'MarkerEdgeColor',[1,1,1]);
h(4) = scatter(NaN,NaN,'MarkerFaceColor','b','MarkerEdgeColor',[1,1,1]);
h(5) = scatter(NaN,NaN,'MarkerFaceColor','r','MarkerEdgeColor',[1,1,1]);
h(6) = scatter(NaN,NaN,'MarkerFaceColor','g','MarkerEdgeColor',[1,1,1]);
hleg=legend('','b-GABA','b-Glut1','b-Glut2','m-GABA','m-Glut1','m-Glut2');

%% bat /mouse graph plot
colors=zeros(length(Dtot),3);

bgaba=logical([zeros(130,1);ones(57,1)]);%contains(n_mbx,'b-');% bat
colors(bgaba,:)=ones(sum(bgaba),3).*([0,0,1]);
colors(~bgaba,:)=ones(sum(~bgaba),3).*([1,0,0]);

% % bglut1=contains(n_mbx,'b-Glut1');
% % colors(bglut1,:)=ones(sum(bglut1),3).*([0,1,0]);
% % bglut2=contains(n_mbx,'b-Glut2');
% % colors(bglut2,:)=ones(sum(bglut2),3).*([0,0,1]);
% % % mouse
% % bgaba=contains(n_mbx(1:130),'GABA');
% % colors(bgaba,:)=ones(sum(bgaba),3).*([1,1,0]);
% % bglut1=contains(n_mbx(1:130),'Glut1');
% % colors(bglut1,:)=ones(sum(bglut1),3).*([0,1,1]);
% % bglut2=contains(n_mbx(1:130),'Glut2');
% % colors(bglut2,:)=ones(sum(bglut2),3).*([1,0,1]);
% plot graph
figure('color','w');
p=plot(G,'Layout','force','NodeColor',colors);
axis equal; axis off;
% p.Marker = 's';
highlight(p,1:187,'MarkerSize',7)
p.NodeFontSize = 5;
h = zeros(3, 1);
hold on
h(1) = scatter(NaN,NaN,'MarkerFaceColor','r','MarkerEdgeColor',[1,1,1]);
h(2) =scatter(NaN,NaN,'MarkerFaceColor','b','MarkerEdgeColor',[1,1,1]);
hleg=legend('','Mouse','bat');

%% cell types graph plot
colors=zeros(length(Dtot),3);
% numarr=string(1:187);
bgaba=logical([zeros(130,1);ones(57,1)]);%contains(n_mbx,'b-');% bat
% colors(bgaba,:)=ones(sum(bgaba),3).*([0,0,1]);
% colors(~bgaba,:)=ones(sum(~bgaba),3).*([1,0,0]);

% % bglut1=contains(n_mbx,'b-Glut1');
% % colors(bglut1,:)=ones(sum(bglut1),3).*([0,1,0]);
% % bglut2=contains(n_mbx,'b-Glut2');
% % colors(bglut2,:)=ones(sum(bglut2),3).*([0,0,1]);
% % % mouse
% % bgaba=contains(n_mbx(1:130),'GABA');
% % colors(bgaba,:)=ones(sum(bgaba),3).*([1,1,0]);
% % bglut1=contains(n_mbx(1:130),'Glut1');
% % colors(bglut1,:)=ones(sum(bglut1),3).*([0,1,1]);
% % bglut2=contains(n_mbx(1:130),'Glut2');
% % colors(bglut2,:)=ones(sum(bglut2),3).*([1,0,1]);
% plot graph
figure('color','w');
p=plot(G,'Layout','force','NodeColor',colors,'Nodelabel',n_mbx);
axis equal; axis off;
% p.Marker = 's';
highlight(p,1:187,'MarkerSize',2)
p.NodeFontSize = 6;
% h = zeros(3, 1);
% hold on
% h(1) = scatter(NaN,NaN,'MarkerFaceColor','r','MarkerEdgeColor',[1,1,1]);
% h(2) =scatter(NaN,NaN,'MarkerFaceColor','b','MarkerEdgeColor',[1,1,1]);
% hleg=legend('','Mouse','bat');
%%

% [Mx,Ix]=sort(u_c(1:130,1:130),2);
% Im=Ix(:,2);
%
% [Mx,Ix]=sort(u_c(131:187,131:187),2);
% Ib=Ix(:,2);
%
% [Mx,Ix]=sort(u_c(1:130,131:187),2);
% Imb=Ix(:,1);
%
% [Mx,Ix]=sort(u_c(131:187,1:130),2);
% Ibm=Ix(:,1);
%