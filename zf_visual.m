%% visual for bat clean after (the newest)
close all
clear all
clc
cond=3;

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
%% scatter tsne with colormap
figure('color','w');



all_medianxy=zeros(length(ut),2);
hold on
%
for x=1:length(ut)
    curr_xy=find(T_cells_tmp==ut(x));
    all_medianxy(x,:)=median([mapped_xy(curr_xy,:)]);
    scatter(mapped_xy(curr_xy,1),mapped_xy(curr_xy,2),10,cmap(x,:),'filled');

end
lgd=legend;
textPosition=calcTextPosInScatterPlot(all_medianxy);
for x=1:length(ut)
    curr_xy=find(T_cells_tmp==ut(x));
    plot([textPosition(x,1) all_medianxy(x,1)], ...
        [textPosition(x,2) all_medianxy(x,2)],'k--','LineWidth',2);
    text(textPosition(x,1), textPosition(x,2), ...
        u_cx(x), 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');%, ...        % 'BackgroundColor',[0.8 0.8 0.8]
    % %          num2str(x) or u_cx(x,:), 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');%, ...        % 'BackgroundColor',[0.8 0.8 0.8]

end


axis off; axis equal;
%% same but with numbers
figure('color','w');



all_medianxy=zeros(length(ut),2);
hold on
%
for x=1:length(ut)
    curr_xy=find(T_cells_tmp==ut(x));
    all_medianxy(x,:)=median([mapped_xy(curr_xy,:)]);
    scatter(mapped_xy(curr_xy,1),mapped_xy(curr_xy,2),10,cmap(x,:),'filled');

end
lgd=legend;
textPosition=calcTextPosInScatterPlot(all_medianxy);
% for x=1:length(ut)
%     curr_xy=find(T_cells_tmp==ut(x));
%     plot([textPosition(x,1) all_medianxy(x,1)], ...
%         [textPosition(x,2) all_medianxy(x,2)],'k--','LineWidth',2);
%     text(textPosition(x,1), textPosition(x,2), ...
%           num2str(x), 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');%, ...        % 'BackgroundColor',[0.8 0.8 0.8]
%
% end

%% show specific gene
% NPT: c/uc
% nupt=table2array(readtable('/data/Technion_analysis/zebrafish/sc_100410/np_tf_list.xlsx'));
% gname=unique(upper(string([nupt(:,1);nupt(:,2);(nupt(:,3))])),'stable');% change to add genes 
% gname(cellfun('isempty',gname)) = []; % remove empty strings
% gname=[gname;"ADCYAP1B";"SLC17A6A";"SLC32A1";"GAD2"];
% TF: c/uc 
genet=readtable('/data/Technion_analysis/zebrafish/sc_100410/TFs.xlsx');
gname=upper(string(table2array(genet(:,1))));
gname(cellfun('isempty',gname)) = []; % remove empty strings
gname=[gname;"FOXP4";"LHX6"];
%%
for g=1:length(gname)

    figure;
    set(gcf,'Color','w')
    fg=gname(g);
    gx=find(string(geneid_all)==fg);
    if isempty(gx)
        disp(fg)
        disp("NOT FOUND!")
    else
        markergene=double(n_data(gx,:)>0);
        % tmpthlow = prctile(markergene(markergene>0),30);
        % tmpthhigh =prctile(markergene(markergene>0),80);
        % markergene(markergene>tmpthhigh) = tmpthhigh;
        % markergene(markergene<tmpthlow) = tmpthlow;
        c_rgb = [1,0,0];
        try
            markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];

            scatter(mapped_xy(:,1),mapped_xy(:,2),30,markergene_color,'.'); % mouse
        end
        axis off; axis tight; axis equal;
        title(fg)
    end
namepdf=char(gname(g));
eval(['export_fig ',namepdf,'.pdf  -nocrop -r 300']);
close all
end