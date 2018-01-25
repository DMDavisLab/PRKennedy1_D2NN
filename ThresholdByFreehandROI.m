
%Find the folder containing the already processed data

DataDirName = uigetdir(pwd,'Choose your data folder');
    
    if DataDirName ~= 0
        
        cd(DataDirName);

        % open the file containing procsettings
        if exist(fullfile(cd, 'ProcSettings.txt'), 'file') == 0
            errordlg('Cannot find ''ProcSettings.txt'' file for your data. You''ll need one to proceed.');
        else
            ProcSettings = LoadProcSettings;
        end
        
%         if isdir('Events')
%             cd('Events');

            DirFullList = dir;                                 % Get the data for the current directory
            FolderIndex = [DirFullList.isdir];                 % Find the index for directories
            FileList = {DirFullList(~FolderIndex).name}';      % Get a list of the files

            for f = 1:length(FileList)
                [~, fname, ext] = fileparts(FileList{f,1});
                FileList{f,2} = fname;
                FileList{f,3} = ext;
            end

            GoodIdx = strcmp(['.',ProcSettings.DataTableFileExt],FileList(:,3));
            FileList(~GoodIdx,:) = [];
            [~, natsortidx] = sort_nat(FileList(:,2),'Ascend');
            FileList = FileList(natsortidx,:);
            
%             cd(DataDirName); % back to the main folder
%         else
%             errordlg('Cannot find the ''Events'' folder for your data.');
%         end

    else
        error('Cancelled?! So rude.');
    end
       
    clear DirFullList FolderIndex f fname ext GoodIdx InfoMessage

% for p = 1:size(FileList,1)
%     
%     %load the image
%     events_img_data = imread(fullfile(pwd,'Events',FileList{p,1}));
%     events_img_fig = imshow(events_img_data);
% %     get(events_img_fig); %, 'Position', get(0,'Screensize')); % Maximize figure.
%     h_evimg_roi = imfreehand(gca);
%     evimg_roi_coords = h_evimg_roi.getPosition();
% 
% 
%     
% end
    
    

% requires data_proc ProcSettings
% data_proc = the processed data, e.g. from the MAT file
ProcSettings.xCol_procd = 1;
ProcSettings.yCol_procd = 2;
ProcSettings.AxisLimits = [0 ProcSettings.ImageSize 0 ProcSettings.ImageSize];
ProcSettings.invlogAFDCol = 8;
ProcSettings.PlotWidth = ProcSettings.AxisLimits(2) - ProcSettings.AxisLimits(1);
ProcSettings.SaveHighDPI = strcat('-r',num2str(ProcSettings.PlotWidth / 10)); % for 10 nm per pixel
ProcSettings.SaveLowDPI = '-r300'; % for quicker saving of preview images
stats_manual = zeros(4,size(FileList,1));      % init stats holder
                                                % row 1 - threshold value for that table's method
                                                % row 2 - total points over threshold
                                                % row 3 - total points in image
                                                % row 4 - percent over threshold

for FileID = 1:size(FileList,1)
    
    % load the mat file
    load(fullfile(pwd,'/Z - MAT Files',[num2str(FileID),'.mat']),'data_proc');
                                                
    data_len = length(data_proc(:,1));     %Get the length of the input vectors.
    data_rand_ind=randperm(data_len);     %Create a vector of random integers from 1 to len
    data_randx=data_proc(data_rand_ind,ProcSettings.xCol_procd);  %Randomize the input vectors.
    data_randy=data_proc(data_rand_ind,ProcSettings.yCol_procd);  %Randomize the input vectors.

    fig_plot_pts_roi = figure('Color','k', ...
                            'Renderer', 'OpenGL', ...
                            'Units', 'inches', ...
                            'PaperUnits', 'inches', ...
                            'PaperSize', [10 10], ...
                            'PaperPositionMode', 'manual', ...
                            'PaperPosition', [0 0 10 10], ...
                            'Visible','on');
    DataEvents_axes =  axes('DataAspectRatio', [1,1,1], ...
                            'Position', [0 0 1 1], ...
                            'Color','k',...
                            'Visible','off'); 
    axis(ProcSettings.AxisLimits);
    axis square image tight
    box('off');
    hold on % don't reset anything when plotting this next thing...

    scatter(data_randx,data_randy,2,[1 0.02 0.2],'.');


    % make a temporary new set of axes to trick matlab in being faster while
    % drawing

    xl = get(gca,'xlim'); 
    yl = get(gca,'ylim'); 
    pos = get(gca,'pos'); 
    xc = get(gca,'xcolor'); 
    yc = get(gca,'ycolor'); 
    xsc = get(gca,'xscale'); 
    ysc = get(gca,'yscale'); 
    xdir = get(gca,'xdir'); 
    ydir = get(gca,'ydir');

    tmpax = axes('pos',pos,...
        'xlim',xl,'ylim',yl,...
        'color','none',...
        'xcolor',xc,'ycolor',yc,...
        'xtick',[],'xticklabel','',...
        'ytick',[],'yticklabel','',...
        'xdir',xdir,'ydir',ydir,...
        'xscale',xsc,'yscale',ysc); 
    axis equal
    axis(ProcSettings.AxisLimits);

    fcn = makeConstrainToRectFcn('imfreehand',get(gca,'XLim'),get(gca,'YLim'));
    hFH = imfreehand(gca,'Closed','true','PositionConstraintFcn',fcn);
    % api.setPositionConstraintFcn(fcn);

    %get coords of imfreehand directly
    pos = hFH.getPosition();

    % Delete temporary set of axes: 
    delete(tmpax)

    %plot the ROI again (to check)
    plot(pos(:,1),pos(:,2),'b','LineWidth',2)
    SaveFileName = [FileList{FileID,2},' - ROI Preview.png'];
    print(fig_plot_pts_roi,'-dpng',ProcSettings.SaveLowDPI,fullfile(pwd,SaveFileName));
    close(fig_plot_pts_roi);
    
    %save the roi
    fnsave = [FileList{FileID,2},' - ROI.mat'];
    save(fullfile(pwd,fnsave),'pos');

    % get events within this region
    [in_roi,on_roi] = inpolygon(data_proc(:,1),data_proc(:,2),pos(:,1),pos(:,2));
%     figure
%     scatter(data_proc(in_roi,1),data_proc(in_roi,2),1,'.');
%     axis equal
%     axis(ProcSettings.AxisLimits);
%     yt = get(gca,'YTick');
%     yt = yt/1000;
%     set(gca,'YTickLabel', sprintf('%d|',yt));
%     xt = get(gca,'XTick');
%     xt = xt/1000;
%     set(gca,'XTickLabel', sprintf('%d|',xt));
%     set(gca,'TickDir','out');

% plot the histogram within the cell

    invlogsumdistanceplot = [];
    xbins = floor(min(data_proc(in_roi,ProcSettings.invlogAFDCol))):0.1:ceil(max(data_proc(in_roi,ProcSettings.invlogAFDCol))); % bin width = 10 (nm)
    [invlogsumdistanceplot(2,:),invlogsumdistanceplot(1,:)] = hist(data_proc(in_roi,ProcSettings.invlogAFDCol),xbins);
    fig_ROI_Hist = figure;
    bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:));
    axis tight
    set(gca,'TickDir','out');
    xlabel(['Inv. Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
    ylabel('Freq.');
    title(['Max Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})-Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
    normallogsumdistance2 = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
    axis([floor(min(data_proc(in_roi,ProcSettings.invlogAFDCol))) ceil(max(data_proc(in_roi,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == normallogsumdistance2)]);
    SaveFileName = [FileList{FileID,2},' - ROI Histogram.png'];
    print(fig_ROI_Hist,'-dpng',ProcSettings.SaveLowDPI,fullfile(pwd,SaveFileName));
    close(fig_ROI_Hist);


% Do a threshold
    fprintf('\t%s','Adaptive threshold > ROI Histogram max');

    Cluster_thr_histmax = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,:)));
%         Cluster_thr_histmax = 3.0;
    SaveToFolder = '';
    SaveFileName = [num2str(FileID),' - Threshold ROI Histogram peak(',num2str(Cluster_thr_histmax),')'];
    RefDataCol = ProcSettings.invlogAFDCol;

    [histmax_pts_over, histmax_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_histmax,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
    num_events_over = size(histmax_pts_over,1);

%     fig_plot_pts = figure;
%     scatter(histmax_pts_under(:,1),histmax_pts_under(:,2),1,'b','.');
%     hold on
%     scatter(histmax_pts_over(:,1),histmax_pts_over(:,2),1,'r','.');
%     axis([min(data_proc(in_roi,ProcSettings.xCol_procd)) max(data_proc(in_roi,ProcSettings.xCol_procd)) min(data_proc(in_roi,ProcSettings.yCol_procd)) max(data_proc(in_roi,ProcSettings.yCol_procd))]);
%     axis square
%     print(fig_plot_pts,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFileName));
%     close(fig_plot_pts);


    % update the summary table
    stats_manual(1,FileID) = Cluster_thr_histmax;           % row 1 - threshold value for that table's method
    stats_manual(2,FileID) = num_events_over;               % row 2 - total points over threshold
    stats_manual(3,FileID) = sum(in_roi);   % row 3 - total points in ROI
    stats_manual(4,FileID) = (num_events_over/sum(in_roi)) * 100;    % row 4 - percent over threshold

    fprintf('\t%s','(done)');
    fprintf('\n');

end

    % format for headers (first row)
    stringy = '%s';
    for g = 1:(size(FileList,1)-1)
        stringy = strcat(stringy,'\t%s'); % replace \t with a comma for csv
    end
    
    % Save ROI stats table
    fn_summary_proc_save = 'RoI Thresholds.txt';
    fid = fopen(fullfile(pwd,fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - by ROIs');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,fn_summary_proc_save),stats_manual,'-append','delimiter','\t','precision','%.3f');    

