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
    
% requires data_proc ProcSettings
% data_proc = the processed data, e.g. from the MAT file
ProcSettings.xCol_procd = 1;
ProcSettings.yCol_procd = 2;
ProcSettings.AxisLimits = [0 ProcSettings.ImageSize 0 ProcSettings.ImageSize];
ProcSettings.invlogAFDCol = 8;
ProcSettings.PlotWidth = ProcSettings.AxisLimits(2) - ProcSettings.AxisLimits(1);
ProcSettings.SaveHighDPI = strcat('-r',num2str(ProcSettings.PlotWidth / 10)); % for 10 nm per pixel
ProcSettings.SaveLowDPI = '-r300'; % for quicker saving of preview images

DoThrByHHMX = true;
DoThrByHMX = true;
DoThrByMMSD = true;
DoThrByMPSD = true;
DoThrByMEDN = true;
DoThrByFIXED = true;

% create containers for threshold summaries
    stats_hhmx = zeros(4,size(FileList,1));     % row 1 - threshold value for that table's method
    stats_mmsd = zeros(4,size(FileList,1));     % row 2 - total points over threshold
    stats_medn = zeros(4,size(FileList,1));     % row 3 - total points in image
    stats_hmx = zeros(4,size(FileList,1));      % row 4 - percent over threshold
    stats_mpsd = zeros(4,size(FileList,1));
    stats_thrfix = zeros(4,size(FileList,1));   % fixed-threshold stats

for FileID = 1:size(FileList,1)
    
    % load the mat file
    load(fullfile(pwd,'/Z - MAT Files',[num2str(FileID),'.mat']),'data_proc');
    load(fullfile(pwd,[num2str(FileID),' - Events - ROI.mat']),'pos');
    
    % get events within this region
    [in_roi,on_roi] = inpolygon(data_proc(:,1),data_proc(:,2),pos(:,1),pos(:,2));
    
    
    % plot the histogram within the cell
    invlogsumdistanceplot = [];
    xbins = floor(min(data_proc(in_roi,ProcSettings.invlogAFDCol))):0.1:ceil(max(data_proc(in_roi,ProcSettings.invlogAFDCol))); % bin width = 10 (nm)
    [invlogsumdistanceplot(2,:),invlogsumdistanceplot(1,:)] = hist(data_proc(in_roi,ProcSettings.invlogAFDCol),xbins);
%     fig_ROI_Hist = figure;
%     bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:));
%     axis tight
%     set(gca,'TickDir','out');
%     xlabel(['Inv. Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
%     ylabel('Freq.');
%     title(['Max Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})-Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
%     normallogsumdistance2 = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
%     axis([floor(min(data_proc(in_roi,ProcSettings.invlogAFDCol))) ceil(max(data_proc(in_roi,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == normallogsumdistance2)]);
% %     SaveFileName = [FileList{FileID,2},' - ROI Histogram.png'];
% %     print(fig_ROI_Hist,'-dpng',ProcSettings.SaveLowDPI,fullfile(pwd,SaveFileName));
%     close(fig_ROI_Hist);
% 
% 
% % Do a threshold
%     fprintf('\t%s','Adaptive threshold > ROI Histogram max');
%     Cluster_thr_histmax = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,:)));
% %   Cluster_thr_fixed = ProcSettings.FixedThreshold;; % fixed threshold
%     SaveToFolder = '';
%     SaveFileName = [num2str(FileID),' - Threshold ROI Histogram peak(',num2str(Cluster_thr_histmax),')'];
%     RefDataCol = ProcSettings.invlogAFDCol;
%     [histmax_pts_over, histmax_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_histmax,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
%     num_events_over = size(histmax_pts_over,1);
% 
% 
%     % update the summary table
%     stats_manual(1,FileID) = Cluster_thr_histmax;           % row 1 - threshold value for that table's method
%     stats_manual(2,FileID) = num_events_over;               % row 2 - total points over threshold
%     stats_manual(3,FileID) = ;   % row 3 - total points in ROI
%     stats_manual(4,FileID) = (num_events_over/sum(in_roi)) * 100;    % row 4 - percent over threshold
% 
%     fprintf('\t%s','(done)');
%     fprintf('\n');

%% Thresholding

        total_events_table = sum(in_roi);
        
        % per-cell data table information
        HeaderFormat = '%s\t%s\t%s';
        fn_stats_save = [FileList{FileID,2},' - Threshold Data (within ROI).txt'];
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_stats_save),'w');
        fprintf(fid,HeaderFormat,'[Threshold Method]','[Threshold Value]','[% In Clusters]');
        fprintf(fid,'\r\n');
        DataFormat = '%s\t%.3f\t%.3f';
                
%====== Adaptive thresholding - Half hist max
    if DoThrByHHMX
            fprintf('%s', [FileList{FileID,2},' - Applying Thresholds:']);
            fprintf('\n');
            fprintf('\t%s','Adaptive threshold > Half histogram max');

            Cluster_thr_halfhistmax = (invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,:))))/2;
            SaveToFolder = fullfile('D - Threshold','Half histogram peak');
            SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_halfhistmax),' (Half hist peak - in ROI)'];
            RefDataCol = ProcSettings.invlogAFDCol;

            [halfhistmax_pts_over, halfhistmax_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_halfhistmax,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
            num_events_over = size(halfhistmax_pts_over,1);
       
            % update the summary table
            stats_hhmx(1,FileID) = Cluster_thr_halfhistmax;
            stats_hhmx(2,FileID) = num_events_over;
            stats_hhmx(3,FileID) = total_events_table;
            stats_hhmx(4,FileID) = (num_events_over/total_events_table) * 100;

            % save results for this cell to a file
            fprintf(fid,DataFormat,'Half hist. peak',Cluster_thr_halfhistmax,stats_hhmx(4,FileID));
            fprintf(fid,'\r\n');
            
            fprintf('\t%s','(done)');
            fprintf('\n');
    end
%====== Adaptive thresholding - Hist max
    if DoThrByHMX
            fprintf('\t%s','Adaptive threshold > Histogram max');
            Cluster_thr_histmax = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,:)));
            SaveToFolder = fullfile('D - Threshold','Histogram peak');
            SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_histmax),' (Histogram peak - in ROI)'];
            RefDataCol = ProcSettings.invlogAFDCol;

            [histmax_pts_over, histmax_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_histmax,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
            num_events_over = size(histmax_pts_over,1);
            
            % update the summary table
            stats_hmx(1,FileID) = Cluster_thr_histmax;
            stats_hmx(2,FileID) = num_events_over;
            stats_hmx(3,FileID) = total_events_table;
            stats_hmx(4,FileID) = (num_events_over/total_events_table) * 100;
            
            % save results for this cell to a file
            fprintf(fid,DataFormat,'Histogram peak',Cluster_thr_histmax,stats_hmx(4,FileID));
            fprintf(fid,'\r\n');
            
            fprintf('\t\t%s','(done)');
            fprintf('\n');
    end
%====== Adaptive thresholding - Mean minus 1*StDev
    if DoThrByMMSD
            fprintf('\t%s','Adaptive threshold > Mean minus StDev');
            Cluster_thr_meanminus1stdev = mean(data_proc(in_roi,ProcSettings.invlogAFDCol)) - std(data_proc(in_roi,ProcSettings.invlogAFDCol));
            SaveToFolder = fullfile('D - Threshold','Mean minus StDev');
            SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_meanminus1stdev),' (Mean - StDev - in ROI)'];
            RefDataCol = ProcSettings.invlogAFDCol;

            [meanminus1stdev_pts_over, meanminus1stdev_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_meanminus1stdev,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
            num_events_over = size(meanminus1stdev_pts_over,1);
       
            % update the summary table
            stats_mmsd(1,FileID) = Cluster_thr_meanminus1stdev;
            stats_mmsd(2,FileID) = num_events_over;
            stats_mmsd(3,FileID) = total_events_table;
            stats_mmsd(4,FileID) = (num_events_over/total_events_table) * 100;

            % save results for this cell to a file
            fprintf(fid,DataFormat,'Mean minus StDev',Cluster_thr_meanminus1stdev,stats_mmsd(4,FileID));    
            fprintf(fid,'\r\n');
            
            fprintf('\t%s','(done)');
            fprintf('\n');
    end        
%====== Adaptive thresholding - Mean plus 1*StDev
    if DoThrByMPSD
            fprintf('\t%s','Adaptive threshold > Mean plus StDev');

            Cluster_thr_meanplus1stdev = mean(data_proc(in_roi,ProcSettings.invlogAFDCol)) + std(data_proc(in_roi,ProcSettings.invlogAFDCol));
            SaveToFolder = fullfile('D - Threshold','Mean plus StDev');
            SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_meanplus1stdev),' (Mean + StDev - in ROI)'];
            RefDataCol = ProcSettings.invlogAFDCol;

            [meanplus1stdev_pts_over, meanplus1stdev_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_meanplus1stdev,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
            num_events_over = size(meanplus1stdev_pts_over,1);
      
            % update the summary table
            stats_mpsd(1,FileID) = Cluster_thr_meanplus1stdev;
            stats_mpsd(2,FileID) = num_events_over;
            stats_mpsd(3,FileID) = total_events_table;
            stats_mpsd(4,FileID) = (num_events_over/total_events_table) * 100;
            
            % save results for this cell to a file
            fprintf(fid,DataFormat,'Mean plus StDev',Cluster_thr_meanplus1stdev,stats_mpsd(4,FileID));
            fprintf(fid,'\r\n');

            fprintf('\t%s','(done)');
            fprintf('\n');
    end        
%====== Adaptive thresholding - Median
    if DoThrByMEDN
            fprintf('\t%s','Adaptive threshold > Median');

            Cluster_thr_median = median(data_proc(in_roi,ProcSettings.invlogAFDCol));
            SaveToFolder = fullfile('D - Threshold','Median');
            SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_median),' (Median - in ROI)'];
            RefDataCol = ProcSettings.invlogAFDCol;

            [median_pts_over, median_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_median,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
            num_events_over = size(median_pts_over,1);
      
            % update the summary table
            stats_medn(1,FileID) = Cluster_thr_median;
            stats_medn(2,FileID) = num_events_over;
            stats_medn(3,FileID) = total_events_table;
            stats_medn(4,FileID) = (num_events_over/total_events_table) * 100;

            % save results for this cell to a file
            fprintf(fid,DataFormat,'Median',Cluster_thr_median,stats_medn(4,FileID));
            fprintf(fid,'\r\n');
            
            fprintf('\t\t\t\t%s','(done)');
            fprintf('\n');
    end
%====== Fixed thresholding - as per ProcSettings
    if DoThrByFIXED
            fprintf('\t%s',['User-specified threshold > ',num2str(ProcSettings.FixedThreshold)]);

            Cluster_thr_fixed = ProcSettings.FixedThreshold;
            SaveToFolder = fullfile('D - Threshold',['User-specified (',num2str(ProcSettings.FixedThreshold),')']);
            SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_fixed),' (user-specified - in ROI)'];
            RefDataCol = ProcSettings.invlogAFDCol;

            [fixedthr_pts_over, fixedthr_pts_under] = thresholder(data_proc(in_roi,:),Cluster_thr_fixed,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
            num_events_over = size(fixedthr_pts_over,1);
           
            % update the summary table
            stats_thrfix(1,FileID) = Cluster_thr_fixed;
            stats_thrfix(2,FileID) = num_events_over;
            stats_thrfix(3,FileID) = total_events_table;
            stats_thrfix(4,FileID) = (num_events_over/total_events_table) * 100;

            % save results for this cell to a file
            fprintf(fid,DataFormat,'Fixed',Cluster_thr_fixed,stats_thrfix(4,FileID));
            fprintf(fid,'\r\n');
            
            fprintf('\t\t%s','(done)');
            fprintf('\n');
    end
    
    % close the per-cell summary text file
    fclose(fid);
    
end

%====== Save Adaptive thresholding summary data

    % format for headers (first row)
    stringy = '%s';
    for g = 1:(size(FileList,1)-1)
        stringy = strcat(stringy,'\t%s'); % replace \t with a comma for csv
    end

    % half hist max
    if DoThrByHHMX
        fn_summary_proc_save = 'Combined Thresholds (in ROI) - Half Histogram Max.txt';
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
        fprintf(fid,'%s','Combined Thresholds - Half Histogram Max');
        fprintf(fid,'\r\n');
        fprintf(fid,stringy,FileList{:,1});
        fprintf(fid,'\r\n');
        fclose(fid);
        dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_hhmx,'-append','delimiter','\t','precision','%.3f');
    end
    % hist max
    if DoThrByHMX
        fn_summary_proc_save = 'Combined Thresholds (in ROI) - Histogram peak.txt';
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
        fprintf(fid,'%s','Combined Thresholds - Histogram peak');
        fprintf(fid,'\r\n');
        fprintf(fid,stringy,FileList{:,1});
        fprintf(fid,'\r\n');
        fclose(fid);
        dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_hmx,'-append','delimiter','\t','precision','%.3f');
    end
    % mean minus stdev
    if DoThrByMMSD
        fn_summary_proc_save = 'Combined Thresholds (in ROI) - Mean minus StDev.txt';
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
        fprintf(fid,'%s','Combined Thresholds - Mean minus StDev');
        fprintf(fid,'\r\n');
        fprintf(fid,stringy,FileList{:,1});
        fprintf(fid,'\r\n');
        fclose(fid);
        dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_mmsd,'-append','delimiter','\t','precision','%.3f');
    end
    % mean plus stdev
    if DoThrByMPSD
        fn_summary_proc_save = 'Combined Thresholds (in ROI) - Mean plus StDev.txt';
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
        fprintf(fid,'%s','Combined Thresholds - Mean plus StDev');
        fprintf(fid,'\r\n');
        fprintf(fid,stringy,FileList{:,1});
        fprintf(fid,'\r\n');
        fclose(fid);
        dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_mpsd,'-append','delimiter','\t','precision','%.3f');
    end
        
    % median
    if DoThrByMEDN
        fn_summary_proc_save = 'Combined Thresholds (in ROI) - Median.txt';
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
        fprintf(fid,'%s','Combined Thresholds - Median');
        fprintf(fid,'\r\n');
        fprintf(fid,stringy,FileList{:,1});
        fprintf(fid,'\r\n');
        fclose(fid);
        dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_medn,'-append','delimiter','\t','precision','%.3f');    
    end
    % user-specified
    if DoThrByFIXED
        fn_summary_proc_save = ['Combined Thresholds (in ROI) - User Specified (',num2str(ProcSettings.FixedThreshold),').txt'];
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
        fprintf(fid,'%s','Combined Thresholds - User Specified (',num2str(ProcSettings.FixedThreshold),')');
        fprintf(fid,'\r\n');
        fprintf(fid,stringy,FileList{:,1});
        fprintf(fid,'\r\n');
        fclose(fid);
        dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_thrfix,'-append','delimiter','\t','precision','%.3f');    
    end