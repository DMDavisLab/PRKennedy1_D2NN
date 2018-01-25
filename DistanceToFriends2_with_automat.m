%%
% This version includes a folder-list option for batch processing

ThisVersion = 0.5;

%% Begin

% housekeeping
home %clean up the command window
rng('shuffle')  % set the random seed
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary'); % disables a warning about a temporary variable in the parallel processing section. This is fine but warnings scare people so I am turning it off.
warning('off','MATLAB:MKDIR:DirectoryExists'); % don't warn about existing data folders

InfoMessage = ['---------------------The Other Clustering Thingy (v',num2str(ThisVersion),')---------------------'];
disp(InfoMessage);

%% Get list of files

% Find the coords file
    DataDirName = uigetdir(pwd,'Choose your data folder');
    
    if DataDirName ~= 0
        
        cd(DataDirName);

        % open the file containing procsettings
        if exist(fullfile(cd, 'ProcSettings.txt'), 'file') == 0
            errordlg('Cannot find ''ProcSettings.txt'' file for your data. You''ll need one to proceed.');
        else
            ProcSettings = LoadProcSettings;
        end
        
        % get list of txt files
        DirFullList = dir;                              % Get the data for the current directory
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
        
    else
        
        error('Cancelled?! So rude.');
        
    end
       
    clear DirFullList FolderIndex f fname ext GoodIdx InfoMessage


% Display some information
    disp('---------------------------------------------------------');
    InfoMessage=[datestr(fix(clock),'HH:MM:SS'),9,'Processing ' num2str(size(FileList,1)) ' data tables...'];
    disp(InfoMessage);
    disp('---------------------------------------------------------');

    mkdir('Events');
    mkdir('Distance Histograms');
    mkdir('A - Distance to NN');
    mkdir('B - Sum of Distances to NN');
    mkdir('C - Inv Log Sum of Distances to NN');
    mkdir('D - Threshold');
    mkdir('X - Tables');
    mkdir('Y - Threshold Stats');
    mkdir('Z - MAT Files');

% create containers for threshold summaries       row 1 - threshold value for that table's method
    stats_hhmx = zeros(4,size(FileList,1));     % row 2 - total points over threshold
    stats_mmsd = zeros(4,size(FileList,1));     % row 3 - total points in image
    stats_medn = zeros(4,size(FileList,1));     % row 4 - percent over threshold

%% Begin Processing data tables

    PreviousRegionTimestamp = 0; % Initialise the timing tracker.
    tic; %start the clock to track processing time

    for n = 1:size(FileList,1)

        CurrentTableFileName=FileList{n,1};

        % error if the data table txt file doesn't exist
        if exist(CurrentTableFileName, 'file') == 0
            ErrorMessage = ['Problem! Cannot find the data table file ', CurrentTableFileName, '.'];
            errordlg(ErrorMessage);
        end

        datatable = importdata(CurrentTableFileName);

        % Check that we have a useable struct array with the data, if not
        % make one from the data we did import.
        if ~isstruct(datatable) 
            data2 = datatable;
            datatable = struct;
            datatable.data = data2;
            
            datatable.colheaders = cell(1,size(data2,2));
            for i=1:size(data2,2)
                datatable.colheaders{i} = ['Data Col. ' num2str(i)];
            end
            clear i data2
        end
        
        % Determine the image boundaries
        if isfield(ProcSettings,'ImageSize');
            ProcSettings.AxisLimits = [0 ProcSettings.ImageSize 0 ProcSettings.ImageSize];
        else
            % try to guess the limits
            minX = floor(min(datatable.data(:,ProcSettings.xCoordsColumn))/1000)*1000;
            minY = floor(min(datatable.data(:,ProcSettings.yCoordsColumn))/1000)*1000;
            maxX = ceil(max(datatable.data(:,ProcSettings.xCoordsColumn))/1000)*1000;
            maxY = ceil(max(datatable.data(:,ProcSettings.yCoordsColumn))/1000)*1000;
            
            % assuming square images were acquired ...
            minXY = min(minX,minY);
            maxXY = max(maxX,maxY);
            ProcSettings.AxisLimits = [minXY maxXY minXY maxXY];

            clear minX minY maxX maxY minXY maxXY
        end
        ProcSettings.PlotWidth = ProcSettings.AxisLimits(2) - ProcSettings.AxisLimits(1);
        ProcSettings.SaveHighDPI = strcat('-r',num2str(ProcSettings.PlotWidth / 10)); % for 10 nm per pixel

        % set up some other stuff that we'll need...
        ProcSettings.ExptTitle = FileList{n,1};

        % randomising the input rows makes matlab render faster for some
        % unknown reason. The rand idx here are kept for later plots too.
        data_len = length(datatable.data(:,ProcSettings.xCoordsColumn));        % Get the length of the input vectors.
        data_rand_ind = randperm(data_len);                                     % Create a vector of random integers from 1 to len
        data_randx = datatable.data(data_rand_ind,ProcSettings.xCoordsColumn);  % Randomize the input vectors.
        data_randy = datatable.data(data_rand_ind,ProcSettings.yCoordsColumn);
        
        % Adaptive region cropping
        data_density = data_len / (ProcSettings.AxisLimits(2) * ProcSettings.AxisLimits(4));
        ProcSettings.TestRegionSize_for_minFF = ceil(sqrt(ProcSettings.FurthestFriendID / data_density));

%====== Plot the events as they are

        DataEvents_fig = figure;
        if ProcSettings.ChannelIDColumn ~= 0
            data_randz=datatable.data(data_rand_ind,ProcSettings.ChannelIDColumn);
            scatter3(data_randx,data_randy,data_randz,1,'.',data_randz);
        %     scatter3(data(:,ProcSettings.xCoordsColumn),data(:,ProcSettings.yCoordsColumn),data(:,ProcSettings.ChannelIDColumn),1,'.',data(:,ProcSettings.ChannelIDColumn));
        else
            scatter(data_randx,data_randy,1,'.');
        %     scatter(data(:,ProcSettings.xCoordsColumn),data(:,ProcSettings.yCoordsColumn),1,'.');
        end
        box on
        axis(ProcSettings.AxisLimits);
        axis image square
        title([ProcSettings.ExptTitle,' : Events']);
        SaveFileName = [FileList{n,2},' - Events.png'];
        print(DataEvents_fig,'-dpng','-r300',fullfile(pwd,'Events',SaveFileName));
        close(DataEvents_fig)

%====== calculate the distances

        [ppFFD_tmp,ppAFD_tmp,DistancesToN] = DTF2ParaFunc(datatable.data,ProcSettings);

        data_proc = horzcat(datatable.data(:,ProcSettings.xCoordsColumn), datatable.data(:,ProcSettings.yCoordsColumn),ppFFD_tmp,ppAFD_tmp);
        ProcSettings.xCol_procd = 1;
        ProcSettings.yCol_procd = 2;
        ProcTblHeaders = {'x','y','FFD(10)','AFD(10)'};
        ProcSettings.FFDCol = size(data_proc,2)-1; % FFD = furthest friend distance
        ProcSettings.AFDCol = size(data_proc,2); % AFD = all friend distance (sum of distances to n friends)

        clear ppAFD_tmp ppFFD_tmp
        
        total_events_table = size(data_proc,1);

%====== Overlay distance histograms

        disthist_max = ceil(max(max(DistancesToN))/10)*10;
        disthist=[];
        for dn = 1:ProcSettings.FurthestFriendID
             disthist(dn,:) = hist(DistancesToN(dn,:),0:0.1:disthist_max);
        end
        disthist_labels = 0:0.1:disthist_max;

        fade_array = linspace(1,0,ProcSettings.FurthestFriendID)';
        fade_r2b = horzcat(fade_array,zeros(ProcSettings.FurthestFriendID,1),flipud(fade_array));

        disthistoverlay = figure;
        p01 = plot(disthist_labels,disthist(1,:),'-','Color',fade_r2b(1,:));
        hold on
        p02 = plot(disthist_labels,disthist(2,:),':','Color',fade_r2b(2,:));
        p03 = plot(disthist_labels,disthist(3,:),':','Color',fade_r2b(3,:));
        p04 = plot(disthist_labels,disthist(4,:),':','Color',fade_r2b(4,:));
        p05 = plot(disthist_labels,disthist(5,:),':','Color',fade_r2b(5,:));
        p06 = plot(disthist_labels,disthist(6,:),':','Color',fade_r2b(6,:));
        p07 = plot(disthist_labels,disthist(7,:),':','Color',fade_r2b(7,:));
        p08 = plot(disthist_labels,disthist(8,:),':','Color',fade_r2b(8,:));
        p09 = plot(disthist_labels,disthist(9,:),':','Color',fade_r2b(9,:));
        p10 = plot(disthist_labels,disthist(10,:),'-','Color',fade_r2b(10,:));
        title([ProcSettings.ExptTitle,' : Distances to n^{th} Neighbour']);
        xlabel('Distance to n^{th} nearest neighbour (nm)');
        ylabel('Frequency');
        legend([p01 p02 p10],{'n=1','n=2-9','n=10'});
        set(gca,'YLim',[0 ceil(max(max(disthist(:,2:end-1)))/100)*100]);
        set(gca,'XLim',[0 50]);
        SaveFileName = [FileList{n,2},' - Distances to nth Neighbour.png'];
        print(disthistoverlay,'-dpng','-r300',fullfile(pwd,'Distance Histograms',SaveFileName));

        clear p01 p02 p03 p04 p05 p06 p07 p08 p09 p10 disthist disthist_labels disthist_max fade_array fade_r2b
        close(disthistoverlay);

%====== Histograms of the single distance to the Nth NN

        plotcollection = figure;        
        
        %Bin the distances to find the most common distance range.        
        distanceplot = [];
        xbins = 0:1:200; % max(data(:,3)); % bin width = 10 (nm)
        [distanceplot(2,:),distanceplot(1,:)] = hist(data_proc(:,ProcSettings.FFDCol),xbins);
        subplot(2,3,1);
        bar(distanceplot(1,:),distanceplot(2,:));
        axis tight
        title(['d(',num2str(ProcSettings.FurthestFriendID),')']);
        % may need to adjust this by 'in cell' or 'out of cell' discriminator
        normaldistance = distanceplot(1,distanceplot(2,:)==max(distanceplot(2,2:end-1))); % ignoring first and last bins
        axis([0 distanceplot(1,end) 0 1.1*distanceplot(2,distanceplot(1,:) == normaldistance)]);

        % Invert distances to n10th in data col 5
        ProcSettings.invFFDCol = size(data_proc,2) + 1;
        data_proc(:,ProcSettings.invFFDCol) = normaldistance ./ data_proc(:,ProcSettings.FFDCol);
        ProcTblHeaders = horzcat(ProcTblHeaders,{'Norm''d FFD(10)'});

        % Distance to the 10th nearest neighbour
        dindexplot = [];
        xbins = 0:0.02:floor(max(data_proc(:,ProcSettings.invFFDCol)));
        [dindexplot(2,:),dindexplot(1,:)] = hist(data_proc(:,ProcSettings.invFFDCol),xbins);
        subplot(2,3,4);
        bar(dindexplot(1,:),dindexplot(2,:));
        axis tight
        title('peakd10 / d(10)');

%====== Colour map of the single distance to the Nth NN

        % may need to adjust this by 'in cell' or 'out of cell' discriminator
        peakdindex = dindexplot(1,dindexplot(2,:)==max(dindexplot(2,2:end-1))); % ignoring first and last bins
        axis([0 dindexplot(1,end) 0 1.1*dindexplot(2,dindexplot(1,:) == peakdindex)]);
        ColourmapMax = 3*peakdindex;
        
        data_randz = data_proc(data_rand_ind,ProcSettings.invFFDCol); % randomised z data using rand idx calculated in the beginning
        
        fig_invFFD = figure;
        scatter(data_randx,data_randy,1,data_randz,'filled');
        axis(ProcSettings.AxisLimits);
        axis square
        caxis([0 ColourmapMax]);
        colorbar
        title([ProcSettings.ExptTitle,' : Distance to n(10)']);
        
        SaveFileName = [FileList{n,2},' - fig_invFFD.png'];
        print(fig_invFFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'A - Distance to NN',SaveFileName));
        close(fig_invFFD)
        
%====== Histograms of the sum of the distances to the Nth NN
        
        set(0, 'CurrentFigure', plotcollection);% send the focus back to the plots

        % Sum of distances to all NN up to Nth NN
        sumdistanceplot = [];
        xbins = 0:5:5000; % bin width = 10 (nm)
        [sumdistanceplot(2,:),sumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.AFDCol),xbins);
        subplot(2,3,2);
        bar(sumdistanceplot(1,:),sumdistanceplot(2,:));
        axis tight
        title('sumd(10)');
        normalsumdistance = sumdistanceplot(1,sumdistanceplot(2,:)==max(sumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([0 sumdistanceplot(1,end) 0 1.1*sumdistanceplot(2,sumdistanceplot(1,:) == normalsumdistance(1))]);

        % Invert sum of distances to Nth NN
        ProcSettings.invAFDCol = size(data_proc,2) + 1;
        data_proc(:,ProcSettings.invAFDCol) = normalsumdistance(1) ./ data_proc(:,ProcSettings.AFDCol);
        ProcTblHeaders = horzcat(ProcTblHeaders,{'Norm''d AFD(10)'});
        
        sindexplot = [];
        xbins = 0:0.01:2;
        [sindexplot(2,:),sindexplot(1,:)] = hist(data_proc(:,ProcSettings.invAFDCol),xbins);
        subplot(2,3,5);
        bar(sindexplot(1,:),sindexplot(2,:));
        axis tight
        title('peaksumd10 / sumd(10)');
        axis([0 2 0 ceil(1.1*max(sindexplot(2,:)))]);
        
%====== Colour map of the sum of the distances to the Nth NN

        % may need to adjust this by 'in cell' or 'out of cell' discriminator
        peaksindex = sindexplot(1,sindexplot(2,:)==max(sindexplot(2,2:end-1))); % ignoring first and last bins
        ColourmapMax = 5*peaksindex;

        data_randz = data_proc(data_rand_ind,ProcSettings.invAFDCol);%randomised z data
        
        fig_invAFD = figure;
        scatter(data_randx,data_randy,1,data_randz,'filled');
        axis(ProcSettings.AxisLimits);
        axis square
        caxis([0 ColourmapMax]);
        colorbar
        title([ProcSettings.ExptTitle,' : Sum of distances to n(1..10)']);
        
        SaveFileName = [FileList{n,2},' - fig_invAFD.png'];
        print(fig_invAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'B - Sum of Distances to NN',SaveFileName));
        close(fig_invAFD);

%====== Histograms of the log of the sum of the distances to the Nth NN

        set(0, 'CurrentFigure', plotcollection);% send the focus back to the plots

        % Log of sum of distances to Nth NN
        ProcSettings.logAFDCol = size(data_proc,2) + 1;
        data_proc(:,ProcSettings.logAFDCol) = log(data_proc(:,ProcSettings.AFDCol));
        ProcTblHeaders = horzcat(ProcTblHeaders,{'Log AFD(10)'});
        
        logsumdistanceplot = [];
        xbins = floor(min(data_proc(:,ProcSettings.logAFDCol))):0.1:ceil(max(data_proc(:,ProcSettings.logAFDCol))); % bin width = 10 (nm)
        [logsumdistanceplot(2,:),logsumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.logAFDCol),xbins);
        subplot(2,3,3);
        bar(logsumdistanceplot(1,:),logsumdistanceplot(2,:));
        axis tight
        title('logsumd(10)');
        normallogsumdistance = logsumdistanceplot(1,logsumdistanceplot(2,:)==max(logsumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.logAFDCol))) ceil(max(data_proc(:,ProcSettings.logAFDCol))) 0 1.1*logsumdistanceplot(2,logsumdistanceplot(1,:) == normallogsumdistance)]);

        % Inverse log of sum of distances to Nth NN (i.e. low distances == high values)
        ProcSettings.invlogAFDCol = size(data_proc,2) + 1;
        MaxLogAFD = max(data_proc(:,ProcSettings.logAFDCol));
        data_proc(:,ProcSettings.invlogAFDCol) = MaxLogAFD - data_proc(:,ProcSettings.logAFDCol);
        ProcTblHeaders = horzcat(ProcTblHeaders,{'Norm''d Log AFD(10)'});

        invlogsumdistanceplot = [];
        xbins = floor(min(data_proc(:,ProcSettings.invlogAFDCol))):0.1:ceil(max(data_proc(:,ProcSettings.invlogAFDCol))); % bin width = 10 (nm)
        [invlogsumdistanceplot(2,:),invlogsumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.invlogAFDCol),xbins);
        subplot(2,3,6);
        bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:));
        axis tight
        title('invlogsumd(10)');
        normallogsumdistance2 = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.invlogAFDCol))) ceil(max(data_proc(:,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == normallogsumdistance2)]);

%====== Colour map of the inv log sum of the distances to the Nth NN

        data_randz = data_proc(data_rand_ind,ProcSettings.invlogAFDCol);
        
        fig_invlogAFD = figure;
        scatter(data_randx,data_randy,1,data_randz,'filled');
        axis(ProcSettings.AxisLimits);
        axis square
        colorbar
        title([ProcSettings.ExptTitle,' : Log of Sum of distances to n(1..10)']);
        
        SaveFileName = [FileList{n,2},' - fig_invlogAFD (Bright Background).png'];
        print(fig_invlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'C - Inv Log Sum of Distances to NN',SaveFileName));
        
        % Apply a dark background and save again
        caxis([(mean(data_proc(:,ProcSettings.invlogAFDCol))-(2 * std(data_proc(:,ProcSettings.invlogAFDCol)))) (mean(data_proc(:,ProcSettings.invlogAFDCol))+(2 * std(data_proc(:,ProcSettings.invlogAFDCol))))]);
        c = colormap;
        set(gca,'Color',c(1,:));
        set(gcf,'InvertHardcopy','off')

        SaveFileName = [FileList{n,2},' - fig_invlogAFD (Dark Background).png'];
        print(fig_invlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'C - Inv Log Sum of Distances to NN',SaveFileName));

        close(fig_invlogAFD)
        
%====== Save the histogram collection

        set(0, 'CurrentFigure', plotcollection);
        SaveFileName = [FileList{n,2},' - NN histogram collection.png'];
        print(plotcollection,'-dpng','-r600',fullfile(pwd,'Distance Histograms',SaveFileName));
        close(plotcollection)
        
        % new histogram plot for InvLogSumDistNN alone
        fig_plot_invlogAFD = figure;
        bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:));
        axis tight
        title([ProcSettings.ExptTitle,' : Inv. Log of Sum of distances to NN_1 to NN_{',num2str(ProcSettings.FurthestFriendID),'}']);
        normallogsumdistance2 = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.invlogAFDCol))) ceil(max(data_proc(:,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == normallogsumdistance2)]);
        SaveFileName = [FileList{n,2},' - InvLogSumDist-NN histogram.png'];
        print(fig_plot_invlogAFD,'-dpng','-r600',fullfile(pwd,'Distance Histograms',SaveFileName));
        close(fig_plot_invlogAFD)
        
%====== Adaptive thresholding - Half hist max

        % Cluster threshold (auto at half max histogram)
        Cluster_thr_halfhistmax = (invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,:))))/2;
        plot_pts_in_halfhistmax = data_proc(data_proc(:,ProcSettings.invlogAFDCol)>=Cluster_thr_halfhistmax,:);
        percent_in_halfhistmax = (size(plot_pts_in_halfhistmax,1)/total_events_table*100);
        
        % update the summary table
        stats_hhmx(1,n) = Cluster_thr_halfhistmax;
        stats_hhmx(2,n) = size(plot_pts_in_halfhistmax,1);
        stats_hhmx(3,n) = total_events_table;
        stats_hhmx(4,n) = percent_in_halfhistmax;
        
        %Get the length of the input vectors.
            data_halfhistmax_len = length(plot_pts_in_halfhistmax(:,1));
        %Create a vector of random integers from 1 to len
            data_halfhistmax_rand_ind=randperm(data_halfhistmax_len);
        %Randomize the input vectors.
            data_halfhistmax_randx=plot_pts_in_halfhistmax(data_halfhistmax_rand_ind,ProcSettings.xCol_procd);
            data_halfhistmax_randy=plot_pts_in_halfhistmax(data_halfhistmax_rand_ind,ProcSettings.yCol_procd);
            data_halfhistmax_randz=plot_pts_in_halfhistmax(data_halfhistmax_rand_ind,ProcSettings.invlogAFDCol);
          
        fig_plot_pts_in_halfhistmax = figure;
        scatter(data_halfhistmax_randx,data_halfhistmax_randy,1,data_halfhistmax_randz,'filled');
        axis equal
        axis(ProcSettings.AxisLimits);

        % block out the below threshold values in the colourbar
        cbar = colormap;
        thold_posn = floor(Cluster_thr_halfhistmax / max(plot_pts_in_halfhistmax(:,ProcSettings.invlogAFDCol)) * 64);
        cbar(1:thold_posn,:) = repmat(cbar(1,:), thold_posn,1);
        colormap(cbar);
        colorbar
        caxis([0 max(plot_pts_in_halfhistmax(:,ProcSettings.invlogAFDCol))]);
        
        title([ProcSettings.ExptTitle,' : Points above threshold=',num2str(Cluster_thr_halfhistmax),' (Half histogram peak)']);
        set(gca,'Color',[0 0 0.561]);
        set(gcf,'InvertHardcopy','off')
        SaveFileName = [FileList{n,2},' - events over threshold=',num2str(Cluster_thr_halfhistmax),' (Half hist peak).png'];
        print(fig_plot_pts_in_halfhistmax,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'D - Threshold',SaveFileName));
        
        close(fig_plot_pts_in_halfhistmax);
        clear data_halfhistmax_len data_halfhistmax_rand_ind data_halfhistmax_randx data_halfhistmax_randy data_halfhistmax_randz cbar thold_posn SaveFileName



%====== Adaptive thresholding - Mean-1*StDev

        % Cluster threshold (mean minus 1 stdev)
        Cluster_thr_meanminus1stdev = mean(data_proc(:,ProcSettings.invlogAFDCol)) - std(data_proc(:,ProcSettings.invlogAFDCol));
        plot_pts_in_meanminus1stdev = data_proc(data_proc(:,ProcSettings.invlogAFDCol)>=Cluster_thr_meanminus1stdev,:);
        percent_in_meanminusstdev = (size(plot_pts_in_meanminus1stdev,1)/total_events_table*100);

        % update the summary table
        stats_mmsd(1,n) = Cluster_thr_meanminus1stdev;
        stats_mmsd(2,n) = size(plot_pts_in_meanminus1stdev,1);
        stats_mmsd(3,n) = total_events_table;
        stats_mmsd(4,n) = percent_in_meanminusstdev;
        
        %Get the length of the input vectors.
            data_meanminus1stdev_len = length(plot_pts_in_meanminus1stdev(:,1));
        %Create a vector of random integers from 1 to len
            data_meanminus1stdev_rand_ind=randperm(data_meanminus1stdev_len);
        %Randomize the input vectors.
            data_meanminus1stdev_randx=plot_pts_in_meanminus1stdev(data_meanminus1stdev_rand_ind,ProcSettings.xCol_procd);
            data_meanminus1stdev_randy=plot_pts_in_meanminus1stdev(data_meanminus1stdev_rand_ind,ProcSettings.yCol_procd);
            data_meanminus1stdev_randz=plot_pts_in_meanminus1stdev(data_meanminus1stdev_rand_ind,ProcSettings.invlogAFDCol);

        fig_plot_pts_in_meanminus1stdev = figure;
        scatter(data_meanminus1stdev_randx,data_meanminus1stdev_randy,1,data_meanminus1stdev_randz,'filled');
        % scatter(plot_pts_in_meanminus1stdev(:,ProcSettings.xCoordsColumn),plot_pts_in_meanminus1stdev(:,ProcSettings.yCoordsColumn),1,plot_pts_in_meanminus1stdev(:,ProcSettings.invlogAFDCol),'filled');
        axis equal
        axis(ProcSettings.AxisLimits);
        
        % block out the below threshold values in the colourbar
        cbar = colormap;
        thold_posn = floor(Cluster_thr_meanminus1stdev / max(plot_pts_in_meanminus1stdev(:,ProcSettings.invlogAFDCol)) * 64);
        cbar(1:thold_posn,:) = repmat(cbar(1,:), thold_posn,1);
        colormap(cbar);
        colorbar
        caxis([0 max(plot_pts_in_meanminus1stdev(:,ProcSettings.invlogAFDCol))]);
        
        title([ProcSettings.ExptTitle,' : Points above threshold=',num2str(Cluster_thr_meanminus1stdev),' (Mean - StDev)']);
        set(gca,'Color',cbar(1,:));
        set(gcf,'InvertHardcopy','off')
        SaveFileName = [FileList{n,2},' - events over threshold=',num2str(Cluster_thr_meanminus1stdev),' (Mean - StDev).png'];
        print(fig_plot_pts_in_meanminus1stdev,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'D - Threshold',SaveFileName));
        
        close(fig_plot_pts_in_meanminus1stdev);
        clear data_meanminus1stdev_len data_meanminus1stdev_rand_ind data_meanminus1stdev_randx data_meanminus1stdev_randy data_meanminus1stdev_randz cbar thold_posn SaveFileName
    
%====== Adaptive thresholding - Median

        % Cluster threshold (median)
        Cluster_thr_median = median(data_proc(:,ProcSettings.invlogAFDCol)); % - std(data(:,ProcSettings.invlogAFDCol));
        plot_pts_in_median = data_proc(data_proc(:,ProcSettings.invlogAFDCol)>=Cluster_thr_median,:);
        percent_in_median = (size(plot_pts_in_median,1)/total_events_table*100);

        
        % update the summary table
        stats_medn(1,n) = Cluster_thr_median;
        stats_medn(2,n) = size(plot_pts_in_median,1);
        stats_medn(3,n) = total_events_table;
        stats_medn(4,n) = percent_in_median;

        %Get the length of the input vectors.
            data_median_len = length(plot_pts_in_median(:,1));
        %Create a vector of random integers from 1 to len
            data_median_rand_ind=randperm(data_median_len);
        %Randomize the input vectors.
            data_median_rand_indx=plot_pts_in_median(data_median_rand_ind,ProcSettings.xCol_procd);
            data_median_rand_indy=plot_pts_in_median(data_median_rand_ind,ProcSettings.yCol_procd);
            data_median_rand_indz=plot_pts_in_median(data_median_rand_ind,ProcSettings.invlogAFDCol);

        fig_plot_pts_in_median = figure;
        scatter(data_median_rand_indx,data_median_rand_indy,1,data_median_rand_indz,'filled');
        % scatter(plot_pts_in_meanminus1stdev(:,ProcSettings.xCoordsColumn),plot_pts_in_meanminus1stdev(:,ProcSettings.yCoordsColumn),1,plot_pts_in_meanminus1stdev(:,ProcSettings.invlogAFDCol),'filled');
        axis equal
        axis(ProcSettings.AxisLimits);
        
        % block out the below threshold values in the colourbar
        cbar = colormap;
        thold_posn = floor(Cluster_thr_median / max(plot_pts_in_median(:,ProcSettings.invlogAFDCol)) * 64);
        cbar(1:thold_posn,:) = repmat(cbar(1,:), thold_posn,1);
        colormap(cbar);
        colorbar
        caxis([0 max(plot_pts_in_median(:,ProcSettings.invlogAFDCol))]);
        
        title([ProcSettings.ExptTitle,' : Points above threshold=',num2str(Cluster_thr_median),' (Median)']);
        set(gca,'Color',[0 0 0.561]);
        set(gcf,'InvertHardcopy','off')
        SaveFileName = [FileList{n,2},' - events over threshold=',num2str(Cluster_thr_median),' (Median).png'];
        print(fig_plot_pts_in_median,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'D - Threshold',SaveFileName));
        
        close(fig_plot_pts_in_median);
        clear data_median_len data_median_rand_ind data_median_rand_indx data_median_rand_indy data_median_rand_indz cbar thold_posn SaveFileName
        
%====== Save processed data tables

        stringy = '%s';
        for g = 1:(length(ProcTblHeaders)-1)
            stringy = strcat(stringy,'\t%s'); % replace \t with a comma for csv
        end
        
        fn_data_proc_save = [FileList{n,2},'_proc.txt'];
        fid = fopen(fullfile(pwd,'X - Tables',fn_data_proc_save),'w');
        fprintf(fid,stringy,ProcTblHeaders{:});
        fprintf(fid,'\r\n');
        fclose(fid);
        dlmwrite(fullfile(pwd,'X - Tables',fn_data_proc_save),data_proc,'-append','delimiter','\t');
            
        fnsave = [FileList{n,2},'.mat'];
        save(fullfile(pwd,'Z - MAT Files',fnsave),'data_proc','ProcTblHeaders','plot_pts_in_halfhistmax','plot_pts_in_meanminus1stdev','plot_pts_in_median');

%====== Save Adaptive thresholding data

        HeaderFormat = '%s\t%s\t%s';
        fn_stats_save = [FileList{n,2},' - Threshold Data.txt'];
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_stats_save),'w');

        fprintf(fid,HeaderFormat,'[Threshold Method]','[Threshold Value]','[% In Clusters]');
        fprintf(fid,'\r\n');
        
        fprintf(fid,HeaderFormat,'Half hist. peak',num2str(Cluster_thr_halfhistmax),num2str(percent_in_halfhistmax));
        fprintf(fid,'\r\n');

        fprintf(fid,HeaderFormat,'Mean minus StDev',num2str(Cluster_thr_meanminus1stdev),num2str(percent_in_meanminusstdev));
        fprintf(fid,'\r\n');
        
        fprintf(fid,HeaderFormat,'Median',num2str(Cluster_thr_median),num2str(percent_in_median));
        fprintf(fid,'\r\n');
        
        fclose(fid);
    
%====== Show progress tracking information

        ThisRegionTimestamp = toc;

        InfoMessage =  [datestr(fix(clock),'HH:MM:SS'),9,'Completed processing Table ',num2str(n),'.'];
        disp(InfoMessage);
        
        InfoMessage =  [9,9,9,sprintf('%0.2f',((ThisRegionTimestamp - PreviousRegionTimestamp)/60)),' minutes, ',num2str(total_events_table),' events, ',num2str(floor(total_events_table/(ThisRegionTimestamp - PreviousRegionTimestamp))),' events per second.'];
        disp(InfoMessage);
        
        InfoMessage = [9,9,9,num2str(round(n / size(FileList,1) * 100)) '% processed.'];
        disp(InfoMessage);
        
        disp('---------------------------------------------------------');
        
        % Update the tracking info to reflect the newly finished region
        PreviousRegionTimestamp = ThisRegionTimestamp;

        % clean up your mess
        clear Cluster_thr_halfhistmax Cluster_thr_meanminus1stdev Cluster_thr_median ColourmapMax c data_density data_len data_proc data_randx data_randy data_randz data_rand_ind datatable dindexplot distanceplot DistancesToN disthistoverlay dn fid fn_data_proc_save fn_stats_save fnsave g HeaderFormat InfoMessage invlogsumdistanceplot logsumdistanceplot MaxLogAFD normaldistance normallogsumdistance normallogsumdistance2 normalsumdistance peakdindex peaksindex percent_in_halfhistmax percent_in_meanminusstdev percent_in_median plot_pts_in_halfhistmax plot_pts_in_meanminus1stdev plot_pts_in_median ProcTblHeaders sindexplot stringy sumdistanceplot total_events_region xbins
        clear DataEvents_fig fig_invAFD fig_invFFD fig_invlogAFD fig_plot_invlogAFD fig_plot_pts_in_halfhistmax fig_plot_pts_in_meanminus1stdev fig_plot_pts_in_median plotcollection

    end % end of this data table processing
    
%% Finish up

%====== Save Adaptive thresholding summary data

    % format for headers (first row)
    stringy = '%s';
    for g = 1:(size(FileList,1)-1)
        stringy = strcat(stringy,'\t%s'); % replace \t with a comma for csv
    end

    % half hist max
    fn_summary_proc_save = 'Combined Thresholds - Half Histogram Max.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Half Histogram Max');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_hhmx,'-append','delimiter','\t');

    % half hist max
    fn_summary_proc_save = 'Combined Thresholds - Mean minus StDev.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Mean minus StDev');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_mmsd,'-append','delimiter','\t');

    % half hist max
    fn_summary_proc_save = 'Combined Thresholds - Median.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Median');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_medn,'-append','delimiter','\t');

%====== Summarise processing time
    
    ExecTime=toc;
    InfoMessage=[datestr(fix(clock),'HH:MM:SS'),9,'Finished processing ' num2str(size(FileList,1)) ' data tables.'];
    disp(InfoMessage);

    InfoMessage=[9,9,9,'Total time: ' sprintf('%0.2f',(ExecTime/60)) ' minutes.'];
    disp(InfoMessage);

    InfoMessage=[9,9,9,'Average time: ' sprintf('%0.2f',(ExecTime/60)/size(FileList,1)) ' minutes per region.'];
    disp(InfoMessage);

    disp('---------------------------------------------------------');

    clear InfoMessage



%% Graveyard


% 
%         DiffsToN = zeros(ProcSettings.FurthestFriendID,size(data,1));
% 
%         DiffsToN(1,:) = DistancesToN(1,:);
% 
%         for r = 2:ProcSettings.FurthestFriendID
%             DiffsToN(r,:) = DistancesToN(r,:) - DistancesToN(r-1,:);
%         end
% 
%         % sum(max(DiffsToN))/numel(DiffsToN);
% 
%         diffshist=[];
%         for n = 1:ProcSettings.FurthestFriendID
%             diffshist(n,:) = hist(DiffsToN(n,:),[0:0.1:10]);
%         end
%         disthist_labels = 0:0.1:10;
% 
%         diffshistoverlay = figure;
%         p1 = plot(disthist_labels,diffshist(1,:),'r');
%         hold on
%         p2 = plot(disthist_labels,diffshist(2,:),'Color',[0.5 0.5 0.5]);
%         plot(disthist_labels,diffshist(3,:),'Color',[0.5 0.5 0.5]);
%         plot(disthist_labels,diffshist(4,:),'Color',[0.5 0.5 0.5]);
%         plot(disthist_labels,diffshist(5,:),'Color',[0.5 0.5 0.5]);
%         plot(disthist_labels,diffshist(6,:),'Color',[0.5 0.5 0.5]);
%         plot(disthist_labels,diffshist(7,:),'Color',[0.5 0.5 0.5]);
%         plot(disthist_labels,diffshist(8,:),'Color',[0.5 0.5 0.5]);
%         plot(disthist_labels,diffshist(9,:),'Color',[0.5 0.5 0.5]);
%         p10 = plot(disthist_labels,diffshist(10,:),'b');
%         title([ProcSettings.ExptTitle,' : Distance to the next neighbour']);
%         legend([p1 p2 p10],{'n=1','n=2-9','n=10'});
%         set(gca,'YLim',[0 ceil(max(max(diffshist(:,2:end-1)))/100)*100]);
%         SaveFileName = [ProcSettings.ExptTitle,' - Distance to next neighbour.png'];
%         print(gcf,'-dpng','-r150',SaveFileName);
% 
        % clear p1 p2 p10 diffshist diffshist_labels
        % close(diffshistoverlay);

% for g = 1:size(data,1)
%     if data(g,3) <= ClusThresh
%         data(g,9) = 1;
%     else
%         data(g,9) = 0;
%     end
% end





% ColourmapMax = 5*normalsumdistance(1);
% figure;
% title('sumd(10)');
% scatter(data(:,ProcSettings.xCoordsColumn),data(:,ProcSettings.yCoordsColumn),10,data(:,ProcSettings.AFDCol),'filled');
% axis square tight
% colorbar
% caxis([0 ColourmapMax]);
% 
%         % a quick randomising routine
%         totalpts = size(data,1);
%         totalarea = (ProcSettings.AxisLimits(2) - ProcSettings.AxisLimits(1))*(ProcSettings.AxisLimits(4) - ProcSettings.AxisLimits(3));
%         density = totalpts / totalarea;
%         testmax = 500;
%         testpoints = floor((testmax*testmax) * density);
%         data_rnd = randi([0 testmax*10],testpoints,2)/10;
%         rndvars = vars;
%         [ppFFD_tmp,ppAFD_tmp,ppDistancesToN] = DTF2ParaFunc(data_rnd,rndvars);
%         data_rnd = horzcat(data_rnd,ppFFD_tmp,ppAFD_tmp);   % data now has cols: x,y,FFD(10),AFD(10)
%         rndDistancesToN = ppDistancesToN;
%         clear ppAFD_tmp ppFFD_tmp ppDistancesToN rndvars testmax testpoints density totalarea totalpts
%         rnd_mediansumdistance = median(data_rnd(:,3));









% 
% 
% % mangle the colorbar into thresholds
% fig_invlogAFD_threshold_colours = figure;
% scatter(data_randx,data_randy,1,data_randz,'filled');
% % scatter(data(:,ProcSettings.xCoordsColumn),data(:,ProcSettings.yCoordsColumn),1,data(:,ProcSettings.invlogAFDCol),'filled');
% axis(ProcSettings.AxisLimits);
% axis square
% colorbar
% 
% % cbar = zeros(64,3);
% % cbarmin = 0;
% % cbarthr0 = 0.3;
% % cbarthr1 = 1.25;
% % cbarthr2 = 2.125;
% % cbarmax = 3.5;
% 
% cbar = zeros(64,3);
% cbarmin = 0;
% cbarthr1 = 1.0;
% cbarthr2 = 2.0;
% cbarthr3 = 3.0;
% cbarmax = 4.0;
% caxis([cbarmin cbarmax]);
% 
% cbar_unit = cbarmax/64;
% cbars1 = 1+floor((cbarthr1 - cbarmin) / cbar_unit);
% cbars2 = floor((cbarthr2 - cbarthr1) / cbar_unit);
% cbars3 = floor((cbarthr3 - cbarthr2) / cbar_unit);
% % cbars4 = floor((cbarmax - cbarthr3) / cbar_unit);
% 
% cbars_thr1 = cbars1;
% cbars_thr2 = cbars_thr1 + cbars2;
% cbars_thr3 = cbars_thr2 + cbars3;
% cbars4 = 64 - cbars_thr3;
% 
% red = [1 0 0];
% green = [0 1 0];
% blue = [0 0 1];
% yellow = [1 1 0];
% magenta = [1 0 1];
% cyan = [0 1 1];
% white = [1 1 1];
% darkblue = [0 0 0.5];
% 
% cbar(1:cbars_thr1,:) = repmat(darkblue,[cbars1 1]);
% cbar(cbars_thr1+1:cbars_thr2,:) = repmat(magenta,[cbars2 1]);
% cbar(cbars_thr2+1:cbars_thr3,:) = repmat(yellow,[cbars3 1]);
% cbar(cbars_thr3+1:end,:) = repmat(cyan,[cbars4 1]);
% 
% colormap(cbar)
% set(gca,'Color','k');
% 
% title([ProcSettings.ExptTitle,' : Log of Sum of distances to n(1..10)']);
% SaveFileName = [ProcSettings.ExptTitle,' - fig_invlogAFD_threshold_colours.png'];
% set(gcf,'InvertHardcopy','off')
% print(fig_invlogAFD_threshold_colours,'-dpng','-r600',SaveFileName);



% %% Save any figure as highres
% %Get the width of the plot to mangle for 1 px/nm
% 
% % send the focus back to the plots
% set(0, 'CurrentFigure', plotcollection);
% SaveHighDPI = strcat('-r',num2str(ProcSettings.PlotWidth / 10));
% % don't force a white background
% % set(gcf,'InvertHardcopy','off')
% % set(gca,'Visible','on');
% SaveFileName = [ProcSettings.ExptTitle,' - plotcollection.png'];
% print('-dpng','-r150',SaveFileName);
% 
% % plotcollection fig_invFFD fig_invAFD fig_invlogAFD fig_invlogAFD_threshold_colours fig_invlogAFD_dark_BG



%% Stats

% % Cluster threshold (manual)
% Cluster_thr = 2.5;


% % Set a bit for events over threshold
%     ProcSettings.clusteredCol = size(data,2) + 1;
% 
%     for d = 1:size(data,1)
%         if data(d,ProcSettings.invlogAFDCol) >= Cluster_thr_meanminus1stdev
%             data(d,ProcSettings.clusteredCol) = 1;
%         else
%             data(d,ProcSettings.clusteredCol) = 0;
%         end
%     end






% %% fwhm
%     f = x;
%     x=1:size(invlogsumdistanceplot,2);
% f1 = invlogsumdistanceplot(2,:)-0.5*max(invlogsumdistanceplot(2,:));
% % bar(invlogsumdistanceplot(1,:),f1(1,:));
% ind = find(f1(1:end-1).*f1(2:end) <=0);
% FWHM = invlogsumdistanceplot(1,ind(2))-invlogsumdistanceplot(1,ind(1));
% 
% %% Threshold by peak
% % peak is normallogsumdistance2
% 
% Cluster_thr_meanminus1stdev = mean(data_proc(:,ProcSettings.invlogAFDCol)) - std(data_proc(:,ProcSettings.invlogAFDCol));
% 
% ProcSettings.clusteredCol = size(data_proc,2) + 1;
% 
% for d = 1:size(data_proc,1)
%     if data_proc(d,ProcSettings.invlogAFDCol) >= Cluster_thr_meanminus1stdev
%         data_proc(d,ProcSettings.clusteredCol) = 1;
%     else
%         data_proc(d,ProcSettings.clusteredCol) = 0;
%     end
% end
% 
% figure;
% scatter(data_proc(:,ProcSettings.xCoordsColumn),data_proc(:,ProcSettings.yCoordsColumn),3,data_proc(:,ProcSettings.clusteredCol),'filled');
% axis square tight
% 
% 
% %% Map interpolation
% 
% % Create the mesh
% tx=(ProcSettings.AxisLimits(1):ProcSettings.GDInterpSpacing:ProcSettings.AxisLimits(2));
% ty=(ProcSettings.AxisLimits(3):ProcSettings.GDInterpSpacing:ProcSettings.AxisLimits(4));
% [XI, YI] = meshgrid(tx, ty);
% 
% % Interpolate the cluster map surface
% if ProcSettings.UseGriddata == true
%     ZI = griddata(data_proc(:,ProcSettings.xCoordsColumn), data_proc(:,ProcSettings.yCoordsColumn), data_proc(:,ProcSettings.invlogAFDCol), XI, YI, 'v4');
% else
%     % Alternate Inerpolation (much faster than griddata)
%     F = scatteredInterpolant(data_proc(:,ProcSettings.xCoordsColumn), data_proc(:,ProcSettings.yCoordsColumn), data_proc(:,ProcSettings.invlogAFDCol),'natural');
%     ZI = F(XI,YI);
% end
% 
% % Generate the plots
% figure('Color',[1 1 1], 'visible', 'off', 'Renderer', 'OpenGL', 'Units', 'pixels');
% axes('Parent',figure,'Layer','top', 'YTick',zeros(1,0),'XTick',zeros(1,0),'DataAspectRatio', [1,1,1],'position',[0,0,1,1]);            
% box('off');
% hold('on');
% set(gcf, 'PaperUnits', 'inches', 'PaperSize', [10 10], 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 10 10]);
% axis(ProcSettings.AxisLimits);
% axis square image
% set(gca, 'Visible','off');
% 
% %plot the cluster contour map
% [ContourArray,ContourMap] = contourf(XI,YI,ZI,100,'LineColor','none', 'Fill','on');
% axis(ProcSettings.AxisLimits);
% caxis([(mean(data_proc(:,ProcSettings.invlogAFDCol))-(2 * std(data_proc(:,ProcSettings.invlogAFDCol)))) (mean(data_proc(:,ProcSettings.invlogAFDCol))+(2 * std(data_proc(:,ProcSettings.invlogAFDCol))))]);
% % 
% % caxis([0 ColourmapMax]);
% 
% % save data arrays
% save(fullfile(cd, 'KS-uncropped-PLL-allevents-invlogsumd10-dataarray.mat'),'data');
% 
% %Save the 2D Contour Array for later processing
% save(fullfile(cd, 'KS-uncropped-Rituximab-allevents-invlogsumd10-contourarray.mat'),'ContourArray');
% 
% %plot points onto the colourmap
% hold on
% PlotPoints = plot(data_proc(:,ProcSettings.xCoordsColumn), data_proc(:,ProcSettings.yCoordsColumn), 'Marker','.','MarkerSize',4,'LineStyle','none','Color',[0 0 0]);
% set(PlotPoints, 'MarkerSize',1);
% 
% %Get the width of the plot to mangle for 1 px/nm
% 
% SaveHighDPI = strcat('-r',num2str(ProcSettings.PlotWidth / 10));
% 
% %write to a temp image
% print('-dpng',SaveHighDPI,'KS Rituxi by cluster.png');
% 
% 
% %% Thresholding
% ColourmapAxes = caxis;
% ThresholdHigh = mean(data_proc(:,ProcSettings.invlogAFDCol)) - (0.5*std(data_proc(:,ProcSettings.invlogAFDCol)));
% ColormapRes = 256;
% ThresholdMap_High = [1,0,0];
% ThresholdMap_Low = [1,1,1];
% ThresholdImageBG = [0,0,0];
% 
% 
% 
% % Change Threshold for Clusters
%     ThresholdHighMap = zeros(ColormapRes,3);
%     CMapThresholdIdx = round(ColormapRes * (ThresholdHigh / ColourmapAxes(2)));
% 
%     for c_ind = 1:CMapThresholdIdx
%         ThresholdHighMap(c_ind,:) = ThresholdMap_Low;
%     end
% 
%     for c_ind2 = CMapThresholdIdx+1:ColormapRes
%         ThresholdHighMap(c_ind2,:) = ThresholdMap_High;
%     end
% 
%     colormap(ThresholdHighMap);
%     set(gca,'color',[0.25 0.25 0.25]);
%     set(gcf,'color',ThresholdImageBG);
%     set(gca, 'Visible','on');
%     set(gcf, 'InvertHardCopy', 'off'); % This stops MATLAB setting a white background
% 
%     %write to a temp image
%     print('-dpng',SaveHighDPI,'testdata threshold map with points.png');
% 

