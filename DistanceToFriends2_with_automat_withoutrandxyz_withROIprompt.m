%%
% This version includes a folder-list option for batch processing

ThisVersion = 0.5;
DoubleCheckAllSettings = true;
FigureVisibility = 'on';


%% Begin

% housekeeping
home %clean up the command window
rng('shuffle')  % set the random seed
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary'); % disables a warning about a temporary variable in the parallel processing section. This is fine but warnings scare people so I am turning it off.
warning('off','MATLAB:MKDIR:DirectoryExists'); % don't warn about existing data folders

InfoMessage = ['---------------------The Other Clustering Thingy (v',num2str(ThisVersion),')---------------------'];
disp(InfoMessage);

%% Get list of files

% Choose the Phase 1 folder
    DataDirName = uigetdir(pwd,'Choose your data folder');
    
    if DataDirName ~= 0
        
        cd(DataDirName);
        
        % open the file containing procsettings
        if exist(fullfile(cd, 'ProcSettings.txt'), 'file') == 0
            errordlg('Cannot find ''ProcSettings.txt'' file for your data. You''ll need one to proceed.');
            error('Cannot find ''ProcSettings.txt'' file for your data. You''ll need one to proceed.');
        else
            ProcSettings = LoadProcSettings;
        end
        
        % open the file containing the Phase 1 files
        if exist(fullfile(DataDirName,'SaveDataFolderName.mat'), 'file') ~= 0
            load(fullfile(DataDirName,'SaveDataFolderName.mat'));
            cd(OutputFolderName);
            load(fullfile(DataDirName,OutputFolderName,'Phase1FileList.mat'));
            InfoMessage = ['Found the output folder: \t', OutputFolderName,'\n'];
            fprintf(InfoMessage);
            clear InfoMessage
            
            % make the rest of the output folders
                mkdir(fullfile(pwd,'Events'));
                mkdir(fullfile(pwd,'Histograms'));
                    mkdir(fullfile(pwd,'Histograms','Distances to Nth NN'));
                    mkdir(fullfile(pwd,'Histograms','Collections'));
                    mkdir(fullfile(pwd,'Histograms','ILSD-NN Histogram'));
                if ProcSettings.DoStrictdNN
                    mkdir(fullfile(pwd,'A - Distance to NN'));
                end
                if ProcSettings.DoSumdNN
                    mkdir(fullfile(pwd,'B - Sum of Distances to NN'));
                end
                if ProcSettings.DoInvLogSumdNN
                    mkdir(fullfile(pwd,'C - Inv Log Sum of Distances to NN'));
                end
                mkdir(fullfile(pwd,'D - Threshold'));
                mkdir(fullfile(pwd,'X - Tables'));
                mkdir(fullfile(pwd,'Y - Threshold Stats'));
        else
            errordlg('Phase 1 folder does not exist! You need to perform Phase 1 operations.');
        end
     
    else
        error('Cancelled?! So rude.');
    end
   

% Display some information
    disp('---------------------------------------------------------');
    InfoMessage=[datestr(fix(clock),'HH:MM:SS'),9,'Processing ' num2str(size(FileList,1)) ' data tables...'];
    disp(InfoMessage);
    disp('---------------------------------------------------------');


%% Phase 2 - Load & check data, load pre-saved ROIs

    PreviousRegionTimestamp = 0; % Initialise the timing tracker.
    mainproc_stopwatch = tic; %start the clock to track processing time

    for FileID = 1:size(FileList,1)

        CurrentTableFileName=FileList{FileID,1};
        InfoMessage =  [datestr(fix(clock),'HH:MM:SS'),9,'[',num2str(FileID),'/',num2str(size(FileList,1)),'] Opening file: ',CurrentTableFileName];
        disp(InfoMessage);
        clear InfoMessage
        
        Phase1Source_filename = [FileList{FileID,2},' - Phase1Files.mat'];
        load(fullfile(DataDirName,OutputFolderName,'Z - MAT Files',Phase1Source_filename));
        clear Phase1Source_filename
        
    %====== Plot the events as they are

        fprintf('%s','Saving image of all events');
            
        [DataEvents_fig, DataEvents_axes] = DoMeAFigure(ProcSettings.AxisLimits,[0 0 0]);
              
        if ProcSettings.ChannelIDColumn ~= 0
            data_randz=datatable.data(data_rand_ind,ProcSettings.ChannelIDColumn);
            EventsPlot = scatter3(data_randx,data_randy,datatable.data(:,ProcSettings.ChannelIDColumn),1,'w.',data_randz);
            set(DataEvents_fig,'InvertHardcopy','off');
            clear data_randz
        else
            EventsPlot = scatter(data_randx,data_randy,1,'w.');
        end

        set(DataEvents_fig,'InvertHardcopy','off');
        SaveFileName = [FileList{FileID,2},' - Events'];
        print(DataEvents_fig,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'Events',[SaveFileName,'.png']));
        clear SaveFileName
        
        fprintf('\t%s','(done)');
        fprintf('\n');

        set(DataEvents_fig,'visible','on');
        set(EventsPlot,'CData',[1 0.02 0.2]);

        %plot the ROI again (to check)
        hold on
        plot(pos(:,1),pos(:,2),'c','LineWidth',2)
        SaveFileName = [FileList{FileID,2},' - ROI Preview'];
        print(DataEvents_fig,'-dpng',ProcSettings.SaveLowDPI,fullfile(pwd,'Events',[SaveFileName,'.png']));
        close(DataEvents_fig);
        
        % get events within this region [InROI_idx, OnROI_idx]
        [InROI_idx,~] = inpolygon(datatable.data(:,ProcSettings.xCoordsColumn),datatable.data(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));

        ROIarea = polyarea(pos(:,1),pos(:,2));
    %     EventDensityInROI = sum(InROI_idx) / ROIarea; % events per area

        clear DataEvents_fig DataEvents_axes SaveFileName tmpax OnROI_idx xl yl xc yc xsc ysc xdir ydir hFH fnsave fcn EventsPlot

        fprintf('\t%s','(done)');
        fprintf('\n');

        
    end
    
   
%% ====== Calculate the distances

% %         % for all events in the entire image table
%         [ppFFD_tmp,ppAFD_tmp,DistancesToN] = DTF2ParaFunc(datatable.data,ProcSettings);
%         data_proc = horzcat(datatable.data(:,ProcSettings.xCoordsColumn), datatable.data(:,ProcSettings.yCoordsColumn),ppFFD_tmp,ppAFD_tmp);

       % It is about twice as fast to only measure for events that are in the ROI. Non-ROI events get a NaN distance result
        [ppFFD_inROI_tmp,ppAFD_inROI_tmp,DistancesToN_inROI] = DTF2ParaFunc(datatable.data(InROI_idx,:),ProcSettings,datatable.data);
        data_proc = horzcat(datatable.data(:,ProcSettings.xCoordsColumn), datatable.data(:,ProcSettings.yCoordsColumn));
        data_proc(:,end+1:end+2) = nan;
        data_proc(InROI_idx,end-1:end) = horzcat(ppFFD_inROI_tmp,ppAFD_inROI_tmp);
        
        ProcSettings.xCol_procd = 1;
        ProcSettings.yCol_procd = 2;
        ProcSettings.FFDCol = size(data_proc,2)-1; % FFD = furthest friend distance
        ProcSettings.AFDCol = size(data_proc,2); % AFD = all friend distance (sum of distances to n friends)

        ProcTblHeaders = {'x','y',['FFD(',num2str(ProcSettings.FurthestFriendID),')'],['AFD(',num2str(ProcSettings.FurthestFriendID),')']};
        horzcat(ProcTblHeaders,{['Norm''d FFD(',num2str(ProcSettings.FurthestFriendID),')']});
        
        clear ppAFD_tmp ppFFD_tmp ppFFD_inROI_tmp ppAFD_inROI_tmp data_proc_inROI
        
    %====== Make a randomised version of this image
         fprintf('%s','Randomising position of events within ROI');
        % generate x and y values to the nearest 10th nm to match the
        % in-ROI density but spread throughout the entire image area
        % We will then adjust the total number of random events until we
        % have the same number in the ROI as the real data.
        rndm_events_reqd = ceil(((ProcSettings.PlotWidth*ProcSettings.PlotWidth)/ROIarea) * sum(InROI_idx));
        datatable.data_rndm = zeros(rndm_events_reqd,size(datatable.data,2));        
        datatable.data_rndm(:,ProcSettings.xCoordsColumn) = randi([ProcSettings.AxisLimits(1)*10, ProcSettings.AxisLimits(2)*10],rndm_events_reqd,1);
        datatable.data_rndm(:,ProcSettings.yCoordsColumn) = randi([ProcSettings.AxisLimits(3)*10, ProcSettings.AxisLimits(4)*10],rndm_events_reqd,1);
        datatable.data_rndm = datatable.data_rndm / 10; % convert vales to nm.

        % find randmised events inside of ROI ... should equal sum(InROI_idx)
        [RndInROI_idx,~] = inpolygon(datatable.data_rndm(:,ProcSettings.xCoordsColumn),datatable.data_rndm(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));
        EventsOffset = sum(RndInROI_idx) - sum(InROI_idx);
        
        % We want the same number of random events within the ROI as the for the data events within the ROI
        % i.e. identical event density. This happens when EventOffset is zero.
        while EventsOffset ~= 0
            if  EventsOffset  < 0
                
                %insufficient events, add more
                EventsToAdd = abs(EventsOffset) * ceil(((ProcSettings.PlotWidth*ProcSettings.PlotWidth)/ROIarea)); % Add events proportional to the ROI area within the image area
                tmp_new_events = zeros(EventsToAdd,size(datatable.data,2));
                tmp_new_events(:,ProcSettings.xCoordsColumn) = randi([ProcSettings.AxisLimits(1)*10, ProcSettings.AxisLimits(2)*10],EventsToAdd,1);
                tmp_new_events(:,ProcSettings.yCoordsColumn) = randi([ProcSettings.AxisLimits(3)*10, ProcSettings.AxisLimits(4)*10],EventsToAdd,1);
                tmp_new_events = tmp_new_events / 10; % convert vales to nm.
                datatable.data_rndm = vertcat(datatable.data_rndm,tmp_new_events);

                [RndInROI_idx,~] = inpolygon(datatable.data_rndm(:,ProcSettings.xCoordsColumn),datatable.data_rndm(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));
                clear EventsToAdd tmp_new_events
                EventsOffset = sum(RndInROI_idx) - sum(InROI_idx);
                
            elseif EventsOffset > 0
                
                % excess events, delete the difference from the list of RndInROI events
                k_idx = find(RndInROI_idx); % indices of datatable.data_rndm that are In ROI.
                EventsToGo = sort(randi([1 size(k_idx,1)],abs(EventsOffset),1)); % pick random indices to remove
%                 figure
%                 scatter(datatable.data_rndm(k_idx,ProcSettings.xCoordsColumn),datatable.data_rndm(k_idx,ProcSettings.yCoordsColumn),5,'r','filled');
%                 hold on
%                 scatter(datatable.data_rndm(k_idx(EventsToGo),ProcSettings.xCoordsColumn),datatable.data_rndm(k_idx(EventsToGo),ProcSettings.yCoordsColumn),5,'y','filled');
%                 axis(ProcSettings.AxisLimits);
                datatable.data_rndm(k_idx(EventsToGo),:) = [];
                [RndInROI_idx,~] = inpolygon(datatable.data_rndm(:,ProcSettings.xCoordsColumn),datatable.data_rndm(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));
                clear k_idx EventsToGo
                EventsOffset = sum(RndInROI_idx) - sum(InROI_idx);
            end
        end
        clear EventsOffset rndm_events_reqd ROIarea pos
              
        [DataEvents_rndroi_fig, DataEvents_rndroi_axes] = DoMeAFigure(ProcSettings.AxisLimits,[0 0 0]);
        scatter(datatable.data_rndm(RndInROI_idx,ProcSettings.xCoordsColumn),datatable.data_rndm(RndInROI_idx,ProcSettings.yCoordsColumn),1,[1 1 1],'.');
        set(DataEvents_rndroi_fig,'InvertHardCopy','off');
        SaveFileName = [FileList{FileID,2},' - Events (In ROI Randomised)'];
        print(DataEvents_rndroi_fig,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'Events',[SaveFileName,'.png']));
        close(DataEvents_rndroi_fig)
        clear SaveFileName DataEvents_rndroi_fig DataEvents_rndroi_axes
        
        fprintf('\t%s','(done)');
        fprintf('\n');
        
%         % measure all the distances for the randomised data (using the entire image field)
%         [ppFFD_rnd_tmp,ppAFD_rnd_tmp,DistancesToN_rndm] = DTF2ParaFunc(datatable.data_rndm,ProcSettings);
%         data_proc_rndm = horzcat(datatable.data_rndm(:,ProcSettings.xCoordsColumn), datatable.data_rndm(:,ProcSettings.yCoordsColumn),ppFFD_rnd_tmp,ppAFD_rnd_tmp);
        
        % measure all the distances for the randomised data
        [ppFFD_rndm_inROI_tmp,ppAFD_rndm_inROI_tmp,DistancesToN_rndm] = DTF2ParaFunc(datatable.data_rndm(RndInROI_idx,:),ProcSettings,datatable.data_rndm);
        data_proc_rndm = horzcat(datatable.data_rndm(:,ProcSettings.xCoordsColumn), datatable.data_rndm(:,ProcSettings.yCoordsColumn));
        data_proc_rndm(:,end+1:end+2) = nan;
        data_proc_rndm(RndInROI_idx,end-1:end) = horzcat(ppFFD_rndm_inROI_tmp,ppAFD_rndm_inROI_tmp);
       
        clear ppAFD_tmp ppFFD_tmp ppFFD_rndm_inROI_tmp ppAFD_rndm_inROI_tmp ppFFD_rnd_tmp ppAFD_rnd_tmp
        clear pos % ROI position information is no longer required (and is save to MAT file)

        % To delete external-ROI points...
%         datatable.data_rndm(~RndInROI_idx,:) = []; 

%   % Simple randomisation of entire image - this results in an unfair redistribution of events within the cell to outside the cell.
%         datatable.data_rndm = zeros(size(datatable.data));
%         datatable.data_rndm(:,ProcSettings.xCoordsColumn) = randi([ProcSettings.AxisLimits(1)*10, ProcSettings.AxisLimits(2)*10],total_events_table,1);
%         datatable.data_rndm(:,ProcSettings.yCoordsColumn) = randi([ProcSettings.AxisLimits(3)*10, ProcSettings.AxisLimits(4)*10],total_events_table,1);
%         datatable.data_rndm = datatable.data_rndm / 10; % convert vales to nm.
        


        
%====== Overlay distance histograms

        fprintf('%s ', 'Saving images:');
        fprintf('\n');

%     
%     % All Events in Image
%         fprintf('\t%s','NN Distance histograms (all events)'); 
%
%         disthist_max = ceil(max(max(DistancesToN))/10)*10;
%         disthist=[];
%         for dn = 1:ProcSettings.FurthestFriendID
%              disthist(dn,:) = hist(DistancesToN(dn,:),0:0.1:disthist_max);
%         end
%         disthist = disthist ./ total_events_table * 100;
%         disthist_labels = 0:0.1:disthist_max;
% 
%         fade_array = linspace(1,0,ProcSettings.FurthestFriendID)';
%         fade_r2b = horzcat(fade_array,zeros(ProcSettings.FurthestFriendID,1),flipud(fade_array));
% 
%         disthistoverlay = figure('Visible',FigureVisibility);
%         p01 = plot(disthist_labels,disthist(1,:),'-','Color',fade_r2b(1,:));
%         hold on
%         p02 = plot(disthist_labels,disthist(2,:),':','Color',fade_r2b(2,:));
%         for p=3:ProcSettings.FurthestFriendID-1
%             plot(disthist_labels,disthist(p,:),':','Color',fade_r2b(p,:));
%         end
%         p_end = plot(disthist_labels,disthist(ProcSettings.FurthestFriendID,:),'-','Color',fade_r2b(10,:));
%         title([ProcSettings.ExptTitle,' : Distances to n^{th} Neighbour']);
%         xlabel('Distance to n^{th} nearest neighbour (nm)');
%         ylabel('Frequency');
%         legend([p01 p02 p_end],{'n=1',['n=2-',num2str(ProcSettings.FurthestFriendID-1)],['n=',num2str(ProcSettings.FurthestFriendID),'']});
% %         set(gca,'YLim',[0 ceil(max(max(disthist(:,2:end-1)))/100)*100]);
%         set(gca,'YLim',[0 4]);
%         set(gca,'XLim',[0 50]); % need to calc the upper limit automagically
%         set(gca,'TickDir','out');
%         xlabel('Distance to n^{th} neighbour (nm)');
%         ylabel('Percent of Events');
%               
%         SaveFileName = [FileList{FileID,2},' - Distances to nth Neighbour'];
%         print(disthistoverlay,'-dpng','-r300',fullfile(pwd,'Histograms','Distances to Nth NN',SaveFileName));
% 
%         fprintf('\t%s','(done)');
%         fprintf('\n');

    % Only Events in ROI
    
        fprintf('\t%s','NN Distance histograms (events within ROI)'); 

        disthist_max_inROI = ceil(max(max(DistancesToN_inROI))/10)*10;
        disthist_in_roi=[];
        for dn = 1:ProcSettings.FurthestFriendID
             disthist_in_roi(dn,:) = hist(DistancesToN_inROI(dn,:),0:0.1:disthist_max_inROI);
        end
        disthist_in_roi = disthist_in_roi ./ sum(InROI_idx) * 100;
        disthist_labels = 0:0.1:disthist_max_inROI;

        fade_array = linspace(1,0,ProcSettings.FurthestFriendID)';
        fade_r2b = horzcat(fade_array,zeros(ProcSettings.FurthestFriendID,1),flipud(fade_array));

        disthistoverlay_inROI = figure('Visible',FigureVisibility);
        p01 = plot(disthist_labels,disthist_in_roi(1,:),'-','Color',fade_r2b(1,:));
        hold on
        p02 = plot(disthist_labels,disthist_in_roi(2,:),':','Color',fade_r2b(2,:));
        for p=3:ProcSettings.FurthestFriendID-1
            plot(disthist_labels,disthist_in_roi(p,:),':','Color',fade_r2b(p,:));
        end
        p_end = plot(disthist_labels,disthist_in_roi(ProcSettings.FurthestFriendID,:),'-','Color',fade_r2b(10,:));
        title([ProcSettings.ExptTitle,' : Distances to n^{th} Neighbour (within ROI)']);
        xlabel('Distance to n^{th} nearest neighbour (nm)');
        ylabel('Frequency');
        legend([p01 p02 p_end],{'n=1',['n=2-',num2str(ProcSettings.FurthestFriendID-1)],['n=',num2str(ProcSettings.FurthestFriendID),'']});
%         set(gca,'YLim',[0 ceil(max(max(disthist_in_roi(:,2:end-1)))/100)*100]);
        set(gca,'YLim',[0 4]);
        set(gca,'XLim',[0 50]); % need to calc the upper limit automagically
        set(gca,'TickDir','out');
        xlabel('Distance to n^{th} neighbour (nm)');
        ylabel('Percent of Events');
        
        SaveFileName = [FileList{FileID,2},' - Distances to nth Neighbour (within ROI)'];
        print(disthistoverlay_inROI,'-dpng','-r300',fullfile(pwd,'Histograms','Distances to Nth NN',SaveFileName));

        fprintf('\t%s','(done)');
        fprintf('\n');
        
    % Randomised in ROI
        
        fprintf('\t%s','NN Distance histograms (randomised events within ROI)'); 

        disthist_rnd_max = ceil(max(max(DistancesToN_rndm))/10)*10;
        disthist_rnd=[];
        for dn = 1:ProcSettings.FurthestFriendID
             disthist_rnd(dn,:) = hist(DistancesToN_rndm(dn,:),0:0.1:disthist_rnd_max);
        end
        disthist_rnd = disthist_rnd ./ sum(InROI_idx) * 100;
        disthist_rnd_labels = 0:0.1:disthist_rnd_max;
        
        fade_array = linspace(1,0,ProcSettings.FurthestFriendID)';
        fade_r2b = horzcat(fade_array,zeros(ProcSettings.FurthestFriendID,1),flipud(fade_array));
        
        disthistoverlay_rand_roi = figure('Visible',FigureVisibility);
        p01 = plot(disthist_rnd_labels,disthist_rnd(1,:),'-','Color',fade_r2b(1,:));
        hold on
        p02 = plot(disthist_rnd_labels,disthist_rnd(2,:),':','Color',fade_r2b(2,:));
        for p=3:ProcSettings.FurthestFriendID-1
            plot(disthist_rnd_labels,disthist_rnd(p,:),':','Color',fade_r2b(p,:));
        end
        p_end = plot(disthist_rnd_labels,disthist_rnd(ProcSettings.FurthestFriendID,:),'-','Color',fade_r2b(10,:));
        title([ProcSettings.ExptTitle,' : Distances to n^{th} Neighbour (randomised to ROI density)']);
        xlabel('Distance to n^{th} nearest neighbour (nm)');
        ylabel('Frequency');
        legend([p01 p02 p_end],{'n=1',['n=2-',num2str(ProcSettings.FurthestFriendID-1)],['n=',num2str(ProcSettings.FurthestFriendID),'']});
%         set(gca,'YLim',[0 ceil(max(max(disthist_rnd(:,2:end-1)))/100)*100]);
        set(gca,'YLim',[0 4]);
        set(gca,'XLim',[0 100]); % need to calc the upper limit automagically
        set(gca,'TickDir','out');
        xlabel('Distance to n^{th} neighbour (nm)');
        ylabel('Percent of Events');
        
        SaveFileName = [FileList{FileID,2},' - Distances to nth Neighbour (randomised to ROI density)'];
        print(disthistoverlay_rand_roi,'-dpng','-r300',fullfile(pwd,'Histograms','Distances to Nth NN',SaveFileName));

%                 % Plot the dNN for both random and real data on the same plot
% 
%                 %         % fit a poisson distribution to the randomised data at dNN
%                 %         [lambdahat,lambdaci] = poissfit(disthist_rnd(ProcSettings.FurthestFriendID,:),0.001); % fit with 99% CI
%                 %         % MLE = lambdahat
%                 %         % CIs = lambdaci
% 
%                         figure;
%                         % plot the randomised data for NN-FFID
%                         rand_plot = bar(disthist_rnd_labels,disthist_rnd(ProcSettings.FurthestFriendID,:),1.0,'FaceColor',[1 0 0],'EdgeColor','none');
%                         pH = arrayfun(@(x) allchild(x),rand_plot);
%                         set(pH,'FaceAlpha',0.1,'EdgeAlpha',0.1);
%                         hold on
% 
%                         % Fit a normal probability distribution function to the randomised data
%                         pd = fitdist(DistancesToN_rndm(10,:)','Kernel','Kernel','epanechnikov');
% 
%                         % Calculate values for plotting
%                         pdfX = pdf(pd,disthist_rnd_labels);
%                 %       plot(disthist_rnd_labels,pdfX,'LineWidth',2,'Color',[0 1 0.5])    % Plot the PDF function as-is
% 
%                         % Normalize the PDF to match the real-data histogram (for overlay)
%                         rangex = max(disthist_rnd_labels) - min(disthist_rnd_labels);       % Finds the range of this data.
%                         binwidth = rangex/size(disthist_rnd_labels,2);                      % Finds the width of each bin.
%                         plotarea = max(disthist_rnd(10,:)) * binwidth;                      % Find the data plot area
%                         pdfarea = max(pdfX) * binwidth;                                     % Find the PDF plot area
%                         y = (plotarea / pdfarea) * pdf(pd,disthist_rnd_labels);             % Expand PDF plot to match data area
%                         hh1 = plot(disthist_rnd_labels,y,'r-','LineWidth',2);               % plot the fitted line
% 
%                         % plot the real data
%                 %         test2 = bar(disthist_labels,disthist(10,:),'-','FillColor',[1 0 0]);
%                         data_plot = bar(disthist_labels,disthist_in_roi(ProcSettings.FurthestFriendID,:),1.0,'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
%                         pH = arrayfun(@(x) allchild(x),data_plot);
%                         set(pH,'FaceAlpha',0.1,'EdgeAlpha',0.1);
% 
%                         axis([0 100 0 1]);
% 
%                         % Normalise the real data histogram to the randomised PDF
%                         disthist_corr = pdf(pd,disthist_labels);
%                         disthist_corr = disthist_corr / max(disthist_corr);
%                 %          test3 = plot(disthist_labels,disthist_corr,'-','Color',[0.5 0.5 0]);
%                         disthist_corrected(ProcSettings.FurthestFriendID,:) = disthist(ProcSettings.FurthestFriendID,:) - (disthist(ProcSettings.FurthestFriendID,:) .* disthist_corr);
%                         test4 = plot(disthist_labels,disthist_corrected(ProcSettings.FurthestFriendID,:),'b-','LineWidth',2);
% 
% 
% 
%                 %+=+=+=+=+=+=+= Recreate the PDF for the AllFD        
%                         DistancesToAFD_rnd_inroi = sum(DistancesToN_inROI,1);
% 
%                         AFDdisthist_rnd_max = ceil(max(max(DistancesToAFD_rnd_inroi))/10)*10;
%                         AFDdisthist_rnd = hist(DistancesToAFD_rnd_inroi,0:0.1:AFDdisthist_rnd_max);
%                         AFDdisthist_rnd = AFDdisthist_rnd ./ sum(InROI_idx) * 100;
%                         AFDdisthist_rnd_labels = 0:0.1:AFDdisthist_rnd_max;
% 
%                         figure;
%                         rand_plot2 = bar(AFDdisthist_rnd_labels,AFDdisthist_rnd,1.0,'FaceColor',[1 0 0],'EdgeColor','none');
%                         pH = arrayfun(@(x) allchild(x),rand_plot2);
%                         set(pH,'FaceAlpha',0.1,'EdgeAlpha',0.1);
%                         axis([0 600 0 0.150]);
%                         hold on
%                         % plot the real data
%                         [sumdistanceplot(2,:),sumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.AFDCol),0:0.1:AFDdisthist_rnd_max);
%                         sumdistanceplot(2,:) = sumdistanceplot(2,:) ./ sum(total_events_table) * 100;
%                         data_plot = bar(sumdistanceplot(1,:),sumdistanceplot(2,:),1.0,'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
%                         pH = arrayfun(@(x) allchild(x),data_plot);
%                         set(pH,'FaceAlpha',0.1,'EdgeAlpha',0.1);
%                 %         
%                 %         
%                 %         % Fit a probability distribution function to the randomised data
%                 %         figure
%                 %         histfit(DistancesToAFD_rnd_inroi,1000,'gamma');
%                 %         pd = fitdist(AFDdisthist_rnd','Normal'); %,'Kernel','Kernel','epanechnikov');
%                 %         
%                 %         % Calculate values for plotting
%                 %         pdfX = pdf(pd,AFDdisthist_rnd_labels);
%                 % %       plot(disthist_rnd_labels,pdfX,'LineWidth',2,'Color',[0 1 0.5])    % Plot the PDF function as-is
%                 % 
%                 %         % Normalize the PDF to match the real-data histogram (for overlay)
%                 %         rangex = max(AFDdisthist_rnd_labels) - min(AFDdisthist_rnd_labels);       % Finds the range of this data.
%                 %         binwidth = rangex/size(AFDdisthist_rnd_labels,2);                      % Finds the width of each bin.
%                 %         plotarea = max(AFDdisthist_rnd) * binwidth;                      % Find the data plot area
%                 %         pdfarea = max(pdfX) * binwidth;                                     % Find the PDF plot area
%                 %         y = (plotarea / pdfarea) * pdf(pd,AFDdisthist_rnd_labels);             % Expand PDF plot to match data area
%                 %         hh1 = plot(AFDdisthist_rnd_labels,y,'g-','LineWidth',2);               % plot the fitted line
%                 % 
%                 %         
%                 %         figure
%                 %         plot(AFDdisthist_rnd_labels,pdfX,'LineWidth',2,'Color',[0 1 0.5])  % PDF as it is
%                 %         
%                 %         
%                 %         
%                 %         % normalise the real data itself
%                 %         distance_penalty = pdf(pd,data_proc(:,ProcSettings.AFDCol));
%                 %         
%                 %         modified_pts = distance_penalty~=0;
%                 %         
%                 %         [modified_pts_fig, modified_pts_axes] = DoMeAFigure(ProcSettings.AxisLimits);
%                 %         set(modified_pts_fig,'visible','on');   
%                 %         modified_ptsPlot = scatter(data_proc(modified_pts,1),data_proc(modified_pts,2),2,'.b');
%                 %         hold on
%                 %         unmodified_ptsPlot = scatter(data_proc(~modified_pts,1),data_proc(~modified_pts,2),2,'.r');
% 
% 
%                 %+=+=+=+=+=+=+=         

        clear dn p p01 p02 p_end  fade_array fade_r2b disthist_labels
        
%         close(disthistoverlay);
%         clear disthist  disthist_max fade_array fade_r2b disthistoverlay disthist_max

        close(disthistoverlay_inROI);
        clear disthistoverlay_inROI disthist_max_inROI disthist_in_roi

        close(disthistoverlay_rand_roi);
        clear disthist_rnd disthist_rnd_labels disthist_rnd_max  disthistoverlay_rand_roi
        
        fprintf('\t%s','(done)');
        fprintf('\n');
        
%====== Histograms of the single distance to the Nth NN

    % create the figure to hold all the distance histograms
    plotcollection = figure('Visible',FigureVisibility);        
        
    if ProcSettings.DoStrictdNN
		
    % for the data
    
        %Bin the distances to find the most common distance range.        
        distanceplot = [];
        xbins = 0:1:200; % max(data(:,3)); % bin width = 10 (nm)
        [distanceplot(2,:),distanceplot(1,:)] = hist(data_proc(:,ProcSettings.FFDCol),xbins);
        subplot(2,3,1);
        bar(distanceplot(1,:),distanceplot(2,:));
        axis tight
        set(gca,'TickDir','out');
        xlabel(['d_{',num2str(ProcSettings.FurthestFriendID),'} (nm)']);
        ylabel('Freq.');
        title(['d_{',num2str(ProcSettings.FurthestFriendID),'}']);
        normaldistance = distanceplot(1,distanceplot(2,:)==max(distanceplot(2,2:end-1))); % ignoring first and last bins
        axis([0 distanceplot(1,end) 0 1.1*distanceplot(2,distanceplot(1,:) == normaldistance)]);
        hold on
            distanceplot_rnd = [];
            [distanceplot_rnd(2,:),distanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.FFDCol),xbins);
            bar(distanceplot_rnd(1,:),distanceplot_rnd(2,:),'FaceColor',[1 0.8 0.8],'EdgeColor','none');
            normaldistance_rnd = distanceplot_rnd(1,distanceplot_rnd(2,:)==max(distanceplot_rnd(2,2:end-1))); % ignoring first and last bins
%         axis([0 distanceplot_rnd(1,end) 0 1.1*distanceplot_rnd(2,distanceplot_rnd(1,:) == normaldistance_rnd)]);

        % Invert distances to n10th in data col 5
        ProcSettings.invFFDCol = size(data_proc,2) + 1;
        ProcTblHeaders = horzcat(ProcTblHeaders,{['Norm''d FFD(',num2str(ProcSettings.FurthestFriendID),')']});

        data_proc(:,ProcSettings.invFFDCol) = normaldistance ./ data_proc(:,ProcSettings.FFDCol);
        data_proc_rndm(:,ProcSettings.invFFDCol) = normaldistance_rnd ./ data_proc_rndm(:,ProcSettings.FFDCol);

        dindexplot = [];
        xbins = 0:0.02:floor(max(data_proc(:,ProcSettings.invFFDCol)));
        [dindexplot(2,:),dindexplot(1,:)] = hist(data_proc(:,ProcSettings.invFFDCol),xbins);
        subplot(2,3,4);
        bar(dindexplot(1,:),dindexplot(2,:));
        axis tight
        set(gca,'TickDir','out');
        xlabel(['Inv. d_{',num2str(ProcSettings.FurthestFriendID),'}']);
        ylabel('Freq.');
        title(['Max d_{',num2str(ProcSettings.FurthestFriendID),'}/d_{',num2str(ProcSettings.FurthestFriendID),'}']);
        axis([0 2 0 ceil(1.1*max(dindexplot(2,:)))]);
        hold on
            dindexplot_rnd = [];
            [dindexplot_rnd(2,:),dindexplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.invFFDCol),xbins);
            bar(dindexplot_rnd(1,:),dindexplot_rnd(2,:),'FaceColor',[1 0.8 0.8],'EdgeColor','none');
    
        
    %====== Colour map of the single distance to the Nth NN

        if ProcSettings.SaveClusterColourImages
            fprintf('\t%s ', 'Single distance to Nth NN');             
            
            peakdindex = dindexplot(1,dindexplot(2,:)==max(dindexplot(2,2:end-1))); % ignoring first and last bins
            axis([0 dindexplot(1,end) 0 1.1*dindexplot(2,dindexplot(1,:) == peakdindex)]);
            ColourmapMax = peakdindex;
            
            data_randz_invFFD = data_proc(data_rand_ind,ProcSettings.invFFDCol); % randomised z data using rand idx calculated in the beginning
        
            [fig_invFFD, axes_invFFD] = DoMeAFigure(ProcSettings.AxisLimits);
            scatter(data_randx,data_randy,1,data_randz_invFFD,'.');
            caxis([0 ColourmapMax]);
            c=colormap;
            SaveFileName = [FileList{FileID,2},' - fig_invFFD'];
            SaveFolder = 'A - Distance to NN';
            print(fig_invFFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
            SaveColourBar(c,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar.png']));
            
            close(fig_invFFD)
            clear fig_invFFD axes_invFFD data_randz_invFFD c
                        
            data_randz_FFD = data_proc(data_rand_ind,ProcSettings.FFDCol); % randomised z data using rand idx calculated in the beginning
        
            [fig_FFD, axes_FFD] = DoMeAFigure(ProcSettings.AxisLimits);
            scatter(data_randx,data_randy,1,data_randz_FFD,'.');
            c=colormap;
            SaveFileName = [FileList{FileID,2},' - fig_FFD'];
            SaveFolder = 'A - Distance to NN';
            print(fig_FFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
            SaveColourBar(c,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar.png']));
            
            close(fig_FFD)
            clear fig_FFD axes_FFD data_randz_FFD
            
        end
        
        clear data_randz_invFFD SaveFileName peakdindex ColourmapMax dindexplot xbins normaldistance distanceplot SaveFolder dindexplot_rnd distanceplot_rnd c normaldistance_rnd
        fprintf('\t%s','(done)');
        fprintf('\n');
        
    end
    
%====== Histograms of the sum of the distances to the Nth NN

    if ProcSettings.DoSumdNN

        % send the focus back to the plots
        set(0, 'CurrentFigure', plotcollection);

        % Sum of distances to all NN up to Nth NN
        sumdistanceplot = [];
        xbins = 0:5:1000; % bin width = 10 (nm)
        [sumdistanceplot(2,:),sumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.AFDCol),xbins);
        subplot(2,3,2);
        bar(sumdistanceplot(1,:),sumdistanceplot(2,:));
        axis tight
        set(gca,'TickDir','out');
        xlabel(['Sum d_{',num2str(ProcSettings.FurthestFriendID),'} (nm)']);
        ylabel('Freq.');
        title(['Sum d_{',num2str(ProcSettings.FurthestFriendID),'}']);
        normalsumdistance = sumdistanceplot(1,sumdistanceplot(2,:)==max(sumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([0 sumdistanceplot(1,end) 0 1.1*sumdistanceplot(2,sumdistanceplot(1,:) == normalsumdistance(1))]);
        hold on
            sumdistanceplot_rnd = [];
            [sumdistanceplot_rnd(2,:),sumdistanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.AFDCol),xbins);
            bar(sumdistanceplot_rnd(1,:),sumdistanceplot_rnd(2,:),'FaceColor',[1 0.8 0.8],'EdgeColor','none');
            normalsumdistance_rnd = sumdistanceplot_rnd(1,sumdistanceplot_rnd(2,:)==max(sumdistanceplot_rnd(2,2:end-1))); % ignoring first and last bins
%             axis([0 sumdistanceplot_rnd(1,end) 0 1.1*sumdistanceplot_rnd(2,sumdistanceplot_rnd(1,:) == normalsumdistance_rnd)]);


        % Invert sum of distances to Nth NN
        ProcSettings.invAFDCol = size(data_proc,2) + 1;
        ProcTblHeaders = horzcat(ProcTblHeaders,{['Norm''d AFD(',num2str(ProcSettings.FurthestFriendID),')']});

        data_proc(:,ProcSettings.invAFDCol) = normalsumdistance(1) ./ data_proc(:,ProcSettings.AFDCol);
        data_proc_rndm(:,ProcSettings.invAFDCol) = normalsumdistance_rnd(1) ./ data_proc_rndm(:,ProcSettings.AFDCol);

        sindexplot = [];
        xbins = 0:0.01:2;
        [sindexplot(2,:),sindexplot(1,:)] = hist(data_proc(:,ProcSettings.invAFDCol),xbins);
        subplot(2,3,5);
        bar(sindexplot(1,:),sindexplot(2,:));
        axis tight
        set(gca,'TickDir','out');
        xlabel(['Inv. Sum d_{',num2str(ProcSettings.FurthestFriendID),'}']);
        ylabel('Freq.');
        title(['Max Sum d_{',num2str(ProcSettings.FurthestFriendID),'}/Sum d_{',num2str(ProcSettings.FurthestFriendID),'}']);
        axis([0 2 0 ceil(1.1*max(sindexplot(2,:)))]);
        hold on
            invsumdistanceplot_rnd = [];
            [invsumdistanceplot_rnd(2,:),invsumdistanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.invAFDCol),xbins);
            bar(invsumdistanceplot_rnd(1,:),invsumdistanceplot_rnd(2,:),'FaceColor',[1 0.8 0.8],'EdgeColor','none');
        

    %====== Colour map of the sum of the distances to the Nth NN

        if ProcSettings.SaveClusterColourImages
            
            fprintf('\t%s','Sum of distances to Nth NN'); 
            
            peaksindex = sindexplot(1,sindexplot(2,:)==max(sindexplot(2,2:end-1))); % ignoring first and last bins
            ColourmapMax = peaksindex;
        
            data_randz_invAFD = data_proc(data_rand_ind,ProcSettings.invAFDCol);%randomised z data

            [fig_invAFD, axes_invAFD] = DoMeAFigure(ProcSettings.AxisLimits);
            scatter(data_randx,data_randy,1,data_randz_invAFD,'.');
            caxis([0 ColourmapMax]);
            SaveFileName = [FileList{FileID,2},' - fig_invAFD'];
            SaveFolder = 'B - Sum of Distances to NN';
            print(fig_invAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
            SaveColourBar(colormap,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar.png']));

            close(fig_invAFD);
            clear fig_invAFD axes_invAFD data_randz_invAFD
            
            data_randz_AFD = data_proc(data_rand_ind,ProcSettings.AFDCol); % randomised z data using rand idx calculated in the beginning
        
%             ColourmapMax = normalsumdistance(1) ./ ColourmapMax;
            [fig_AFD, axes_AFD] = DoMeAFigure(ProcSettings.AxisLimits);
            scatter(data_randx,data_randy,1,data_randz_AFD,'.');
%             caxis([0 ColourmapMax]);
            SaveFileName = [FileList{FileID,2},' - fig_AFD'];
            SaveFolder = 'B - Sum of Distances to NN';
            print(fig_AFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
            SaveColourBar(colormap,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar.png']));
            
            close(fig_AFD)
            clear fig_AFD axes_AFD data_randz_AFD
            
        end
        
        clear data_randz_invAFD SaveFileName ColourmapMax peaksindex sindexplot xbins normalsumdistance sumdistanceplot SaveFolder invsumdistanceplot_rnd normalsumdistance_rnd sumdistanceplot_rnd
        fprintf('\t%s',' (done)');
        fprintf('\n');
        
    end
        
%====== Histograms of the log of the sum of the distances to the Nth NN

    if ProcSettings.DoInvLogSumdNN

        % Log of sum of distances to Nth NN
        ProcSettings.logAFDCol = size(data_proc,2) + 1;
        data_proc(:,ProcSettings.logAFDCol) = log(data_proc(:,ProcSettings.AFDCol));
        data_proc_rndm(:,ProcSettings.logAFDCol) = log(data_proc_rndm(:,ProcSettings.AFDCol));
        ProcTblHeaders = horzcat(ProcTblHeaders,{['Log AFD(',num2str(ProcSettings.FurthestFriendID),')']});
                
        % send the focus back to the plots
        set(0, 'CurrentFigure', plotcollection);

        logsumdistanceplot = [];
        xbins = floor(min(data_proc(:,ProcSettings.logAFDCol))):0.01:ceil(max(data_proc(:,ProcSettings.logAFDCol))); % bin width = 10 (nm)
        [logsumdistanceplot(2,:),logsumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.logAFDCol),xbins);
        subplot(2,3,3);
        bar(logsumdistanceplot(1,:),logsumdistanceplot(2,:));
        axis tight
        set(gca,'TickDir','out');
        xlabel(['Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
        ylabel('Freq.');
        title(['Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
        normallogsumdistance = logsumdistanceplot(1,logsumdistanceplot(2,:)==max(logsumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.logAFDCol))) ceil(max(data_proc(:,ProcSettings.logAFDCol))) 0 1.1*logsumdistanceplot(2,logsumdistanceplot(1,:) == normallogsumdistance)]);
        hold on
            logsumdistanceplot_rnd = [];    
            [logsumdistanceplot_rnd(2,:),logsumdistanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.logAFDCol),xbins);
            bar(logsumdistanceplot_rnd(1,:),logsumdistanceplot_rnd(2,:),'FaceColor',[1 0.8 0.8],'EdgeColor','none');
%   this isn't required here -->          normallogsumdistance_rnd = logsumdistanceplot_rnd(1,logsumdistanceplot_rnd(2,:)==max(logsumdistanceplot_rnd(2,2:end-1))); % ignoring first and last bins
%             axis([0 logsumdistanceplot_rnd(1,end) 0 1.1*logsumdistanceplot_rnd(2,logsumdistanceplot_rnd(1,:) == normallogsumdistance_rnd)]);

% Inverse log of sum of distances to Nth NN (i.e. low distances == high values)
        ProcSettings.invlogAFDCol = size(data_proc,2) + 1;
        ProcTblHeaders = horzcat(ProcTblHeaders,{['Norm''d Log AFD(',num2str(ProcSettings.FurthestFriendID),')']});

        MaxLogAFD = max(data_proc(:,ProcSettings.logAFDCol));
        data_proc(:,ProcSettings.invlogAFDCol) = MaxLogAFD - data_proc(:,ProcSettings.logAFDCol);
        
        MaxLogAFD_rnd = max(data_proc_rndm(:,ProcSettings.logAFDCol));
        data_proc_rndm(:,ProcSettings.invlogAFDCol) = MaxLogAFD_rnd - data_proc_rndm(:,ProcSettings.logAFDCol);
        
        invlogsumdistanceplot = [];
        xbins = floor(min(data_proc(:,ProcSettings.invlogAFDCol))):0.01:ceil(max(data_proc(:,ProcSettings.invlogAFDCol))); % bin width = 10 (nm)
        [invlogsumdistanceplot(2,:),invlogsumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.invlogAFDCol),xbins);
        subplot(2,3,6);
        bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:));
        axis tight
        set(gca,'TickDir','out');
        xlabel(['Inv. Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
        ylabel('Freq.');
        title(['Max Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})-Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
        normallogsumdistance2 = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.invlogAFDCol))) ceil(max(data_proc(:,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == normallogsumdistance2)]);
        hold on
            invlogsumdistanceplot_rnd = [];    
            [invlogsumdistanceplot_rnd(2,:),invlogsumdistanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.invlogAFDCol),xbins);
            bar(invlogsumdistanceplot_rnd(1,:),invlogsumdistanceplot_rnd(2,:),'FaceColor',[1 0.8 0.8],'EdgeColor','none');

    %====== Colour map of the inv log sum of the distances to the Nth NN

        fprintf('\t%s ', 'Inv log sum of distances to Nth NN'); 
        
        data_randz_invlogAFD = data_proc(data_rand_ind,ProcSettings.invlogAFDCol);
        
        if ProcSettings.SaveClusterColourImages
            [fig_invlogAFD, axes_invlogAFD] = DoMeAFigure(ProcSettings.AxisLimits);
            scatter(data_randx,data_randy,1,data_randz_invlogAFD,'.');
            SaveFileName = [FileList{FileID,2},' - fig_invlogAFD (Bright Background)'];
            SaveFolder = 'C - Inv Log Sum of Distances to NN';
            print(fig_invlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
            c = colormap;
            uncorrcmapmax = caxis;
            SaveColourBar(c,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar.png'])); % save current colormap
            set(axes_invlogAFD,'Color',c(1,:));
            set(fig_invlogAFD,'Color',c(1,:),'InvertHardcopy','off')
            SaveFileName = [FileList{FileID,2},' - fig_invlogAFD (Dark Background)'];
            print(fig_invlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
        
            close(fig_invlogAFD); 
            clear fig_invlogAFD axes_invlogAFD data_randz_invlogAFD
        end
               
        clear data_randz_invlogAFD SaveFileName c normallogsumdistance xbins SaveFolder normallogsumdistance2 c
        fprintf('\t%s','(done)');
        fprintf('\n');
        
    

                %====== Histograms of the random-normalised log of the sum of the distances to the Nth NN

                            if ProcSettings.NormaliseToDensity

                            % Find the peak of the ROI-randomised data
                                AFDdisthist_rnd_max = ceil(max(max(data_proc_rndm(:,ProcSettings.AFDCol)))/10)*10;
                                AFDdisthist_rnd = hist(data_proc_rndm(:,ProcSettings.AFDCol),0:0.1:AFDdisthist_rnd_max);
                                AFDdisthist_rnd = (AFDdisthist_rnd ./ sum(RndInROI_idx)) * 100;
                                AFDdisthist_rnd_labels = 0:0.1:AFDdisthist_rnd_max;

                            % Plot AFD for randomised and real-data
                                plot_correction_collection = figure;
                                subplot(1,2,1);
                                rand_plot2 = bar(AFDdisthist_rnd_labels,AFDdisthist_rnd,1.0,'FaceColor',[1 0 0],'EdgeColor','none');
                                pH = arrayfun(@(x) allchild(x),rand_plot2);
                                set(pH,'FaceAlpha',0.1,'EdgeAlpha',0.1);
                                hold on
                            % plot the real data
                                [sumdistanceplot_data(2,:),sumdistanceplot_data(1,:)] = hist(data_proc(InROI_idx,ProcSettings.AFDCol),0:0.1:AFDdisthist_rnd_max);
                                sumdistanceplot_data(2,:) = sumdistanceplot_data(2,:) ./ sum(total_events_table) * 100;
                                data_plot = bar(sumdistanceplot_data(1,:),sumdistanceplot_data(2,:),1.0,'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
                                pH = arrayfun(@(x) allchild(x),data_plot);
                                set(pH,'FaceAlpha',0.1,'EdgeAlpha',0.1);
                                
                            % set comfortable axis limits
                                axis tight
                                tmp_x_axis_max = max([floor(0.8*max(AFDdisthist_rnd_labels)/50)*50 , floor(0.7*max(sumdistanceplot_data(1,2:end-1))/50)*50]);
                                tmp_y_axis_max = max([floor(1.2*max(AFDdisthist_rnd)*50)/50 , floor(1.2*max(sumdistanceplot_data(2,2:end-1))*50)/50]);
                                axis([0 tmp_x_axis_max 0 tmp_y_axis_max]);
                                clear tmp_x_axis_max tmp_y_axis_max
                                
                                
                             % Option A - Fit a probability distribution function to the randomised data
%                                 histfit(data_rndm_proc(:,ProcSettings.AFDCol),1000,'gamma');
%                                 set(gca,'YLim',[0 2000]);
%                                 set(gca,'XLim',[0 600]);
% 
%                            % Option B - using normfit
%                                 [mu,sigma,muci,sigmaci] = normfit(data_proc_rndm(:,ProcSettings.AFDCol));
%                                 data_rnd_normdist = normpdf(AFDdisthist_rnd_labels,mu,sigma);
%                                 figure
%                                 plot(AFDdisthist_rnd_labels,data_rnd_normdist,'LineWidth',2,'Color',[0 1 0.5])  % PDF as it is
%                                 % Normalize the PDF to match the real-data histogram (for overlay)
%                                 rangex = max(AFDdisthist_rnd_labels) - min(AFDdisthist_rnd_labels);       % Finds the range of this data.
%                                 binwidth = rangex/size(AFDdisthist_rnd_labels,2);                      % Finds the width of each bin.
%                                 plotarea = max(AFDdisthist_rnd) * binwidth;                      % Find the data plot area
%                                 pdfarea = max(data_rnd_normdist) * binwidth;                            % Find the PDF plot area
%                                 y = (plotarea / pdfarea) * normpdf(AFDdisthist_rnd_labels,mu,sigma);           % Expand PDF plot to match data area
%                                 plot(AFDdisthist_rnd_labels,y,'g-','LineWidth',2);                      % plot the fitted line
                                
                                
                            % Option C - using an Epanechnikov kernel dist fit
                                pdistn2 = fitdist(data_proc_rndm(:,ProcSettings.AFDCol),'Kernel','Kernel','epanechnikov');
                                data_rnd_normdist2 = pdf(pdistn2,AFDdisthist_rnd_labels); % Calculate values for plotting a different PDF
%                                 plot(AFDdisthist_rnd_labels,data_rnd_normdist2,'LineWidth',2,'Color',[0 1 0.5])  % PDF as it is

                            % Normalize the PDF to match the real-data histogram (for overlay)
                                rangex = max(AFDdisthist_rnd_labels) - min(AFDdisthist_rnd_labels);       % Finds the range of this data.
                                binwidth = rangex/size(AFDdisthist_rnd_labels,2);                      % Finds the width of each bin.
                                plotarea = max(AFDdisthist_rnd) * binwidth;                      % Find the data plot area

                                pdfarea2 = max(data_rnd_normdist2) * binwidth;                          % Find the PDF plot area
                                y = (plotarea / pdfarea2) * pdf(pdistn2,AFDdisthist_rnd_labels);        % Expand PDF plot to match data area
                                hold on
                                plot(AFDdisthist_rnd_labels,y,'g-','LineWidth',2);                      % plot the fitted line
                                set(gca,'TickDir','out');
                                xlabel(['Sum distance to NN_{',num2str(ProcSettings.FurthestFriendID),'}']);
                                ylabel('% events');
                                title(['Sum of distances to NN_{',num2str(ProcSettings.FurthestFriendID),'}']);

                            %Find the peak of the fitted curve and convert it to an ILSD value
                                [Epanechnikov_peak_val, Epanechnikov_peak_idx] = max(data_rnd_normdist2);
                                Epanechnikov_peak_AFD = AFDdisthist_rnd_labels(1,Epanechnikov_peak_idx);
                                RandInROI_CorrFactor = MaxLogAFD_rnd - log(Epanechnikov_peak_AFD);
                                
                            % Offset the Inverse Log of sum of distances to Nth NN ot the peak of the randomised data
                                ProcSettings.CorrinvlogAFDCol = size(data_proc,2) + 1;
                                ProcTblHeaders = horzcat(ProcTblHeaders,{['Corrected inv Log AFD(',num2str(ProcSettings.FurthestFriendID),')']});
                                
                                data_proc(:,ProcSettings.CorrinvlogAFDCol) = data_proc(:,ProcSettings.invlogAFDCol) - RandInROI_CorrFactor;
                                data_proc_rndm(:,ProcSettings.CorrinvlogAFDCol) = data_proc_rndm(:,ProcSettings.invlogAFDCol) - RandInROI_CorrFactor;

                                subplot(1,2,2);
                                corrinvlogsumdistanceplot = [];
                                xbins = floor(min(data_proc(:,ProcSettings.CorrinvlogAFDCol))):0.01:ceil(max(data_proc(:,ProcSettings.CorrinvlogAFDCol))); % bin width = 10 (nm)
                                [corrinvlogsumdistanceplot(2,:),corrinvlogsumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.CorrinvlogAFDCol),xbins);
                                bar(corrinvlogsumdistanceplot(1,:),corrinvlogsumdistanceplot(2,:));
                                axis tight
                                set(gca,'TickDir','out');
                                xlabel(['Corrected inv.Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
                                ylabel('Freq.');
                                title(['Corrected inv Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
                                normalcorrinvlogsumdistance = corrinvlogsumdistanceplot(1,corrinvlogsumdistanceplot(2,:)==max(corrinvlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
                                axis([floor(min(data_proc(:,ProcSettings.CorrinvlogAFDCol))) ceil(max(data_proc(:,ProcSettings.CorrinvlogAFDCol))) 0 1.1*corrinvlogsumdistanceplot(2,corrinvlogsumdistanceplot(1,:) == normalcorrinvlogsumdistance)]);
                                hold on
                                    corrinvlogsumdistanceplot_rnd = [];    
                                    [corrinvlogsumdistanceplot_rnd(2,:),corrinvlogsumdistanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.CorrinvlogAFDCol),50);
                                    bar(corrinvlogsumdistanceplot_rnd(1,:),corrinvlogsumdistanceplot_rnd(2,:),'FaceColor',[1 0.8 0.8],'EdgeColor','none');
                                set(gca,'XLim',([floor(min(data_proc_rndm(:,ProcSettings.CorrinvlogAFDCol))) ceil(max(data_proc(:,ProcSettings.CorrinvlogAFDCol)))])); % fix axis to show the corrected histogram positioning including the randomised data
                                
                                SaveFileName = [FileList{FileID,2},' - NN correction histograms to randomised in ROI'];
                                print(plot_correction_collection,'-dpng','-r600',fullfile(pwd,'Histograms','Collections',[SaveFileName,'.png']));
                                close(plot_correction_collection)
                                clear plot_correction_collection xbins Epanechnikov_peak_AFD pdfarea2 y pdistn2 ...
                                      data_rnd_normdist2 pdfarea2 rangex binwidth plotarea corrinvlogsumdistanceplot_rnd ...
                                      corrinvlogsumdistanceplot AFDdisthist_rnd_labels AFDdisthist_rnd_max AFDdisthist_rnd ...
                                      data_plot Epanechnikov_peak_idx Epanechnikov_peak_val rand_plot2 pH normalcorrinvlogsumdistance ...
                                      sumdistanceplot_data MaxLogAFD MaxLogAFD_rnd 
                                  
                            %====== Colour map of the randomised-corrected inv log sum of the distances to the Nth NN

                                fprintf('\t%s ', 'Inv log sum of distances to Nth NN, normalised to random'); 

                                if ProcSettings.SaveClusterColourImages
                                    
                                % save the corrected image for the ROI data
 
                                    data_randz_corrinvlogAFD = data_proc(data_rand_ind,ProcSettings.CorrinvlogAFDCol);
                                    
                                    [fig_corrinvlogAFD, axes_corrinvlogAFD] = DoMeAFigure(ProcSettings.AxisLimits);
                                    scatter(data_randx,data_randy,1,data_randz_corrinvlogAFD,'.');
                                    corr_caxis = caxis;
                                    corr_caxis(1) = 0;
                                    caxis(corr_caxis);
                                    SaveFileName = [FileList{FileID,2},' - fig_corr_invlogAFD (Bright Background)'];
                                    SaveFolder = 'C - Inv Log Sum of Distances to NN';
                                    print(fig_corrinvlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
                                    
                                    c = colormap;                           
                                    set(axes_corrinvlogAFD,'Color',c(1,:));
                                    set(fig_corrinvlogAFD,'Color',c(1,:),'InvertHardcopy','off')
                                    SaveFileName = [FileList{FileID,2},' - fig_corr_invlogAFD (Dark Background)'];
                                    print(fig_corrinvlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));

                                    SaveColourBar(c,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar.png'])); % save current colormap

%                                % match the colour scaling of the uncorrected image - to check that it has actually adjusted something...
%                                     caxis(uncorrcmapmax); 
%                                     SaveFileName = [FileList{FileID,2},' - fig_corr_invlogAFD (Dark Background, matched colour scale to uncorr)'];
%                                     print(fig_corrinvlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
%                                     c = colormap;
%                                     SaveColourBar(c,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar(matched to uncorr).png'])); % save current colormap
%                                     clear uncorrcmapmax

                                    close(fig_corrinvlogAFD);
                                    clear fig_corrinvlogAFD axes_corrinvlogAFD data_randz_corrinvlogAFD c data_randz_corrinvlogAFD
                                    
                                    
                                % save the corrected image for the ROI randomised data
%                                     data_randz_corrinvlogAFD_rndm = ;
                                    
                                    [fig_corrinvlogAFD_rnd, axes_corrinvlogAFD_rnd] = DoMeAFigure(ProcSettings.AxisLimits);
%                                     set(gcf,'visible','on');
                                    scatter(data_proc_rndm(RndInROI_idx,ProcSettings.xCol_procd),data_proc_rndm(RndInROI_idx,ProcSettings.yCol_procd),1,data_proc_rndm(RndInROI_idx,ProcSettings.CorrinvlogAFDCol),'.');
%                                     scatter(data_randx,data_randy,1,data_randz_corrinvlogAFD_rndm,'.',data_proc_rndm(data_rand_ind,ProcSettings.CorrinvlogAFDCol));
                                    caxis(corr_caxis); % match colour axis of randomised data to match that of the real data
                                    SaveFileName = [FileList{FileID,2},' - fig_corr_invlogAFD (Randomised - Bright Background)'];
                                    SaveFolder = 'C - Inv Log Sum of Distances to NN';
                                    print(fig_corrinvlogAFD_rnd,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
                                    
                                    c = colormap;
                                    set(axes_corrinvlogAFD_rnd,'Color',c(1,:));
                                    set(fig_corrinvlogAFD_rnd,'Color',c(1,:),'InvertHardcopy','off')
                                    SaveFileName = [FileList{FileID,2},' - fig_corr_invlogAFD (Randomised - Dark Background)'];
                                    print(fig_corrinvlogAFD_rnd,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
                                    
                                    close(fig_corrinvlogAFD_rnd);
                                    clear fig_corrinvlogAFD_rnd axes_corrinvlogAFD_rnd data_randz_corrinvlogAFD c corr_caxis
                                    
                                end

                                clear data_randz_invlogAFD SaveFileName c normallogsumdistance normallogsumdistance2 xbins ...
                                      MaxLogAFD SaveFolder uncorrcmapmax data_randz_corrinvlogAFD_rndm  ...
                                      invlogsumdistanceplot_rnd logsumdistanceplot logsumdistanceplot_rnd logsumdistanceplot logsumdistanceplot_rnd
                                fprintf('\t%s','(done)');
                                fprintf('\n');

                            end

    end
    
%====== Save the histogram collection

        fprintf('\t%s','Histogram summaries'); 

        set(0, 'CurrentFigure', plotcollection);
        SaveFileName = [FileList{FileID,2},' - NN histogram collection'];
        print(plotcollection,'-dpng','-r600',fullfile(pwd,'Histograms','Collections',[SaveFileName,'.png']));
        close(plotcollection)
        clear plotcollection
        
    % new histogram plot for InvLogSumDistNN alone
        fig_plot_invlogAFD = figure('Visible','off');
        bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:));
        axis tight
        set(gca,'TickDir','out');
        title([ProcSettings.ExptTitle,' : Inv. Log of Sum of distances to NN_{1}–NN_{',num2str(ProcSettings.FurthestFriendID),'}']);
        normallogsumdistance2 = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.invlogAFDCol))) ceil(max(data_proc(:,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == normallogsumdistance2)]);
        SaveFileName = [FileList{FileID,2},' - InvLogSumDist-NN histogram'];
        print(fig_plot_invlogAFD,'-dpng','-r600',fullfile(pwd,'Histograms','ILSD-NN Histogram',[SaveFileName,'.png']));
        close(fig_plot_invlogAFD)
        clear fig_plot_invlogAFD normallogsumdistance2 invlogsumdistanceplot
        
        fprintf('\t%s', '(done)');
        fprintf('\n');

%% Thresholding

        if ProcSettings.DoThresholds
            
            % create containers for threshold summaries
            stats_hhmx = zeros(4,size(FileList,1));     % row 1 - threshold value for that table's method
            stats_mmsd = zeros(4,size(FileList,1));     % row 2 - total points over threshold
            stats_medn = zeros(4,size(FileList,1));     % row 3 - total points in image
            stats_hmx = zeros(4,size(FileList,1));      % row 4 - percent over threshold
            stats_mpsd = zeros(4,size(FileList,1));
            stats_thrfix = zeros(4,size(FileList,1));   % fixed-threshold stats

        %====== Adaptive thresholding - Half hist max

                fprintf('%s', 'Applying Thresholds:');
                fprintf('\n');
                fprintf('\t%s','Adaptive threshold > Half histogram max');

                Cluster_thr_halfhistmax = (invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,:))))/2;
                SaveToFolder = fullfile('D - Threshold','Half histogram peak');
                SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_halfhistmax),' (Half hist peak)'];
                RefDataCol = ProcSettings.invlogAFDCol;

                [halfhistmax_pts_over, halfhistmax_pts_under] = thresholder(data_proc,Cluster_thr_halfhistmax,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
                num_events_over = size(halfhistmax_pts_over,1);

            % update the summary table
                stats_hhmx(1,FileID) = Cluster_thr_halfhistmax;
                stats_hhmx(2,FileID) = num_events_over;
                stats_hhmx(3,FileID) = total_events_table;
                stats_hhmx(4,FileID) = (num_events_over/total_events_table) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');

        %====== Adaptive thresholding - Hist max

                fprintf('\t%s','Adaptive threshold > Histogram max');
                Cluster_thr_histmax = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,:)));
                SaveToFolder = fullfile('D - Threshold','Histogram peak');
                SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_histmax),' (Histogram peak)'];
                RefDataCol = ProcSettings.invlogAFDCol;

                [histmax_pts_over, histmax_pts_under] = thresholder(data_proc,Cluster_thr_histmax,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
                num_events_over = size(histmax_pts_over,1);

            % update the summary table
                stats_hmx(1,FileID) = Cluster_thr_histmax;
                stats_hmx(2,FileID) = num_events_over;
                stats_hmx(3,FileID) = total_events_table;
                stats_hmx(4,FileID) = (num_events_over/total_events_table) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');

        %====== Adaptive thresholding - Mean minus 1*StDev

                fprintf('\t%s','Adaptive threshold > Mean minus StDev');
                Cluster_thr_meanminus1stdev = mean(data_proc(:,ProcSettings.invlogAFDCol)) - std(data_proc(:,ProcSettings.invlogAFDCol));
                SaveToFolder = fullfile('D - Threshold','Mean minus StDev');
                SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_meanminus1stdev),' (Mean - StDev)'];
                RefDataCol = ProcSettings.invlogAFDCol;

                [meanminus1stdev_pts_over, meanminus1stdev_pts_under] = thresholder(data_proc,Cluster_thr_meanminus1stdev,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
                num_events_over = size(meanminus1stdev_pts_over,1);

            % update the summary table
                stats_mmsd(1,FileID) = Cluster_thr_meanminus1stdev;
                stats_mmsd(2,FileID) = num_events_over;
                stats_mmsd(3,FileID) = total_events_table;
                stats_mmsd(4,FileID) = (num_events_over/total_events_table) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');

        %====== Adaptive thresholding - Mean plus 1*StDev

                fprintf('\t%s','Adaptive threshold > Mean plus StDev');

                Cluster_thr_meanplus1stdev = mean(data_proc(:,ProcSettings.invlogAFDCol)) + std(data_proc(:,ProcSettings.invlogAFDCol));
                SaveToFolder = fullfile('D - Threshold','Mean plus StDev');
                SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_meanplus1stdev),' (Mean + StDev)'];
                RefDataCol = ProcSettings.invlogAFDCol;

                [meanplus1stdev_pts_over, meanplus1stdev_pts_under] = thresholder(data_proc,Cluster_thr_meanplus1stdev,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
                num_events_over = size(meanplus1stdev_pts_over,1);

            % update the summary table
                stats_mpsd(1,FileID) = Cluster_thr_meanplus1stdev;
                stats_mpsd(2,FileID) = num_events_over;
                stats_mpsd(3,FileID) = total_events_table;
                stats_mpsd(4,FileID) = (num_events_over/total_events_table) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');

        %====== Adaptive thresholding - Median

                fprintf('\t%s','Adaptive threshold > Median');

                Cluster_thr_median = median(data_proc(:,ProcSettings.invlogAFDCol));
                SaveToFolder = fullfile('D - Threshold','Median');
                SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_median),' (Median)'];
                RefDataCol = ProcSettings.invlogAFDCol;

                [median_pts_over, median_pts_under] = thresholder(data_proc,Cluster_thr_median,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
                num_events_over = size(median_pts_over,1);

            % update the summary table
                stats_medn(1,FileID) = Cluster_thr_median;
                stats_medn(2,FileID) = num_events_over;
                stats_medn(3,FileID) = total_events_table;
                stats_medn(4,FileID) = (num_events_over/total_events_table) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');

        %====== Fixed thresholding - as per ProcSettings

                fprintf('\t%s','User-specified threshold > ',num2str(ProcSettings.FixedThreshold));

                Cluster_thr_fixed = ProcSettings.FixedThreshold;
                SaveToFolder = fullfile('D - Threshold',['User-specified (',num2str(ProcSettings.FixedThreshold),')']);
                SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_fixed),' (user-specified)'];
                RefDataCol = ProcSettings.invlogAFDCol;

                [fixedthr_pts_over, fixedthr_pts_under] = thresholder(data_proc,Cluster_thr_fixed,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
                num_events_over = size(fixedthr_pts_over,1);

            % update the summary table
                stats_thrfix(1,FileID) = Cluster_thr_fixed;
                stats_thrfix(2,FileID) = num_events_over;
                stats_thrfix(3,FileID) = total_events_table;
                stats_thrfix(4,FileID) = (num_events_over/total_events_table) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');

        %====== Save Adaptive thresholding data

                HeaderFormat = '%s\t%s\t%s';
                fn_stats_save = [FileList{FileID,2},' - Threshold Data.txt'];
                fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_stats_save),'w');

                fprintf(fid,HeaderFormat,'[Threshold Method]','[Threshold Value]','[% In Clusters]');
                fprintf(fid,'\r\n');

                DataFormat = '%s\t%.3f\t%.3f';

                fprintf(fid,DataFormat,'Half hist. peak',Cluster_thr_halfhistmax,stats_hhmx(4,FileID));
                fprintf(fid,'\r\n');

                fprintf(fid,DataFormat,'Histogram peak',Cluster_thr_histmax,stats_hmx(4,FileID));
                fprintf(fid,'\r\n');

                fprintf(fid,DataFormat,'Mean minus StDev',Cluster_thr_meanminus1stdev,stats_mmsd(4,FileID));
                fprintf(fid,'\r\n');

                fprintf(fid,DataFormat,'Mean plus StDev',Cluster_thr_meanplus1stdev,stats_mpsd(4,FileID));
                fprintf(fid,'\r\n');

                fprintf(fid,DataFormat,'Median',Cluster_thr_median,stats_medn(4,FileID));
                fprintf(fid,'\r\n');

                fprintf(fid,DataFormat,'Fixed',Cluster_thr_fixed,stats_thrfix(4,FileID));
                fprintf(fid,'\r\n');

                fclose(fid);
                clear fid
                
                fprintf('%s', '(done)');
                fprintf('\n');
        end

%% ====== Save processed data tables
        if ProcSettings.SaveMATFiles
            fprintf('%s', 'Saving data and summaries');

            stringy = '%s';
            for g = 1:(length(ProcTblHeaders)-1)
                stringy = strcat(stringy,'\t%s'); % replace \t with a comma for csv
            end

            fn_data_proc_save = [FileList{FileID,2},'_proc.txt'];
            fid = fopen(fullfile(pwd,'X - Tables',fn_data_proc_save),'w');
            fprintf(fid,stringy,ProcTblHeaders{:});
            fprintf(fid,'\r\n');
            fclose(fid);
            dlmwrite(fullfile(pwd,'X - Tables',fn_data_proc_save),data_proc,'-append','delimiter','\t');

            fnsave = [FileList{FileID,2},'.mat'];
            if ProcSettings.DoThresholds
                save(fullfile(pwd,'Z - MAT Files',fnsave),...
                    'data_proc',...
                    'ProcTblHeaders',...
                    'halfhistmax_pts_over',...
                    'halfhistmax_pts_under',...
                    'histmax_pts_over',...
                    'histmax_pts_under',...
                    'meanminus1stdev_pts_over',...
                    'meanminus1stdev_pts_under',...
                    'meanplus1stdev_pts_over',...
                    'meanplus1stdev_pts_under',...
                    'median_pts_over',...
                    'median_pts_under'...
                    );
            else
                save(fullfile(pwd,'Z - MAT Files',fnsave),...
                    'data_proc',...
                    'ProcTblHeaders'...
                    );
            end
            
            clear fid fn_data_proc_save g stringy
            
            fprintf('%s', '(done)');
            fprintf('\n');
        end
    
%====== Show progress tracking information

        ThisRegionTimestamp = toc(mainproc_stopwatch);

        InfoMessage =  [datestr(fix(clock),'HH:MM:SS'),9,'Completed processing Table ',num2str(FileID),'.'];
        disp(InfoMessage);
        
        InfoMessage =  [9,9,9,sprintf('%0.2f',((ThisRegionTimestamp - PreviousRegionTimestamp)/60)),' minutes, ',num2str(total_events_table),' events, ',num2str(floor(total_events_table/(ThisRegionTimestamp - PreviousRegionTimestamp))),' events per second.'];
        disp(InfoMessage);
        
        InfoMessage = [9,9,9,num2str(round(FileID / size(FileList,1) * 100)) '% processed.'];
        disp(InfoMessage);
        
        disp('---------------------------------------------------------');
        
        % Update the tracking info to reflect the newly finished region
        PreviousRegionTimestamp = ThisRegionTimestamp;

%         % clean up your mess
%         clear Cluster_thr_halfhistmax Cluster_thr_meanminus1stdev Cluster_thr_median ColourmapMax c data_density data_len data_proc data_randx data_randy data_randz data_rand_ind datatable dindexplot distanceplot DistancesToN disthistoverlay dn fid fn_data_proc_save fn_stats_save fnsave g HeaderFormat InfoMessage invlogsumdistanceplot logsumdistanceplot MaxLogAFD normaldistance normallogsumdistance normallogsumdistance2 normalsumdistance peakdindex peaksindex percent_in_halfhistmax percent_in_meanminusstdev percent_in_median plot_pts_in_halfhistmax plot_pts_in_meanminus1stdev plot_pts_in_median ProcTblHeaders sindexplot stringy sumdistanceplot total_events_region xbins
%         clear invlogsumdistanceplot DataEvents_fig fig_invAFD fig_invFFD fig_invlogAFD fig_plot_invlogAFD fig_plot_pts_in_halfhistmax fig_plot_pts_in_meanminus1stdev fig_plot_pts_in_median plotcollection

    end % end of this data table processing
    
%% Finish up - Overall Summaries

%====== Save Adaptive thresholding summary data

if ProcSettings.DoThresholds
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
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_hhmx,'-append','delimiter','\t','precision','%.3f');
    
    % hist max
    fn_summary_proc_save = 'Combined Thresholds - Histogram peak.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Histogram peak');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_hmx,'-append','delimiter','\t','precision','%.3f');
    
    % mean minus stdev
    fn_summary_proc_save = 'Combined Thresholds - Mean minus StDev.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Mean minus StDev');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_mmsd,'-append','delimiter','\t','precision','%.3f');

    % mean plus stdev
    fn_summary_proc_save = 'Combined Thresholds - Mean plus StDev.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Mean plus StDev');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_mpsd,'-append','delimiter','\t','precision','%.3f');
    
    % median
    fn_summary_proc_save = 'Combined Thresholds - Median.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Median');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_medn,'-append','delimiter','\t','precision','%.3f');    

    
    % user-specified
    fn_summary_proc_save = ['Combined Thresholds - User Specified (',num2str(ProcSettings.FixedThreshold),').txt'];
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - User Specified (',num2str(ProcSettings.FixedThreshold),')');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_thrfix,'-append','delimiter','\t','precision','%.3f');
    
    clear fid fn_summary_proc_save g stringy
end
%====== Summarise processing time
    
    ExecTime=toc(mainproc_stopwatch);
    InfoMessage=[datestr(fix(clock),'HH:MM:SS'),9,'Finished processing ' num2str(size(FileList,1)) ' data tables.'];
    disp(InfoMessage);

    InfoMessage=[9,9,9,'Total time: ' sprintf('%0.2f',(ExecTime/60)) ' minutes.'];
    disp(InfoMessage);

    InfoMessage=[9,9,9,'Average time: ' sprintf('%0.2f',(ExecTime/60)/size(FileList,1)) ' minutes per data table.'];
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
% scatter(data(:,ProcSettings.xCoordsColumn),data(:,ProcSettings.yCoordsColumn),10,data(:,ProcSettings.AFDCol),'.');
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
% scatter(data_randx,data_randy,1,data_randz,'.');
% % scatter(data(:,ProcSettings.xCoordsColumn),data(:,ProcSettings.yCoordsColumn),1,data(:,ProcSettings.invlogAFDCol),'.');
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
% scatter(data_proc(:,ProcSettings.xCoordsColumn),data_proc(:,ProcSettings.yCoordsColumn),3,data_proc(:,ProcSettings.clusteredCol),'.');
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

