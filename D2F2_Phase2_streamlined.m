%%
% This version includes a folder-list option for batch processing

ThisVersion = 0.5;
DoubleCheckAllSettings = true;
FigureVisibility = 'on';
load('CustomColorMap.mat'); % from the folder with this script in it. Gets around MATLABs constant tinkering with LUTs.

ColorPlots_Data = [0.34 0.34 0.97]; %Deep blue RGB 87 87 249 ... light blue RGB 192 192 255
ColorPlots_Rndm = [0.97 0.25 0.25]; %Deep red RGB 249 64 64 ... light red RGB 255 224 224

SimpleThresholdFactor = 0.75;

OverwriteExisting = false;  % Don't overwrite files that already exist; the saves time if you restart and don't want to wait while images are saved.


%% Begin

% housekeeping
home %clean up the command window
rng('shuffle')  % set the random seed
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary'); % disables a warning about a temporary variable in the parallel processing section. This is fine but warnings scare people so I am turning it off.
warning('off','MATLAB:MKDIR:DirectoryExists'); % don't warn about existing data folders

% initiate garbage collection to free up memory; this can help improve stability.
heapMaxMemory = java.lang.Runtime.getRuntime.maxMemory;
heapTotalMemory = java.lang.Runtime.getRuntime.totalMemory;
heapFreeMemory = java.lang.Runtime.getRuntime.freeMemory;
if(heapFreeMemory < (heapTotalMemory*0.01))
    java.lang.Runtime.getRuntime.gc;
end
clear heapTotalMemory heapFreeMemory

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
    if exist(fullfile(DataDirName,'Phase1FileList.mat'), 'file') ~= 0
        load(fullfile(DataDirName,'Phase1FileList.mat'));
        if strfind(pwd,OutputFolderName) ~=0
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
            errordlg('This folder doesn''t match the description in SDFN.mat');
            error('This folder doesn''t match the description in SDFN.mat');
        end

    else
        errordlg('Phase 1 folder information file does not exist! You need to perform Phase 1 operations.');
    end

else
    error('Cancelled?! So rude.');
end

% correct the 8 bit RGB colours to fractional
 if sum(ProcSettings.ColourEventsOver)+sum(ProcSettings.ColourEventsUnder)+sum(ProcSettings.ColourBackground) > 1
    ProcSettings.ColourEventsOver = ProcSettings.ColourEventsOver ./ 255;
    ProcSettings.ColourEventsUnder = ProcSettings.ColourEventsUnder ./ 255;
    ProcSettings.ColourBackground = ProcSettings.ColourBackground ./ 255;
 end
 
ProcSettings.SaveLowDPI = '-r600'; % for quicker saving of preview images
ProcSettings.SaveHighDPI = '-r1800'; % for quicker saving of preview images


if DoubleCheckAllSettings
    prompt = {'xCoord Column:','yCoord Column:','Channel Col:','Furthest Neighbour:','Output Image Scale (nm/px):'};
    dlg_title = 'Check settings to apply...';
    num_lines = 1;
    defaults = {num2str(ProcSettings.xCoordsColumn), ...
                num2str(ProcSettings.yCoordsColumn), ...
                num2str(ProcSettings.ChannelIDColumn), ...
                num2str(ProcSettings.FurthestFriendID), ...
                num2str(ProcSettings.ImageSize / str2double(strrep(ProcSettings.SaveHighDPI,'-r',''))) ...
                };

    answer = inputdlg(prompt,dlg_title,num_lines,defaults);

    if ~isempty(answer)
        ProcSettings.xCoordsColumn = str2double(answer(1,1));
        ProcSettings.yCoordsColumn = str2double(answer(2,1));
        ProcSettings.ChannelIDColumn = str2double(answer(3,1));
        ProcSettings.FurthestFriendID = str2double(answer(4,1));
        ProcSettings.SaveHighDPI = ['-r',num2str(ProcSettings.ImageSize / str2double(answer{5,1}))];
    else
        error('Cancelled?! So rude.');
    end

    DoubleCheckAllSettings = false; % we only want to do this once
    clear prompt dlg_title num_lines defaults answer
end 

% Display some information
disp('---------------------------------------------------------');
InfoMessage=[datestr(fix(clock),'HH:MM:SS'),9,'Processing ' num2str(size(FileList,1)) ' data tables...'];
disp(InfoMessage);
disp('---------------------------------------------------------');

% quick stats container for simple threshold
stats_corr2rand = zeros(6,size(FileList,1));   % threshold stats
stats_ekovcorr2rand = zeros(6,size(FileList,1));   % threshold stats
stats_thrfix = zeros(4,size(FileList,1));   % fixed-threshold stats
stats_thrtheo = zeros(4,size(FileList,1));   % fixed-threshold stats


%% Phase 2 - Load & check data, load pre-saved ROIs

PreviousRegionTimestamp = 0; % Initialise the timing tracker.
mainproc_stopwatch = tic; %start the clock to track processing time

for FileID = 1:size(FileList,1)

    CurrentTableFileName=FileList{FileID,1};
    InfoMessage =  [datestr(fix(clock),'HH:MM:SS'),9,'[',num2str(FileID),'/',num2str(size(FileList,1)),'] Opening file: ',CurrentTableFileName];
    disp(InfoMessage);
    clear InfoMessage

    Phase1Source_filename = [FileList{FileID,2},' - Phase1Files.mat'];
    load(fullfile(DataDirName,'Z - MAT Files',Phase1Source_filename));
    clear Phase1Source_filename
    
    % Update ProcSettings with new variables for this image
    ProcSettings.AxisLimits = ProcSettingsLocal.AxisLimits;
    ProcSettings.PlotWidth = ProcSettingsLocal.PlotWidth;
    ProcSettings.SaveHighDPI = ProcSettingsLocal.SaveHighDPI;
    ProcSettings.ExptTitle = ProcSettingsLocal.ExptTitle;

    % Events that lie exactly on the image boundaries regularly crash matlab during print() function
    % to get around this, any events lying exactly on the image boundaries are offset by 0.1 (nm)
    % It is anticipated that events this close to the edge are not useful events anyway in SMLM
    % images as their detection is probably spurious in the first place -- incomplete PSFs (i.e.
    % overlapping the image sensor edge) are unreliably detected this close to the edge of the camera
    
    % fix x boundary events
    badidx = find(datatable.data(:,ProcSettings.xCoordsColumn)==0);
    datatable.data(badidx,ProcSettings.xCoordsColumn) = 0.1;
    badidx = find(datatable.data(:,ProcSettings.xCoordsColumn)==ProcSettings.ImageSize);
    datatable.data(badidx,ProcSettings.xCoordsColumn) = ProcSettings.ImageSize - 0.1;

    % fix y boundary events
    badidx = find(datatable.data(:,ProcSettings.yCoordsColumn)==0);
    datatable.data(badidx,ProcSettings.yCoordsColumn) = 0.1;
    badidx = find(datatable.data(:,ProcSettings.yCoordsColumn)==ProcSettings.ImageSize);
    datatable.data(badidx,ProcSettings.yCoordsColumn) = ProcSettings.ImageSize - 0.1;
    
    ROI_stats = cell(19,2);
    ROI_stats(1,:) = {'FileID',FileID};
    
%==========================================================================
%	Save an image of the events along with the ROI
%==========================================================================

    
   
    % Make a set of random indices for plotting (renders much faster for an unknown matlabby reason...
    data_rand_ind = randperm(size(datatable.data,1));                           % Create a vector of random integers from 1 to len
    data_randx = datatable.data(data_rand_ind,ProcSettings.xCoordsColumn);  % Randomize the x input vectors.
    data_randy = datatable.data(data_rand_ind,ProcSettings.yCoordsColumn);  % Randomize the y input vectors.


    SaveFileName = [FileList{FileID,2},' - Events'];
    if OverwriteExisting || ~exist(fullfile(pwd,'Events',[SaveFileName,'.png']),'file')
        fprintf('%s','Saving image of all events');
        
        % Plot all the events and save an image
        [DataEvents_fig, DataEvents_axes] = DoMeAFigure(ProcSettings.AxisLimits,[0 0 0]);
        if ProcSettings.ChannelIDColumn ~= 0
            data_randz=datatable.data(data_rand_ind,ProcSettings.ChannelIDColumn);
            EventsPlot = scatter3(data_randx,data_randy,datatable.data(:,ProcSettings.ChannelIDColumn),1,'w.',data_randz);
            clear data_randz
        else
            EventsPlot = scatter(data_randx,data_randy,1,'w.');
        end
        axis(ProcSettings.AxisLimits); % reassert image axes in case RoI strayed outside
        set(DataEvents_fig,'InvertHardcopy','off');
        print(DataEvents_fig,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'Events',[SaveFileName,'.png']));

        close(DataEvents_fig);
        clear DataEvents_fig DataEvents_axes EventsPlot
        
        fprintf('\t%s','(done)');
        fprintf('\n');
    end
        clear SaveFileName


    %plot the ROI and save another image

    SaveFileName = [FileList{FileID,2},' - ROI Preview'];
    if OverwriteExisting || ~exist(fullfile(pwd,'Events',[SaveFileName,'.png']),'file')

        fprintf('%s','Saving image of all events with ROI');
        
        [DataROIEvents_fig, DataROIEvents_axes] = DoMeAFigure(ProcSettings.AxisLimits,[0 0 0]);

        RoIClosed = pos;
        RoIClosed(end+1,:) = pos(1,:); % repeat first coord to close the ROI shape
        plot(RoIClosed(:,1),RoIClosed(:,2),'c','LineWidth',2);

        if ProcSettings.ChannelIDColumn ~= 0
            data_randz=datatable.data(data_rand_ind,ProcSettings.ChannelIDColumn);
            EventsROIPlot = scatter3(data_randx,data_randy,datatable.data(:,ProcSettings.ChannelIDColumn),1,'w.',data_randz);
            clear data_randz
        else
            EventsROIPlot = scatter(data_randx,data_randy,1,'w.');
        end

        set(EventsROIPlot,'CData',[1 0.02 0.2]);
        set(DataROIEvents_fig,'InvertHardcopy','off');

        axis(ProcSettings.AxisLimits); % reassert image axes in case RoI strayed outside

        print(DataROIEvents_fig,'-dpng',ProcSettings.SaveLowDPI,fullfile(pwd,'Events',[SaveFileName,'.png']));
        close(DataROIEvents_fig);
        
        fprintf('\t%s','(done)');
        fprintf('\n');
    end
        
    % Find the events within the region [InROI_idx, OnROI_idx]
    [InROI_idx,~] = inpolygon(datatable.data(:,ProcSettings.xCoordsColumn),datatable.data(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));
    ROIarea = polyarea(pos(:,1),pos(:,2));
    
    ROI_stats(2,:) = {'Area';ROIarea};
    ROI_stats(3,:) = {'Events Inside RoI';numel(InROI_idx)};
    
    data_proc = zeros(sum(InROI_idx),4);
    data_proc(:,1) = datatable.data(InROI_idx,ProcSettings.xCoordsColumn);
    data_proc(:,2) = datatable.data(InROI_idx,ProcSettings.yCoordsColumn);
    ProcSettings.xCol_procd = 1;
    ProcSettings.yCol_procd = 2;
    
    % TheUnwanted are the events outside of the ROI. We keep a record of these for later plots.
    TheUnwanted = horzcat(datatable.data(~InROI_idx,ProcSettings.xCoordsColumn),datatable.data(~InROI_idx,ProcSettings.yCoordsColumn));
    TheUnwanted_rand_ind = randperm(size(TheUnwanted,1));                         % Create a vector of random integers from 1 to len
    TheUnwanted_randx = TheUnwanted(TheUnwanted_rand_ind,ProcSettings.xCol_procd);  % Randomize the x input vectors.
    TheUnwanted_randy = TheUnwanted(TheUnwanted_rand_ind,ProcSettings.yCol_procd);  % Randomize the y input vectors.
    ROI_stats(4,:) = {'Events Outside RoI';size(TheUnwanted,1)};
    clear TheUnwanted_rand_ind TheUnwanted
    
    clear RoIClosed EventsROIPlot DataROIEvents_fig DataROIEvents_axes SaveFileName tmpax OnROI_idx xl yl xc yc xsc ysc xdir ydir hFH fnsave fcn EventsPlot

   
    
%             %==========================================================================
%             %	 Calculate the distances between events with local randomisation
%             %==========================================================================
%             % The array 'data_proc' holds only events within the ROI but distances are 
%             % processed relative to the entire image. It is faster to only measure for 
%             % events that are in the ROI.
% 
%                 ProcSettings.LocalRandRadius = 300;
% 
%                 ProcTimeStopwatch = tic;
%                 QueryXYtemp = horzcat(datatable.data(:,ProcSettings.xCoordsColumn),datatable.data(:,ProcSettings.yCoordsColumn));
%                 [ppFFD_inROI_tmp,ppAFD_inROI_tmp,DistancesToN_inROI,ppFFD_localrand,ppAFD_localrand,DistancesToN_localrand] = DTF2ParaFunc3(data_proc,ProcSettings,QueryXYtemp);
%                 DMPRocTime = toc(ProcTimeStopwatch);
%                 ROI_stats(18,:) = {'Distance Matrix with Local Randomisation Processing Time';DMPRocTime};
%                 clear DMPRocTime ProcTimeStopwatch
% 
%                 data_proc(:,3) = ppFFD_inROI_tmp;
%                 data_proc(:,4) = ppAFD_inROI_tmp;
% 
%                 data_proc_rndm(:,ProcSettings.FFDCol) = ppFFD_localrand;
%                 data_proc_rndm(:,ProcSettings.AFDCol) = ppAFD_localrand;
%                 % DistancesToN_localrand
% 
%                 ProcSettings.FFDCol = 3;	% FFD = furthest friend distance
%                 ProcSettings.AFDCol = 4;	% AFD = all friend distance (sum of distances to n friends)
% 
%             %     clear ppAFD_tmp ppFFD_tmp ppFFD_inROI_tmp ppAFD_inROI_tmp data_proc_inROI QueryXYtemp
%             % 
%             %     ProcTblHeaders = {'x','y',['FFD(',num2str(ProcSettings.FurthestFriendID),')'],['AFD(',num2str(ProcSettings.FurthestFriendID),')']};
%             %     horzcat(ProcTblHeaders,{['Norm''d FFD(',num2str(ProcSettings.FurthestFriendID),')']});
%             %     
%             %     % Make a set of random indices for plotting (renders much faster for an unknown matlabby reason...)
%             %     data_proc_rand_ind = randperm(size(data_proc,1));                         % Create a vector of random integers from 1 to len
%             %     data_proc_randx = data_proc(data_proc_rand_ind,ProcSettings.xCol_procd);  % Randomize the x input vectors.
%             %     data_proc_randy = data_proc(data_proc_rand_ind,ProcSettings.yCol_procd);  % Randomize the y input vectors.

    



%==========================================================================
%	 Calculate the distances between events
%==========================================================================
% The array 'data_proc' holds only events within the ROI but distnaces are 
% processed relative to the entire image. It is faster to only measure for 
% events that are in the ROI.

    ProcTimeStopwatch = tic;
    QueryXYtemp = horzcat(datatable.data(:,ProcSettings.xCoordsColumn),datatable.data(:,ProcSettings.yCoordsColumn));
    [ppFFD_inROI_tmp,ppAFD_inROI_tmp,DistancesToN_inROI] = DTF2ParaFunc2(data_proc,ProcSettings,QueryXYtemp);
    DMPRocTime = toc(ProcTimeStopwatch);
    ROI_stats(18,:) = {'Distance Matrix Processing Time';DMPRocTime};
    clear DMPRocTime ProcTimeStopwatch

    data_proc(:,3) = ppFFD_inROI_tmp;
    data_proc(:,4) = ppAFD_inROI_tmp;
    ProcSettings.FFDCol = 3;	% FFD = furthest friend distance
    ProcSettings.AFDCol = 4;	% AFD = all friend distance (sum of distances to n friends)

    clear ppAFD_tmp ppFFD_tmp ppFFD_inROI_tmp ppAFD_inROI_tmp data_proc_inROI QueryXYtemp

    ProcTblHeaders = {'x','y',['FFD(',num2str(ProcSettings.FurthestFriendID),')'],['AFD(',num2str(ProcSettings.FurthestFriendID),')']};
    horzcat(ProcTblHeaders,{['Norm''d FFD(',num2str(ProcSettings.FurthestFriendID),')']});
    
    % Make a set of random indices for plotting (renders much faster for an unknown matlabby reason...)
    data_proc_rand_ind = randperm(size(data_proc,1));                         % Create a vector of random integers from 1 to len
    data_proc_randx = data_proc(data_proc_rand_ind,ProcSettings.xCol_procd);  % Randomize the x input vectors.
    data_proc_randy = data_proc(data_proc_rand_ind,ProcSettings.yCol_procd);  % Randomize the y input vectors.

%==========================================================================
%	 Randomise the position of the events within the ROI.
%==========================================================================
% By generating x and y values to the nearest 10th nm to match the in-ROI 
% density but spread throughout the entire image area. We will then adjust
% the total number of random events until we have the same number in the 
% ROI as the real data.

    fprintf('%s','Randomising position of events within ROI');
    rndm_events_reqd = ceil(((ProcSettings.PlotWidth*ProcSettings.PlotWidth)/ROIarea) * sum(InROI_idx)); % Randomise over entire image area, matching density within ROI
    datatable.data_rndm = zeros(rndm_events_reqd,size(datatable.data,2));        
    datatable.data_rndm(:,ProcSettings.xCoordsColumn) = randi([(ProcSettings.AxisLimits(1)*10)+1, (ProcSettings.AxisLimits(2)*10)+1],rndm_events_reqd,1);
    datatable.data_rndm(:,ProcSettings.yCoordsColumn) = randi([(ProcSettings.AxisLimits(3)*10)+1, (ProcSettings.AxisLimits(4)*10)+1],rndm_events_reqd,1);
    datatable.data_rndm = datatable.data_rndm / 10; % convert values back to nm.

    % find randmised events inside of ROI and calculate how many events we are
    % short or how many are in excess. Add/subtract and repeat until we have an
    % equal number of ranomised events within the ROI, i.e. EventOffset is zero.
    [RndInROI_idx,~] = inpolygon(datatable.data_rndm(:,ProcSettings.xCoordsColumn),datatable.data_rndm(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));
    EventsOffset = sum(RndInROI_idx) - sum(InROI_idx);

    while EventsOffset ~= 0
        if  EventsOffset  < 0  %insufficient events, add more
            
            EventsToAdd = abs(EventsOffset) * ceil(((ProcSettings.PlotWidth*ProcSettings.PlotWidth)/ROIarea)); % Add events proportional to the ROI area within the image area
            tmp_new_events = zeros(EventsToAdd,size(datatable.data,2));
            tmp_new_events(:,ProcSettings.xCoordsColumn) = randi([ProcSettings.AxisLimits(1)*10, ProcSettings.AxisLimits(2)*10],EventsToAdd,1);
            tmp_new_events(:,ProcSettings.yCoordsColumn) = randi([ProcSettings.AxisLimits(3)*10, ProcSettings.AxisLimits(4)*10],EventsToAdd,1);
            tmp_new_events = tmp_new_events / 10; % convert vales to nm.
            datatable.data_rndm = vertcat(datatable.data_rndm,tmp_new_events);

            [RndInROI_idx,~] = inpolygon(datatable.data_rndm(:,ProcSettings.xCoordsColumn),datatable.data_rndm(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));
            clear EventsToAdd tmp_new_events
            EventsOffset = sum(RndInROI_idx) - sum(InROI_idx);
            
        elseif EventsOffset > 0 % excess events, delete the difference from the list of RndInROI events
            
            k_idx = find(RndInROI_idx); % indices of datatable.data_rndm that are In ROI.
            EventsToGo = sort(randi([1 size(k_idx,1)],abs(EventsOffset),1)); % pick random indices to remove
            datatable.data_rndm(k_idx(EventsToGo),:) = [];
            
            [RndInROI_idx,~] = inpolygon(datatable.data_rndm(:,ProcSettings.xCoordsColumn),datatable.data_rndm(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));
            clear k_idx EventsToGo
            EventsOffset = sum(RndInROI_idx) - sum(InROI_idx);
            
        end
    end
    data_proc_rndm = zeros(sum(RndInROI_idx),4);
    data_proc_rndm(:,ProcSettings.xCol_procd) = datatable.data_rndm(RndInROI_idx,ProcSettings.xCoordsColumn);
    data_proc_rndm(:,ProcSettings.yCol_procd) = datatable.data_rndm(RndInROI_idx,ProcSettings.yCoordsColumn);

    data_proc_rndm_randx = data_proc_rndm(data_proc_rand_ind,ProcSettings.xCol_procd);  % Randomize the x input vectors.
    data_proc_rndm_randy = data_proc_rndm(data_proc_rand_ind,ProcSettings.yCol_procd);  % Randomize the y input vectors.

    clear EventsOffset rndm_events_reqd
    
% Save an image of the ROI with randomly placed events inside

    SaveFileName = [FileList{FileID,2},' - Events (In ROI Randomised)'];
    if OverwriteExisting || ~exist(fullfile(pwd,'Events',[SaveFileName,'.png']),'file')

        [DataEvents_rndroi_fig, DataEvents_rndroi_axes] = DoMeAFigure(ProcSettings.AxisLimits,[0 0 0]);
        scatter(data_proc_rndm_randx,data_proc_rndm_randy,1,[1 1 1],'.');
        axis(ProcSettings.AxisLimits);
        set(DataEvents_rndroi_fig,'InvertHardCopy','off');

        print(DataEvents_rndroi_fig,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'Events',[SaveFileName,'.png']));
        close(DataEvents_rndroi_fig)
        
        clear SaveFileName DataEvents_rndroi_fig DataEvents_rndroi_axes
    
    end
    
    fprintf('\t%s','(done)');
    fprintf('\n');
    clear SaveFileName

        
% measure distances for the randomised data (within ROI) relative to the 
% entire randomised data set. Array data_proc_rndm holds only in-ROI
% events, as for the regular data.
    QueryXYtemp = horzcat(datatable.data_rndm(:,ProcSettings.xCoordsColumn),datatable.data_rndm(:,ProcSettings.yCoordsColumn));
    [ppFFD_rndm_inROI_tmp,ppAFD_rndm_inROI_tmp,DistancesToN_rndm] = DTF2ParaFunc2(data_proc_rndm,ProcSettings,QueryXYtemp);    
    data_proc_rndm(:,ProcSettings.FFDCol) = ppFFD_rndm_inROI_tmp;
    data_proc_rndm(:,ProcSettings.AFDCol) = ppAFD_rndm_inROI_tmp;

    clear QueryXYtemp ppAFD_tmp ppFFD_tmp ppFFD_rndm_inROI_tmp ppAFD_rndm_inROI_tmp ppFFD_rnd_tmp ppAFD_rnd_tmp DTF2vars
    clear pos % ROI position information is no longer required (and is saved to MAT file already)

%==========================================================================
%	 Histograms: Single distances to Nth NN, overlaid
%==========================================================================

    fprintf('%s ', 'Saving images:');
    fprintf('\n');
  
% Single distances for data (in ROI)
    fprintf('\t%s','NN Distance histograms (events within ROI)'); 
    
    disthist_max_inROI = ceil(max(max(DistancesToN_inROI))/10)*10;
    disthist_in_roi=[];
    for dn = 1:ProcSettings.FurthestFriendID
         disthist_in_roi(dn,:) = hist(DistancesToN_inROI(dn,:),0:0.1:disthist_max_inROI);
    end
%     disthist_in_roi = disthist_in_roi ./ sum(InROI_idx) * 100;
    disthist_labels = 0:0.1:disthist_max_inROI;
    
    xaxismax = 10*ceil(0.2*median(DistancesToN_inROI(ProcSettings.FurthestFriendID,:)));
    yaxismax = ceil(1.1*max(max(disthist_in_roi)));

    fade_array = linspace(1,0,ProcSettings.FurthestFriendID)';
    fade_r2b = horzcat(fade_array,zeros(ProcSettings.FurthestFriendID,1),flipud(fade_array));

    disthistoverlay_inROI = figure('Visible',FigureVisibility);
    p01 = plot(disthist_labels,disthist_in_roi(1,:),'-','Color',fade_r2b(1,:));
    hold on
    if ProcSettings.FurthestFriendID > 1
        p02 = plot(disthist_labels,disthist_in_roi(2,:),':','Color',fade_r2b(2,:));
        for p=3:ProcSettings.FurthestFriendID-1
            plot(disthist_labels,disthist_in_roi(p,:),':','Color',fade_r2b(p,:));
        end
        p_end = plot(disthist_labels,disthist_in_roi(ProcSettings.FurthestFriendID,:),'-','Color',fade_r2b(ProcSettings.FurthestFriendID,:));
        legend([p01 p02 p_end],{'n=1',['n=2-',num2str(ProcSettings.FurthestFriendID-1)],['n=',num2str(ProcSettings.FurthestFriendID),'']});
    end
    title([ProcSettings.ExptTitle,' : Distances to n^{th} Neighbour (within ROI)']);
    xlabel('Distance to n^{th} nearest neighbour (nm)');
    ylabel('Events');
    %set(gca,'YLim',[0 ceil(max(max(disthist_in_roi(:,2:end-1)))/100)*100]);
    set(gca,'YLim',[0 yaxismax]);
    set(gca,'XLim',[0 xaxismax]); % should calc the upper limit automagically
    set(gca,'TickDir','out');
    xlabel('Distance to n^{th} neighbour (nm)');
    ylabel('Events');

    SaveFileName = [FileList{FileID,2},' - Distances to nth Neighbour (within ROI)'];
    print(disthistoverlay_inROI,'-dpng','-r300',fullfile(pwd,'Histograms','Distances to Nth NN',SaveFileName));

    fprintf('\t%s','(done)');
    fprintf('\n');
        
% Single distances for randomised (in ROI)
    fprintf('\t%s','NN Distance histograms (randomised events within ROI)'); 

    disthist_rnd_max = ceil(max(max(DistancesToN_rndm))/10)*10;
    disthist_rnd=[];
    for dn = 1:ProcSettings.FurthestFriendID
         disthist_rnd(dn,:) = hist(DistancesToN_rndm(dn,:),0:0.1:disthist_rnd_max);
    end
%     disthist_rnd = disthist_rnd ./ sum(InROI_idx) * 100;
    disthist_rnd_labels = 0:0.1:disthist_rnd_max;

    xaxismax = 10*ceil(0.2*median(DistancesToN_rndm(ProcSettings.FurthestFriendID,:)));
    yaxismax = ceil(1.1*max(max(disthist_rnd)));

    fade_array = linspace(1,0,ProcSettings.FurthestFriendID)';
    fade_r2b = horzcat(fade_array,zeros(ProcSettings.FurthestFriendID,1),flipud(fade_array));

    disthistoverlay_rand_roi = figure('Visible',FigureVisibility);
    p01 = plot(disthist_rnd_labels,disthist_rnd(1,:),'-','Color',fade_r2b(1,:));
    hold on
    if ProcSettings.FurthestFriendID > 1
        p02 = plot(disthist_rnd_labels,disthist_rnd(2,:),':','Color',fade_r2b(2,:));
        for p=3:ProcSettings.FurthestFriendID-1
            plot(disthist_rnd_labels,disthist_rnd(p,:),':','Color',fade_r2b(p,:));
        end
        p_end = plot(disthist_rnd_labels,disthist_rnd(ProcSettings.FurthestFriendID,:),'-','Color',fade_r2b(ProcSettings.FurthestFriendID,:));
        legend([p01 p02 p_end],{'n=1',['n=2-',num2str(ProcSettings.FurthestFriendID-1)],['n=',num2str(ProcSettings.FurthestFriendID),'']});
    end
    title([ProcSettings.ExptTitle,' : Distances to n^{th} Neighbour (randomised to ROI density)']);
    %set(gca,'YLim',[0 ceil(max(max(disthist_rnd(:,2:end-1)))/100)*100]);
    set(gca,'YLim',[0 yaxismax]);
    set(gca,'XLim',[0 xaxismax]); % should calc the upper limit automagically
    set(gca,'TickDir','out');
    xlabel('Distance to n^{th} neighbour (nm)');
    ylabel('Events');

    SaveFileName = [FileList{FileID,2},' - Distances to nth Neighbour (randomised to ROI density)'];
    print(disthistoverlay_rand_roi,'-dpng','-r300',fullfile(pwd,'Histograms','Distances to Nth NN',SaveFileName));
    
    fprintf('\t%s','(done)');
    fprintf('\n');

% Clean up    
    clear dn p p01 p02 p_end  fade_array fade_r2b disthist_labels xaxismax yaxismax
        
    close(disthistoverlay_inROI);
    clear disthistoverlay_inROI disthist_max_inROI disthist_in_roi

    close(disthistoverlay_rand_roi);
    clear disthist_rnd disthist_rnd_labels disthist_rnd_max  disthistoverlay_rand_roi

%==========================================================================
%	 Log Sum Distance to NN : Images and histogram subplots
%==========================================================================

    % Calculate the log_e for the sum-distance values
        ProcSettings.logAFDCol = size(data_proc,2) + 1;
        data_proc(:,ProcSettings.logAFDCol) = log(data_proc(:,ProcSettings.AFDCol));
        data_proc_rndm(:,ProcSettings.logAFDCol) = log(data_proc_rndm(:,ProcSettings.AFDCol));
        ProcTblHeaders = horzcat(ProcTblHeaders,{['Log AFD(',num2str(ProcSettings.FurthestFriendID),')']});
    % Invert values (more clustered = higher value)
        ProcSettings.invlogAFDCol = size(data_proc,2) + 1;
        ProcTblHeaders = horzcat(ProcTblHeaders,{['Norm''d Log AFD(',num2str(ProcSettings.FurthestFriendID),')']});
        
        MaxLogAFD = max(data_proc(:,ProcSettings.logAFDCol));
        MaxLogAFD_rnd = max(data_proc_rndm(:,ProcSettings.logAFDCol));
        
        data_proc(:,ProcSettings.invlogAFDCol) = MaxLogAFD_rnd - data_proc(:,ProcSettings.logAFDCol);
        data_proc_rndm(:,ProcSettings.invlogAFDCol) = MaxLogAFD_rnd - data_proc_rndm(:,ProcSettings.logAFDCol);
        
    % Histograms - inverted values
        % for the data
        plotcollection = figure('Visible',FigureVisibility);
        invlogsumdistanceplot = [];
        xbins = floor(min(data_proc(:,ProcSettings.invlogAFDCol))):0.01:ceil(max(data_proc(:,ProcSettings.invlogAFDCol))); % bin width = 10 (nm)
        [invlogsumdistanceplot(2,:),invlogsumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.invlogAFDCol),xbins);
%         subplot(1,2,2);
        bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:),1.0,'FaceColor',ColorPlots_Data,'EdgeColor',ColorPlots_Data,'LineWidth',0.1);
        alpha(0.5);
        
        axis tight
        set(gca,'TickDir','out');
        xlabel(['Inv. Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
        ylabel('Events');
        title(['Max Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})-Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
        invlogsumdistance_rnd = max(invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1)))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.invlogAFDCol))) ceil(max(data_proc(:,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == invlogsumdistance_rnd)]);
        histogram_axis = axis;
        hold on
        
        % for the randomised data
        invlogsumdistanceplot_rnd = [];    
        [invlogsumdistanceplot_rnd(2,:),invlogsumdistanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.invlogAFDCol),xbins);
        bar(invlogsumdistanceplot_rnd(1,:),invlogsumdistanceplot_rnd(2,:),1.0,'FaceColor',ColorPlots_Rndm,'EdgeColor','none');
        alpha(0.5);
        
        %Save the histogram
        suptitle(ProcSettings.ExptTitle); % add an overall figure title
        SaveFileName = [FileList{FileID,2},' - NN histogram collection'];
        print(plotcollection,'-dpng','-r600',fullfile(pwd,'Histograms','Collections',[SaveFileName,'.png']));
        close(plotcollection)
        clear plotcollection

% % % %
% % % %   simple normalisation
% % % %

        ProcSettings.CorrectedILSD = size(data_proc,2) + 1;
        ProcTblHeaders = horzcat(ProcTblHeaders,{['Corrected ILSD(',num2str(ProcSettings.FurthestFriendID),')']});

        PeakLogAFD_rnd = median(data_proc_rndm(:,ProcSettings.invlogAFDCol)); % take the peak median
        ROI_stats(5,:) = {'PeakLogAFD_rnd (Offset)';PeakLogAFD_rnd};

        data_proc(:,ProcSettings.CorrectedILSD) = data_proc(:,ProcSettings.invlogAFDCol) - PeakLogAFD_rnd;
        data_proc_rndm(:,ProcSettings.CorrectedILSD) = data_proc_rndm(:,ProcSettings.invlogAFDCol) - PeakLogAFD_rnd;

        RangeLogAFD_rnd = max(data_proc_rndm(:,ProcSettings.CorrectedILSD)) - min(data_proc_rndm(:,ProcSettings.CorrectedILSD));
        ClusterThresholdVal = 0.5*SimpleThresholdFactor*RangeLogAFD_rnd;% the cluster threshold
        DispersedThresholdVal = -0.5*SimpleThresholdFactor*RangeLogAFD_rnd;% the dispersal threshold
        ROI_stats(6,:) = {'Threshold Factor';SimpleThresholdFactor};
        ROI_stats(7,:) = {'ClusterThresholdVal';ClusterThresholdVal};
        ROI_stats(8,:) = {'DispersedThresholdVal';DispersedThresholdVal};
                
    % for the data
        CorrectedILSD_plot = figure;
        CorrILSDplot = []; 
        xbins = floor(min(data_proc(:,ProcSettings.CorrectedILSD))):0.01:ceil(max(data_proc(:,ProcSettings.CorrectedILSD))); % bin width = 10 (nm)
        [CorrILSDplot(2,:),CorrILSDplot(1,:)] = hist(data_proc(:,ProcSettings.CorrectedILSD),xbins);
        %subplot(2,3,6);
        bar(CorrILSDplot(1,:),CorrILSDplot(2,:),1.0,'FaceColor',ColorPlots_Data,'EdgeColor',ColorPlots_Data,'LineWidth',0.1);
        alpha(0.5)
        axis tight
        set(gca,'TickDir','out'); 
        xlabel(['Corrected Inv. Log Sum d_{',num2str(ProcSettings.FurthestFriendID),'}']);
        ylabel('Events');
        title(['Corrected ILSD{',num2str(ProcSettings.FurthestFriendID),'}): ',FileList{FileID,2}]);
        CorrILSD_valmax = max(CorrILSDplot(1,CorrILSDplot(2,:)==max(CorrILSDplot(2,2:end-1)))); % ignoring first and last bins
        axis([floor(min(data_proc(:,ProcSettings.CorrectedILSD))) ceil(max(data_proc(:,ProcSettings.CorrectedILSD))) 0 1.1*CorrILSDplot(2,CorrILSDplot(1,:) == CorrILSD_valmax)]);
        hold on
        
    % for the randomised data
        CorrILSDplot_rnd = [];    
        [CorrILSDplot_rnd(2,:),CorrILSDplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.CorrectedILSD),xbins);
        bar(CorrILSDplot_rnd(1,:),CorrILSDplot_rnd(2,:),1.0,'FaceColor',ColorPlots_Rndm,'EdgeColor','none');
        alpha(0.5)
        axesmax = axis;
        ClusterThresholdLine = line([ClusterThresholdVal ClusterThresholdVal],[0 axesmax(4)],'Color','m','LineWidth',2,'LineStyle',':');
        DispersalThresholdLine = line([DispersedThresholdVal DispersedThresholdVal],[0 axesmax(4)],'Color','c','LineWidth',2,'LineStyle',':');
        
        SaveFileName = [FileList{FileID,2},' - Corrected ILSD with thresholds'];
        print(CorrectedILSD_plot,'-dpng','-r600',fullfile(pwd,'Histograms',SaveFileName));
        close(CorrectedILSD_plot)
        clear CorrectedILSD_plot ClusterThresholdLine DispersalThresholdLine SaveFileName

    %fig by simple threshold
        [fig_corrected_test, axes_corrected_test] = DoMeAFigure(ProcSettings.AxisLimits);
        set(fig_corrected_test,'visible','off');
        
        [EventsClus_idx,~] = find(data_proc(:,ProcSettings.CorrectedILSD) >= ClusterThresholdVal);
        [EventsRndm_idx,~] = find(data_proc(:,ProcSettings.CorrectedILSD) > DispersedThresholdVal & data_proc(:,ProcSettings.CorrectedILSD) < ClusterThresholdVal);
        [EventsDisp_idx,~] = find(data_proc(:,ProcSettings.CorrectedILSD) <= DispersedThresholdVal);
             
        hold on
        scatter(TheUnwanted_randx,TheUnwanted_randy,1,[0.8 0.8 0.8],'x');
        scatter(data_proc(EventsDisp_idx,ProcSettings.xCol_procd),data_proc(EventsDisp_idx,ProcSettings.yCol_procd),1,'b.');
        scatter(data_proc(EventsRndm_idx,ProcSettings.xCol_procd),data_proc(EventsRndm_idx,ProcSettings.yCol_procd),1,'g.');       
        scatter(data_proc(EventsClus_idx,ProcSettings.xCol_procd),data_proc(EventsClus_idx,ProcSettings.yCol_procd),1,'m.');
        axis(ProcSettings.AxisLimits);
        SaveFileName = [FileList{FileID,2},' - Thold by Rand Distn'];
        print(fig_corrected_test,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'D - Threshold',SaveFileName));
        close(fig_corrected_test)
        clear fig_corrected_test axes_corrected_test SaveFileName
        
    % some quick stats
        stats_corr2rand(1,FileID) = ClusterThresholdVal;
        stats_corr2rand(2,FileID) = DispersedThresholdVal;
        stats_corr2rand(3,FileID) = size(EventsClus_idx,1);
        stats_corr2rand(4,FileID) = size(EventsRndm_idx,1);
        stats_corr2rand(5,FileID) = size(EventsDisp_idx,1);
        stats_corr2rand(6,FileID) = size(data_proc,1);
        
        ROI_stats(9,:) = {'Clustered Events';size(EventsClus_idx,1)};
        ROI_stats(10,:) = {'Random-like Events'; size(EventsRndm_idx,1)};
        ROI_stats(11,:) = {'Dispersed Events';size(EventsDisp_idx,1)};
                
        SaveToFolder = fullfile('D - Threshold');
        SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(ClusterThresholdVal),' (by Simple Thresholding)'];
        RefDataCol = ProcSettings.CorrectedILSD;
        [simplethr_pts_over, simplethr_pts_under] = thresholder(data_proc,ClusterThresholdVal,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);

    % Colour map of the inverted log sum-distances to the Nth NN (for the real data)
      
        if ProcSettings.SaveClusterColourImages
            fprintf('\t%s ', ['Normalised Inverted Log Sum Distance to NN(',num2str(ProcSettings.FurthestFriendID),')']); 

        % Data: inverted log sum-distance to NN
            [fig_invlogAFD, axes_invlogAFD] = DoMeAFigure(ProcSettings.AxisLimits);
            data_proc_randz = data_proc(data_proc_rand_ind,ProcSettings.CorrectedILSD); % randomised z data using rand idx calculated in the beginning    
            scatter(TheUnwanted_randx,TheUnwanted_randy,1,[0.8 0.8 0.8],'x');
            scatter(data_proc_randx,data_proc_randy,5,data_proc_randz,'.');
            colormap(CustomColorMap);
            corr_caxis = caxis;
            corr_caxis(1) = -1;  % Pin the colour scale to -1 minimum
            caxis(corr_caxis);
            SaveFileName = [FileList{FileID,2},' - fig_Norm2Rand_InvLogAFD (Bright Background)'];
            SaveFolder = 'C - Inv Log Sum of Distances to NN';
            print(fig_invlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
        
            % change background to dark colour
            c = colormap;
            %uncorrcmapmax = caxis;
            set(axes_invlogAFD,'Color',c(1,:));
            set(fig_invlogAFD,'Color',c(1,:),'InvertHardcopy','off')
            SaveFileName = [FileList{FileID,2},' - fig_Norm2Rand_InvLogAFD (Dark Background)'];
            print(fig_invlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
        
            SaveColourBar(c,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar'])); % save current colormap
            
            close(fig_invlogAFD); 
            clear fig_invlogAFD axes_invlogAFD data_proc_randz
        end
               
        clear data_randz_invlogAFD SaveFileName c normallogsumdistance xbins SaveFolder normallogsumdistance2 c
        fprintf('\t%s','(done)');
        fprintf('\n');

%==========================================================================
%	 Correction of ILSD values to the ROI density
%    Using ILSD peak value from randomised (density equiv.) data
%==========================================================================
        if ProcSettings.NormaliseToDensity
    
        % Histograms for the density-correction of the inverted log sum-distances to the Nth NN

        % Find the peak of the ROI-randomised data
            AFDdisthist_rnd_max = ceil(max(max(data_proc_rndm(:,ProcSettings.AFDCol)))/10)*10;
            AFDdisthist_rnd = hist(data_proc_rndm(:,ProcSettings.AFDCol),0:0.1:AFDdisthist_rnd_max);
            %AFDdisthist_rnd = (AFDdisthist_rnd ./ sum(RndInROI_idx)) * 100;
            AFDdisthist_rnd_labels = 0:0.1:AFDdisthist_rnd_max;

        % Histogram of the raw AFD values for randomised & real data
            % AFD randomised data (first, to go beneath the real)
            plot_correction_collection = figure;
            subplot(1,2,1);
            rand_plot2 = bar(AFDdisthist_rnd_labels,AFDdisthist_rnd,1.0,'FaceColor',[1 0.875 0.875],'EdgeColor','none');
            alpha(0.5)
            hold on
            
            % AFD real data
            [sumdistanceplot_data(2,:),sumdistanceplot_data(1,:)] = hist(data_proc(:,ProcSettings.AFDCol),0:0.1:AFDdisthist_rnd_max);
%             sumdistanceplot_data(2,:) = sumdistanceplot_data(2,:) ./ sum(size(datatable.data,1)) * 100;
            data_plot = bar(sumdistanceplot_data(1,:),sumdistanceplot_data(2,:),1.0,'FaceColor',ColorPlots_Data,'EdgeColor',ColorPlots_Data,'LineWidth',0.1);
            alpha(0.5)
            % set comfortable axis limits
            axis tight
            tmp_x_axis_max = max([ceil(max(AFDdisthist_rnd_labels)/50)*50 , ceil(max(sumdistanceplot_data(1,2:end-1))/50)*50]);
            tmp_y_axis_max = max([ceil(1.2*max(AFDdisthist_rnd)*50)/50 , floor(1.2*max(sumdistanceplot_data(2,2:end-1))*50)/50]);
            axis([0 tmp_x_axis_max 0 tmp_y_axis_max]);
            clear tmp_x_axis_max tmp_y_axis_max
        
            % Fit a curve using an Epanechnikov kernel dist fit
            pdistn2 = fitdist(data_proc_rndm(:,ProcSettings.AFDCol),'Kernel','Kernel','epanechnikov');
            data_rnd_normdist2 = pdf(pdistn2,AFDdisthist_rnd_labels); % Calculate values for plotting a different PDF

            % Normalize the PDF to match the real-data histogram (for overlay)
            rangex = max(AFDdisthist_rnd_labels) - min(AFDdisthist_rnd_labels);     % Finds the range of this data.
            binwidth = rangex/size(AFDdisthist_rnd_labels,2);                       % Finds the width of each bin.
            plotarea = max(AFDdisthist_rnd) * binwidth;                             % Find the data plot area
            pdfarea2 = max(data_rnd_normdist2) * binwidth;                          % Find the PDF plot area
            y = (plotarea / pdfarea2) * pdf(pdistn2,AFDdisthist_rnd_labels);        % Expand PDF plot to match data area
            hold on
            
            % Add normalised fitted curve to the histogram
            plot(AFDdisthist_rnd_labels,y,'-','Color',ColorPlots_Rndm,'LineWidth',2);   % plot the fitted line
            set(gca,'TickDir','out');
            xlabel(['Sum distance to NN_{',num2str(ProcSettings.FurthestFriendID),'} (nm)']);
            ylabel('Events');
            title(['Sum of distances to NN_{',num2str(ProcSettings.FurthestFriendID),'}']);

        %Find the peak of the fitted curve and convert it to an ILSD value
            [Epanechnikov_peak_val, Epanechnikov_peak_idx] = max(data_rnd_normdist2);
            Epanechnikov_peak_AFD = AFDdisthist_rnd_labels(1,Epanechnikov_peak_idx);
            RandInROI_CorrFactor = MaxLogAFD_rnd - log(Epanechnikov_peak_AFD);
            
            ROI_stats(12,:) = {'Epanechnikov CorrFactor';RandInROI_CorrFactor};

        % Offset the Inverse Log of sum of distances to Nth NN ot the peak of the randomised data
            ProcSettings.EpanechnikovInvlogAFDCol = size(data_proc,2) + 1;
            ProcTblHeaders = horzcat(ProcTblHeaders,{['Epanechnikov corrected InvLogAFD(',num2str(ProcSettings.FurthestFriendID),')']});

            data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol) = data_proc(:,ProcSettings.invlogAFDCol) - RandInROI_CorrFactor;
            data_proc_rndm(:,ProcSettings.EpanechnikovInvlogAFDCol) = data_proc_rndm(:,ProcSettings.invlogAFDCol) - RandInROI_CorrFactor;
         
        % Histogram of the corrected, inverted log-sum-distances to the Nth NN
            subplot(1,2,2);
            % the corrected data
            corrinvlogsumdistanceplot = [];
            xbins = floor(min(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol))):0.01:ceil(max(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol))); % bin width = 10 (nm)
            [corrinvlogsumdistanceplot(2,:),corrinvlogsumdistanceplot(1,:)] = hist(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol),xbins);
            bar(corrinvlogsumdistanceplot(1,:),corrinvlogsumdistanceplot(2,:),1.0,'FaceColor',ColorPlots_Data,'EdgeColor',ColorPlots_Data,'LineWidth',0.1);
            alpha(0.5)
            axis tight
            set(gca,'TickDir','out');
            xlabel(['Corrected inv.Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
            ylabel('Events');
            title(['Corrected inv Log(Sum d_{',num2str(ProcSettings.FurthestFriendID),'})']);
            normalcorrinvlogsumdistance = max(corrinvlogsumdistanceplot(1,corrinvlogsumdistanceplot(2,:)==max(corrinvlogsumdistanceplot(2,2:end-1)))); % ignoring first and last bins
            axis([floor(min(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol))) ceil(max(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol))) 0 1.1*corrinvlogsumdistanceplot(2,corrinvlogsumdistanceplot(1,:) == normalcorrinvlogsumdistance)]);
            hold on
            
            % the corrected randomised data
            corrinvlogsumdistanceplot_rnd = [];    
            [corrinvlogsumdistanceplot_rnd(2,:),corrinvlogsumdistanceplot_rnd(1,:)] = hist(data_proc_rndm(:,ProcSettings.EpanechnikovInvlogAFDCol),50);
            bar(corrinvlogsumdistanceplot_rnd(1,:),corrinvlogsumdistanceplot_rnd(2,:),1.0,'FaceColor',[1 0.875 0.875],'EdgeColor','none');
            set(gca,'XLim',([floor(min(data_proc_rndm(:,ProcSettings.EpanechnikovInvlogAFDCol))) ceil(max(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol)))])); % fix axis to show the corrected histogram positioning including the randomised data

            suptitle(ProcSettings.ExptTitle); % add an overall figure title
            SaveFileName = [FileList{FileID,2},' - Epanechnikov corrected histograms in ROI'];
            print(plot_correction_collection,'-dpng','-r600',fullfile(pwd,'Histograms','Collections',[SaveFileName,'.png']));
            
            close(plot_correction_collection)
            clear plot_correction_collection xbins Epanechnikov_peak_AFD pdfarea2 y pdistn2 ...
                  data_rnd_normdist2 pdfarea2 rangex binwidth plotarea corrinvlogsumdistanceplot_rnd ...
                  AFDdisthist_rnd_labels AFDdisthist_rnd_max AFDdisthist_rnd ...
                  data_plot Epanechnikov_peak_idx Epanechnikov_peak_val rand_plot2 pH normalcorrinvlogsumdistance ...
                  sumdistanceplot_data MaxLogAFD MaxLogAFD_rnd
        
        % Lone histogram plot for Corrected InvLogSumDistNN alone
            fig_plot_corrinvlogAFD = figure('Visible','off');
            bar(corrinvlogsumdistanceplot(1,:),corrinvlogsumdistanceplot(2,:),1.0);
            axis tight
            set(gca,'TickDir','out');
            title([ProcSettings.ExptTitle,' : Epanechnikov Corrected Inv. Log of Sum of distances to NN_{1}–NN_{',num2str(ProcSettings.FurthestFriendID),'}']);
            %corrinvlogsumdistance = corrinvlogsumdistanceplot(1,corrinvlogsumdistanceplot(2,:)==max(corrinvlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
            axis(histogram_axis);
            %axis([0 ceil(max(data_proc(:,ProcSettings.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == invlogsumdistance)]); % match axis of uncorrected data
            SaveFileName = [FileList{FileID,2},' - Epanechnikov Corrected InvLogSumDist-NN histogram'];
            print(fig_plot_corrinvlogAFD,'-dpng','-r600',fullfile(pwd,'Histograms','ILSD-NN Histogram',[SaveFileName,'.png']));
            close(fig_plot_corrinvlogAFD)
            clear fig_plot_corrinvlogAFD invlogsumdistance
      
        % Colour map of the corrected inverted log sum-distances to the Nth NN
            if ProcSettings.SaveClusterColourImages
                fprintf('\t%s ', 'Inv log sum of distances to Nth NN, normalised to random by Epanechnikov'); 
                        
            % Data: corrected inverted log sum-distance to NN
                [fig_corrinvlogAFD, axes_corrinvlogAFD] = DoMeAFigure(ProcSettings.AxisLimits);
                data_proc_randz = data_proc(data_proc_rand_ind,ProcSettings.EpanechnikovInvlogAFDCol); % randomised z data using rand idx calculated in the beginning    
                scatter(TheUnwanted_randx,TheUnwanted_randy,1,[0.8 0.8 0.8],'x');
                scatter(data_proc_randx,data_proc_randy,5,data_proc_randz,'.');
                colormap(CustomColorMap);
                corr_caxis = caxis;
                corr_caxis(1) = -1;  % Pin the colour scale to -1 minimum
                caxis(corr_caxis);
                
                SaveFileName = [FileList{FileID,2},' - fig_Norm2Epanechnikov_InvlogAFD (Bright Background)'];
                SaveFolder = 'C - Inv Log Sum of Distances to NN';
                print(fig_corrinvlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));

                % change background to dark colour
                c = colormap;                           
                set(axes_corrinvlogAFD,'Color',c(1,:));
                set(fig_corrinvlogAFD,'Color',c(1,:),'InvertHardcopy','off')
                SaveFileName = [FileList{FileID,2},' - fig_Norm2Epanechnikov_InvlogAFD (Dark Background)'];
                print(fig_corrinvlogAFD,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));

                SaveColourBar(c,caxis,fullfile(pwd,SaveFolder,[SaveFileName,'_colourbar'])); % save current colormap
                
                close(fig_corrinvlogAFD);
                clear fig_corrinvlogAFD axes_corrinvlogAFD data_randz_corrinvlogAFD c data_proc_randz
                
            % Randomised: corrected inverted log sum-distance to NN
                [fig_corrinvlogAFD_rnd, axes_corrinvlogAFD_rnd] = DoMeAFigure(ProcSettings.AxisLimits);
                data_proc_rndm_randz = data_proc_rndm(data_proc_rand_ind,ProcSettings.EpanechnikovInvlogAFDCol); % randomised z data using rand idx calculated in the beginning    
                scatter(TheUnwanted_randx,TheUnwanted_randy,1,[0.8 0.8 0.8],'x');
                scatter(data_proc_rndm_randx,data_proc_rndm_randy,5,data_proc_rndm_randz,'.');
                colormap(CustomColorMap);
                
                caxis(corr_caxis); % match colour axis of randomised data to match that of the real data
                
                SaveFileName = [FileList{FileID,2},' - fig_Norm2Epanechnikov_InvlogAFD (Randomised - Bright Background)'];
                SaveFolder = 'C - Inv Log Sum of Distances to NN';
                print(fig_corrinvlogAFD_rnd,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));

                % change background to dark colour
                c = colormap;
                set(axes_corrinvlogAFD_rnd,'Color',c(1,:));
                set(fig_corrinvlogAFD_rnd,'Color',c(1,:),'InvertHardcopy','off')
                SaveFileName = [FileList{FileID,2},' - fig_Norm2Epanechnikov_InvlogAFD (Randomised - Dark Background)'];
                print(fig_corrinvlogAFD_rnd,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,SaveFolder,[SaveFileName,'.png']));
                                    
                close(fig_corrinvlogAFD_rnd);
                clear fig_corrinvlogAFD_rnd axes_corrinvlogAFD_rnd data_proc_rndm_randz c corr_caxis
                                    
            end
            
            % epanchnikov threshold
            RangeLogAFD_rnd = max(data_proc_rndm(:,ProcSettings.CorrectedILSD)) - min(data_proc_rndm(:,ProcSettings.CorrectedILSD));

            EpanechnikovClusterThresholdVal = 0.5*SimpleThresholdFactor*RangeLogAFD_rnd;% the cluster threshold
            EpanechnikovDispersedThresholdVal = -0.5*SimpleThresholdFactor*RangeLogAFD_rnd;% the dispersal threshold
            ROI_stats(13,:) = {'EpanechnikovClusterThresholdVal';EpanechnikovClusterThresholdVal};
            ROI_stats(14,:) = {'EpanechnikovDispersedThresholdVal';EpanechnikovDispersedThresholdVal};

        %fig by ekov threshold
            [fig_corrected_ekov, axes_corrected_ekov] = DoMeAFigure(ProcSettings.AxisLimits);
            set(fig_corrected_ekov,'visible','off');

            [EkovEventsClus_idx,~] = find(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol) >= EpanechnikovClusterThresholdVal);
            [EkovEventsRndm_idx,~] = find(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol) > EpanechnikovDispersedThresholdVal & data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol) < EpanechnikovClusterThresholdVal);
            [EkovEventsDisp_idx,~] = find(data_proc(:,ProcSettings.EpanechnikovInvlogAFDCol) <= EpanechnikovDispersedThresholdVal);

            hold on
            scatter(TheUnwanted_randx,TheUnwanted_randy,1,[0.8 0.8 0.8],'x');
            scatter(data_proc(EkovEventsDisp_idx,ProcSettings.xCol_procd),data_proc(EkovEventsDisp_idx,ProcSettings.yCol_procd),1,'b.');
            scatter(data_proc(EkovEventsRndm_idx,ProcSettings.xCol_procd),data_proc(EkovEventsRndm_idx,ProcSettings.yCol_procd),1,'g.');       
            scatter(data_proc(EkovEventsClus_idx,ProcSettings.xCol_procd),data_proc(EkovEventsClus_idx,ProcSettings.yCol_procd),1,'m.');
            axis(ProcSettings.AxisLimits);
            SaveFileName = [FileList{FileID,2},' - Thold by Rand Distn and Epanechnikov Fit'];
            print(fig_corrected_ekov,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,'D - Threshold',SaveFileName));
            close(fig_corrected_ekov)
            clear fig_corrected_ekov axes_corrected_ekov SaveFileName

        % some quick stats
            stats_ekovcorr2rand(1,FileID) = EpanechnikovClusterThresholdVal;
            stats_ekovcorr2rand(2,FileID) = EpanechnikovDispersedThresholdVal;
            stats_ekovcorr2rand(3,FileID) = size(EkovEventsClus_idx,1);
            stats_ekovcorr2rand(4,FileID) = size(EkovEventsRndm_idx,1);
            stats_ekovcorr2rand(5,FileID) = size(EkovEventsDisp_idx,1);
            stats_ekovcorr2rand(6,FileID) = size(data_proc,1);

            ROI_stats(15,:) = {'Epanechnikov Clustered Events';size(EkovEventsClus_idx,1)};
            ROI_stats(16,:) = {'Epanechnikov Random-like Events'; size(EkovEventsRndm_idx,1)};
            ROI_stats(17,:) = {'Epanechnikov Dispersed Events';size(EkovEventsDisp_idx,1)};

            SaveToFolder = fullfile('D - Threshold');
            SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(EpanechnikovClusterThresholdVal),' (by Epanechnikov Thresholding)'];
            RefDataCol = ProcSettings.EpanechnikovInvlogAFDCol;
            [simplethr_pts_over, simplethr_pts_under] = thresholder(data_proc,EpanechnikovClusterThresholdVal,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);

            clear data_randz_invlogAFD SaveFileName c normallogsumdistance normallogsumdistance2 xbins ...
                  MaxLogAFD SaveFolder uncorrcmapmax data_randz_corrinvlogAFD_rndm  ...
                  invlogsumdistanceplot_rnd logsumdistanceplot logsumdistanceplot_rnd logsumdistanceplot logsumdistanceplot_rnd
            
            fprintf('\t%s','(done)');
            fprintf('\n');

        end

    

%% Thresholding

        if ProcSettings.DoThresholds
			
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
                stats_thrfix(3,FileID) = size(data_proc,1);
                stats_thrfix(4,FileID) = (num_events_over/size(data_proc,1)) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');

                
        %====== Theoretical thresholding - from randdata standard curve
                TheoreticalThreshold = (0.09057*log10((size(data_proc,1)/ROIarea))) + 0.3693;

                fprintf('\t%s','Theoretical threshold > ',num2str(TheoreticalThreshold));
                
                SaveToFolder = fullfile('D - Threshold',['Theoretical (',num2str(TheoreticalThreshold),')']);
                SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(TheoreticalThreshold),' (user-specified)'];
                RefDataCol = ProcSettings.CorrectedILSD;

                [theothr_pts_over, theothr_pts_under] = thresholder(data_proc,TheoreticalThreshold,RefDataCol,SaveToFolder,SaveFileName,ProcSettings);
                num_events_over = size(theothr_pts_over,1);

            % update the summary table
                stats_thrtheo(1,FileID) = TheoreticalThreshold;
                stats_thrtheo(2,FileID) = num_events_over;
                stats_thrtheo(3,FileID) = size(data_proc,1);
                stats_thrtheo(4,FileID) = (num_events_over/size(data_proc,1)) * 100;

                fprintf('\t%s','(done)');
                fprintf('\n');
                
        %====== Save Adaptive thresholding data

                HeaderFormat = '%s\t%s\t%s';
                fn_stats_save = [FileList{FileID,2},' - Threshold Data.txt'];
                fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_stats_save),'w');

                fprintf(fid,HeaderFormat,'[Threshold Method]','[Threshold Value]','[% In Clusters]');
                fprintf(fid,'\r\n');

                DataFormat = '%s\t%.3f\t%.3f';

                fprintf(fid,DataFormat,'Fixed',Cluster_thr_fixed,stats_thrfix(4,FileID));
                fprintf(fid,'\r\n');
                
                fprintf(fid,DataFormat,'Auto By Peak Random',ClusterThresholdVal,(stats_corr2rand(3,FileID)/stats_corr2rand(6,FileID))*100);
                fprintf(fid,'\r\n');

                fprintf(fid,DataFormat,'Auto By Epanechnikov Fit',EpanechnikovClusterThresholdVal,(stats_ekovcorr2rand(3,FileID)/stats_ekovcorr2rand(6,FileID))*100);
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
                    'ProcTblHeaders'...
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
        
        ROI_stats(19,:) = {'Processing Time';ThisRegionTimestamp - PreviousRegionTimestamp};

        %====== Save summary stats data
        fn_stats_save = [FileList{FileID,2},' - ROI Summary Data.txt'];
        fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_stats_save),'w');
        DataFormat = '%s\t%.3f\r\n';
        for r = 1:size(ROI_stats,1)
            fprintf(fid,DataFormat,ROI_stats{r,:});
        end
        fclose(fid);
        clear fid

        InfoMessage =  [datestr(fix(clock),'HH:MM:SS'),9,'Completed processing Table ',num2str(FileID),'.'];
        disp(InfoMessage);
        
        InfoMessage =  [9,9,9,sprintf('%0.2f',((ThisRegionTimestamp - PreviousRegionTimestamp)/60)),' minutes, ',num2str(size(datatable.data,1)),' events, ',num2str(floor(size(datatable.data,1)/(ThisRegionTimestamp - PreviousRegionTimestamp))),' events per second.'];
        disp(InfoMessage);
        
        InfoMessage = [9,9,9,num2str(round(FileID / size(FileList,1) * 100)) '% processed.'];
        disp(InfoMessage);
        
        disp('---------------------------------------------------------');
        
        % Update the tracking info to reflect the newly finished region
        PreviousRegionTimestamp = ThisRegionTimestamp;

end % end of this data table processing
    
%% Finish up - Overall Summaries

%====== Save Adaptive thresholding summary data

if ProcSettings.DoThresholds
    % format for headers (first row)
    stringy = '%s';
    for g = 1:(size(FileList,1)-1)
        stringy = strcat(stringy,'\t%s'); % replace \t with a comma for csv
    end

    % user-specified
    fn_summary_proc_save = ['Combined Thresholds - User Specified (',num2str(ProcSettings.FixedThreshold),').txt'];
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - User Specified (',num2str(ProcSettings.FixedThreshold),')');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_thrfix,'-append','delimiter','\t','precision','%.3f');
    
   
    % Corrected to Randomised (by peak)
    fn_summary_proc_save = 'Combined Thresholds - Corrected to Randomised by median peak.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Corrected to Randomised (Simple Peak)');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_corr2rand,'-append','delimiter','\t','precision','%.3f');

    
    % Corrected to Randomised (by fitted peak)
    fn_summary_proc_save = 'Combined Thresholds - Corrected to Randomised by fitted peak.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Corrected to Randomised (Epanechnikov)');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_ekovcorr2rand,'-append','delimiter','\t','precision','%.3f');

    % Corrected to CSR_theoretical
    fn_summary_proc_save = 'Combined Thresholds - Theoretical thresold.txt';
    fid = fopen(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),'w');
    fprintf(fid,'%s','Combined Thresholds - Theoretical');
    fprintf(fid,'\r\n');
    fprintf(fid,stringy,FileList{:,1});
    fprintf(fid,'\r\n');
    fclose(fid);
    dlmwrite(fullfile(pwd,'Y - Threshold Stats',fn_summary_proc_save),stats_thrtheo,'-append','delimiter','\t','precision','%.3f');

    
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
    
% eof
