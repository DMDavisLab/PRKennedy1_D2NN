%
%           Use this as 'Phase 1' to load data files, verify settings,
%           and selection a Region of Interest (ROI) containing your
%           cell or whatever you are looking at ;-)
%
ThisVersion = 0.5;
DoubleCheckAllSettings = true;
FigureVisibility = 'on';


%% Begin

% housekeeping
home %clean up the command window
rng('shuffle')  % set the random seed
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary'); % disables a warning about a temporary variable in the parallel processing section. This is fine but warnings scare people so I am turning it off.
warning('off','MATLAB:MKDIR:DirectoryExists'); % don't warn about existing data folders

InfoMessage = ['---------------------Phase 1 - ROI Selection (v',num2str(ThisVersion),')---------------------'];
disp(InfoMessage);
clear InfoMessage ThisVersion

%% Get list of files

% Find the coords file
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
        
        % get list of txt files
        DirFullList = dir;                              % Get the data for the current directory
        FolderIndex = [DirFullList.isdir];              % Find the index for directories
        FileList = {DirFullList(~FolderIndex).name}';   % Get a list of the files
        
        % Remove unwanted files from the list
        badFiles = [];                              
        for f = 1:length(FileList)
            [~, fname, ext] = fileparts(FileList{f,1});
            FileList{f,2} = fname;
            FileList{f,3} = ext;
            
            if ~strcmp(ext,['.',ProcSettings.DataTableFileExt])
                badFiles(end+1,1) = f;
            end
            if ~isempty(strfind(fname,'ProcSettings')) || ~isempty(strfind(fname,'coords')) || ~isempty(strfind(fname,'notes'))
                badFiles(end+1,1) = f;
            end
        end
        FileList(badFiles,:) = [];
        
        % Sort what remains
        [~, natsortidx] = sort_nat(FileList(:,2),'Ascend');
        FileList = FileList(natsortidx,:);
        
    else
        
        error('Cancelled?! So rude.');
        
    end
       
    clear badFiles DirFullList ext f fname FolderIndex natsortidx
 
% Display some information
    disp('---------------------------------------------------------');
    InfoMessage=[datestr(fix(clock),'HH:MM:SS'),9,'Preparing ' num2str(size(FileList,1)) ' data tables.'];
    disp(InfoMessage);
    disp('---------------------------------------------------------');
    clear InfoMessage

%% Set up a results folder
    
% Check for an existing data folder
if exist(fullfile(DataDirName,'SaveDataFolderName.mat'), 'file') ~= 0
%     load(fullfile(DataDirName,'SaveDataFolderName.mat'));
    errordlg('Phase 1 folder already exists! You need to move to Phase 2.');
else
    alphanums = ['a':'z' 'A':'Z' '0':'9'];
    randname = alphanums(randi(numel(alphanums),[1 5]));
    OutputFolderName = ['D2F2-(',randname,')'];
    mkdir(fullfile(DataDirName,OutputFolderName));
    cd(fullfile(DataDirName,OutputFolderName));
    mkdir(fullfile(pwd,'Z - MAT Files'));
end
   
% correct the 8 bit RGB colours to fractional
 if sum(ProcSettings.ColourEventsOver)+sum(ProcSettings.ColourEventsUnder)+sum(ProcSettings.ColourBackground) > 1
    ProcSettings.ColourEventsOver = ProcSettings.ColourEventsOver ./ 255;
    ProcSettings.ColourEventsUnder = ProcSettings.ColourEventsUnder ./ 255;
    ProcSettings.ColourBackground = ProcSettings.ColourBackground ./ 255;
 end
 
ProcSettings.SaveLowDPI = '-r600'; % for quicker saving of preview images

% clean up
clear alphanums randname

    for FileID = 1:size(FileList,1)
    
    % Load the file
        CurrentTableFileName=FileList{FileID,1};
        InfoMessage =  [datestr(fix(clock),'HH:MM:SS'),9,'[',num2str(FileID),'/',num2str(size(FileList,1)),'] ',CurrentTableFileName];
        disp(InfoMessage);
        clear InfoMessage
        
        % error if the data table file doesn't exist
        if exist(fullfile(DataDirName,CurrentTableFileName), 'file') == 0
            ErrorMessage = ['Problem! Cannot find the data table file ', CurrentTableFileName, '.'];
            errordlg(ErrorMessage);
            clear ErrorMessage
        end

        datatable = importdata(fullfile(DataDirName,CurrentTableFileName));

        % Check that we have a useable struct array with the data
        % if not make one from the data we did import.
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
        
        % correct for any scaling, i.e. convert distances to real units (nm)
        if ProcSettings.DataTableScale ~= 1
            datatable.data(:,ProcSettings.xCoordsColumn) = datatable.data(:,ProcSettings.xCoordsColumn) * ProcSettings.DataTableScale;
            datatable.data(:,ProcSettings.yCoordsColumn) = datatable.data(:,ProcSettings.yCoordsColumn) * ProcSettings.DataTableScale;
        end
        
        % Determine the image boundaries from ProcSettings or try to guess
        % them from the imported table
        if isfield(ProcSettings,'ImageSize');
            ProcSettings.AxisLimits = [0 ProcSettings.ImageSize 0 ProcSettings.ImageSize];
        else
            % try to guess the limits from the image
            minX = floor(min(datatable.data(:,ProcSettings.xCoordsColumn))/1000)*1000;
            minY = floor(min(datatable.data(:,ProcSettings.yCoordsColumn))/1000)*1000;
            maxX = ceil(max(datatable.data(:,ProcSettings.xCoordsColumn))/1000)*1000;
            maxY = ceil(max(datatable.data(:,ProcSettings.yCoordsColumn))/1000)*1000;
            
            % assume images are square ...
            minXY = min(minX,minY);
            maxXY = max(maxX,maxY);
            ProcSettings.AxisLimits = [minXY maxXY minXY maxXY];

            clear minX minY maxX maxY minXY maxXY
        end
        
        ProcSettings.PlotWidth = ProcSettings.AxisLimits(2) - ProcSettings.AxisLimits(1);
        ProcSettings.SaveHighDPI = strcat('-r',num2str(ProcSettings.PlotWidth / 10)); % for 10 nm per pixel
        
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
        
        % set up some other stuff that we'll need for the current image table...
        ProcSettings.ExptTitle = FileList{FileID,2};
        
        % Adaptive region cropping
        data_density = size(datatable.data,1) / (ProcSettings.AxisLimits(2) * ProcSettings.AxisLimits(4));
        ProcSettings.TestRegionSize_for_minFF = ceil(sqrt(ProcSettings.FurthestFriendID / data_density));
        clear data_density
        
        

    %====== Plot the events
        % randomising the input rows makes matlab render faster for some
        % unknown reason. The rand idx here are kept for later plots.
        data_rand_ind = randperm(size(datatable.data,1));                           % Create a vector of random integers from 1 to len
        data_randx = datatable.data(data_rand_ind,ProcSettings.xCoordsColumn);  % Randomize the x input vectors.
        data_randy = datatable.data(data_rand_ind,ProcSettings.yCoordsColumn);  % Randomize the y input vectors.

        [DataEvents_fig, DataEvents_axes] = DoMeAFigure(ProcSettings.AxisLimits,[0 0 0]);
              
        if ProcSettings.ChannelIDColumn ~= 0
            data_randz=datatable.data(data_rand_ind,ProcSettings.ChannelIDColumn);
            EventsPlot = scatter3(data_randx,data_randy,datatable.data(:,ProcSettings.ChannelIDColumn),1,'w.',data_randz);
            set(DataEvents_fig,'InvertHardcopy','off');
        else
            EventsPlot = scatter(data_randx,data_randy,1,'w.');
        end

        

    %====== Ask user for a ROI

        
        fprintf('%s\n','Please select a Region of Interest...');

        set(DataEvents_fig,'visible','on');
        set(DataEvents_fig,'units','normalized','position',[0 0 1 1]);
        set(EventsPlot,'CData',[1 0.02 0.2]);
        set(DataEvents_fig,'Color','k');

        % make a temporary new set of axes to trick matlab in being faster while
        % drawing
        xl = get(gca,'xlim'); 
        yl = get(gca,'ylim'); 
        posn = get(gca,'pos'); 
        xc = get(gca,'xcolor'); 
        yc = get(gca,'ycolor'); 
        xsc = get(gca,'xscale'); 
        ysc = get(gca,'yscale'); 
        xdir = get(gca,'xdir'); 
        ydir = get(gca,'ydir');


        ROI_unsaved = true;

        while ROI_unsaved
            tmpax = axes('pos',posn,...
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

            % Construct a questdlg with three options
            CheckHappiness = questdlg('Are you happy with your ROI?', ...
                'Happiness Check', ...
                'Yes (save and continue)','No (reselect ROI)','Yes (save and continue)');
            % Handle response
            switch CheckHappiness
                case 'Yes (save and continue)'
                    disp('ROI is satisfacrory. Saving ROI and moving to next image...')
                    % Delete temporary set of axes: 
                    delete(tmpax)
                    close(DataEvents_fig)
                    clear DataEvents_fig DataEvents_axes tmpax xl yl xc yc xsc ysc xdir ydir hFH fnsave fcn EventsPlot

                    % save these image-specific settings for later (so that we don't have to overwrite ProcSettings)
                    ProcSettingsLocal.AxisLimits = ProcSettings.AxisLimits;
                    ProcSettingsLocal.PlotWidth = ProcSettings.PlotWidth;
                    ProcSettingsLocal.SaveHighDPI = ProcSettings.SaveHighDPI;
                    ProcSettingsLocal.ExptTitle = ProcSettings.ExptTitle;
                    ProcSettingsLocal.TestRegionSize_for_minFF = ProcSettings.TestRegionSize_for_minFF;

                    %save everything
                    fnsave = [FileList{FileID,2},' - Phase1Files.mat'];
                    save(fullfile(pwd,'Z - MAT Files',fnsave),...
                        'datatable',...
                        'pos',...
                        'ProcSettingsLocal'...
                        );
                    
                    ROI_unsaved = false; % we have saved the ROI, exit the while-loop
                    
                    fprintf('\t%s','(done)');
                    fprintf('\n');

                case 'No (reselect ROI)'
                    delete(tmpax)
                    clear fcn hFH pos
            end
        end
    end
    
%save FileList in the newly created folder. This serves to mark a folder is
%'ready' for Phase2 processing of the data using the ROIs.
save(fullfile(pwd,'Phase1FileList.mat'),'FileList','OutputFolderName');
copyfile(fullfile(DataDirName,'ProcSettings.txt'),fullfile(pwd,'ProcSettings.txt'));

fprintf('\n');
fprintf('%s\n','Finished! Now you can do Phase 2.');
