% Generate randomised data field containing different densities
BoundaryAvoidance = 0.1; % keep random events this far from the edge to avoid crashing.
ImageBounds = [0 18 0 18]; % in um

% density_options = [...  % events per um^2, typical cell (~2500 events/um)
%     1, ...
%     3, ...
%     10, ...
%     30, ...
%     100, ...
%     300, ...
%     1000, ...
%     3000, ...
%     10000, ...          %
%     30000, ...          % typical cell ranges
%     100000, ...         %
%     300000, ...
%     1000000 ...         % one event per nm^2
%     ];

density_options = [...  % events per um^2, typical cell (~2500 events/um)
    1,1,1,1,1,1,1, ...
    3,3,3,3,3,3,3, ...
    10,10,10,10,10,10,10, ...
    30,30,30,30,30,30,30, ...
    100,100,100,100,100,100,100, ...
    300,300,300,300,300,300,300, ...
    1000,1000,1000,1000,1000,1000,1000, ...
    3000,3000,3000,3000,3000,3000,3000, ...
    10000,10000,10000,10000,10000,10000,10000, ...          %
    30000,30000,30000,30000,30000,30000,30000 ...          % typical cell ranges
    ];

pos = [1000,1000;1000,17000;17000,17000;17000,10000];

home %clean up the command window
rng('shuffle')  % reset the random seed
image_area = (ImageBounds(2)-ImageBounds(1))*(ImageBounds(4)-ImageBounds(3)); % um^2

datatable = struct;

DataDirName = uigetdir(pwd,'Choose your data folder');

if DataDirName ~= 0
    cd(DataDirName);
end
alphanums = ['a':'z' 'A':'Z' '0':'9'];
randname = alphanums(randi(numel(alphanums),[1 5]));
OutputFolderName = ['D2F2-(',randname,')'];
mkdir(fullfile(DataDirName,OutputFolderName));
cd(fullfile(DataDirName,OutputFolderName));
mkdir(fullfile(pwd,'Z - MAT Files'));


for d = 1:size(density_options,2)
   
    total_events = image_area * density_options(d);
    
    tmp_image_table = zeros(total_events,2); 
    tmp_image_table(:,1) = randi([(ImageBounds(1)*10000), (ImageBounds(2)*10000)],total_events,1);
    tmp_image_table(:,2) = randi([(ImageBounds(3)*10000), (ImageBounds(4)*10000)],total_events,1);
    tmp_image_table = tmp_image_table / 10; % convert values back to nm.
    
    % fix x boundary events
    badidx = find(tmp_image_table(:,1)==0);
    tmp_image_table(badidx,1) = 0.1;
    badidx = find(tmp_image_table(:,1)==ImageBounds(2)*10000);
    tmp_image_table(badidx,1) = (ImageBounds(2)*10000) - 0.1;

    % fix y boundary events
    badidx = find(tmp_image_table(:,2)==0);
    tmp_image_table(badidx,2) = 0.1;
    badidx = find(tmp_image_table(:,2)==ImageBounds(4)*10000);
    tmp_image_table(badidx,2) = (ImageBounds(4)*10000) - 0.1;
    
    datatable.data = tmp_image_table;
    datatable.colheaders = {'X (nm)','Y (nm)'};
    
    table_UID = alphanums(randi(numel(alphanums),[1 5]));

    ProcSettingsLocal = struct;
    ProcSettingsLocal.ExptTitle = ['Random Density=',num2str(density_options(d)),' per um2 (UID ',table_UID,')'];
    ProcSettingsLocal.AxisLimits = ImageBounds * 1000; % in nm
    ProcSettingsLocal.PlotWidth = ProcSettingsLocal.AxisLimits(2);
    ProcSettingsLocal.SaveHighDPI = strcat('-r',num2str(ProcSettingsLocal.PlotWidth / 10)); % for 10 nm per pixel

    fnsave = [ProcSettingsLocal.ExptTitle,' - Phase1Files.mat'];
    save(fullfile(pwd,'Z - MAT Files',fnsave),...
        'datatable',...
        'pos',...
        'ProcSettingsLocal'...
        );
   
    tmp_image_fname = [ProcSettingsLocal.ExptTitle,'.txt'];
    dlmwrite(tmp_image_fname,tmp_image_table,'Delimiter','\t');
    
end

% get list of the newly created files
DirFullList = dir(fullfile(pwd,'Z - MAT Files'));                              % Get the data for the current directory
FolderIndex = [DirFullList.isdir];              % Find the index for directories
FileList = {DirFullList(~FolderIndex).name}';   % Get a list of the files
for f = 1:length(FileList)
    FileList{f,1} = strrep(FileList{f,1},' - Phase1Files.mat','.txt');
    [~, fname, ext] = fileparts(FileList{f,1});
    FileList{f,2} = fname;
    FileList{f,3} = ext;
end
[~, natsortidx] = sort_nat(FileList(:,2),'Ascend');
FileList = FileList(natsortidx,:);
            
% make a phase1files list
save(fullfile(pwd,'Phase1FileList.mat'),'FileList','OutputFolderName');

% make a procsettings
ProcSettingsGenerated = cell(20,3);

ProcSettingsGenerated(1,:) = {'xCoordsColumn',':','1'};
ProcSettingsGenerated(2,:) = {'yCoordsColumn',':','2'};
ProcSettingsGenerated(3,:) = {'ChannelIDColumn',':','None'};
ProcSettingsGenerated(4,:) = {'FurthestFriendID',':','10'};
ProcSettingsGenerated(5,:) = {'DataTableFileExt',':','txt'};
ProcSettingsGenerated(6,:) = {'DataDelimiter',':','tab'};
ProcSettingsGenerated(7,:) = {'DataTableScale',':','1'};
ProcSettingsGenerated(8,:) = {'ImageSize',':','18000'};
ProcSettingsGenerated(9,:) = {'SaveClusterColourImages',':','True'};
ProcSettingsGenerated(10,:) = {'DoStrictdNN',':','False'};
ProcSettingsGenerated(11,:) = {'DoSumdNN',':','False'};
ProcSettingsGenerated(12,:) = {'DoInvLogSumdNN',':','True'};
ProcSettingsGenerated(13,:) = {'NormaliseToDensity',':','True'};
ProcSettingsGenerated(14,:) = {'DoThresholds',':','true'};
ProcSettingsGenerated(15,:) = {'SaveThresholdImages',':','true'};
ProcSettingsGenerated(16,:) = {'ColourEventsOver',':','[241 215 0]'};
ProcSettingsGenerated(17,:) = {'ColourEventsUnder',':','[46 146 166]'};
ProcSettingsGenerated(18,:) = {'ColourBackground',':','[0 0 0]'};
ProcSettingsGenerated(19,:) = {'FixedThreshold',':','1.0'};
ProcSettingsGenerated(20,:) = {'SaveMATFiles',':','true'};

    fileID = fopen(fullfile(pwd,'ProcSettings.txt'),'w');
    fprintf(fileID,'%s','# Automatically Generated ProcSettings Files for CSR Data');
    fprintf(fileID,'\r\n');
    formatSpec = '%s\t%s\t%s\r\n';
    nrows = size(ProcSettingsGenerated,1);
    for row = 1:nrows
        fprintf(fileID,formatSpec,ProcSettingsGenerated{row,:});
    end
    fclose(fileID);
