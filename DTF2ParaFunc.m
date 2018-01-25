function [ppFFD_tmp,ppAFD_tmp,ppDistancesToN] = DTF2ParaFunc(data,vars,data2)

% data = reference coordinates table
% data2 = query coordinates table
xCoordsColumn = vars(1);
yCoordsColumn = vars(2);
FurthestFriendID = vars(3);

if nargin < 3
    data2 = data;
    InfoMessage = ['Examining ',num2str(size(data,1)),' events in 1 channel:'];
else
    InfoMessage = ['Examining ',num2str(size(data,1)),' (subject) and ',num2str(size(data2,1)),' (reference) events in 2 channels:'];
end

ProgressStepping = floor(size(data,1) / 50);

% Start a parpool for different versions of MATLAB
labversion = strsplit(version,{'(',')'},'CollapseDelimiters',true); 
if strcmp(labversion(2),'R2014a') || strcmp(labversion(2),'R2014b')
    parpool;
elseif strcmp(labversion(2),'R2013a') || strcmp(labversion(2),'R2013b')
    DistanceMeasurers = matlabpool;
end

ppDistancesToN = zeros(FurthestFriendID,size(data,1));

disp(InfoMessage);
% fprintf(['\n' repmat('.',1,m) '\n\n']);
fprintf([9,repmat('.',1,50),'\n',9,'\n']);

parfor_stopwatch = tic;

parfor eventID1 = 1:size(data,1);
    x1 = data(eventID1,xCoordsColumn);
    y1 = data(eventID1,yCoordsColumn);

    if numel(vars) == 4
        ppTestRegionSize = vars(4);
    else
        ppTestRegionSize = 15;
    end

    UseableRegion = false; % TestRegion is unsuitable until it has been checked for suitability (below)
   
    while ~UseableRegion
        
        TestRegionBounds = [x1-ppTestRegionSize x1+ppTestRegionSize y1-ppTestRegionSize y1+ppTestRegionSize];
        TestRegion = RegionCropper(data2,TestRegionBounds,[xCoordsColumn yCoordsColumn]);
        
        if size(TestRegion,1) >= FurthestFriendID + 1 % plus 1 to include self
            UseableRegion = true;
        else
            ppTestRegionSize = 2*ppTestRegionSize; % If there aren't enough points then enlarge the region and try again
        end
        
    end

    ppGeometryTemp = zeros(size(TestRegion,1),1);
    
    for eventID2 = 1:size(TestRegion,1)
        x2 = TestRegion(eventID2,xCoordsColumn);
        y2 = TestRegion(eventID2,yCoordsColumn);
        
        deltaX = x1-x2;
        deltaY = y1-y2;
        
        ppEucDist = sqrt((deltaX*deltaX) + (deltaY*deltaY));
        ppGeometryTemp(eventID2,1) = ppEucDist;
    end
       
    ppGeometryTemp = sort(ppGeometryTemp(:,1));
    ppGeometryTemp(1,:) = []; %delete the distance-to-self
    
    ppDistancesToN(:,eventID1) = ppGeometryTemp(1:FurthestFriendID,1);
%     ppIDsToN(:,eventID1) = ppGeometryTemp(1:vars.FurthestFriendID,1);
    
    ppFurthestFriendDistance = ppGeometryTemp(FurthestFriendID,1);

    ppAllFriendDistance = sum(ppGeometryTemp(1:FurthestFriendID,1));
    
    ppFFD_tmp(eventID1,1) = ppFurthestFriendDistance;
    ppAFD_tmp(eventID1,1) = ppAllFriendDistance;

%     data(eventID1,vars.FFDCol) = FurthestFriendDistance;
%     data(eventID1,vars.AFDCol) = AllFriendDistance;
    
    if rem(eventID1,ProgressStepping) == 0
        fprintf('\b|\n');

%         %Display some information
%         t = getCurrentTask();
%         TrackingMessage=[datestr(fix(clock),'HH:MM:SS'),9,'CPU-' num2str(t.ID),' completed event-row ',num2str(eventID1),' in ',vars.ExptTitle,'.'];
%         disp(TrackingMessage);

    end
        
end

proc_time = toc(parfor_stopwatch);
proc_hours = floor(proc_time / 3600);
proc_time = proc_time -  proc_hours * 3600;
proc_mins = floor(proc_time / 60);
proc_secs = floor(10 * (proc_time -  proc_mins * 60)) / 10;

EndMessage = ['Table examined in ',num2str(proc_hours),' hours, ',num2str(proc_mins),' minutes, and ',num2str(proc_secs),' seconds (',num2str(floor((size(data,1) / proc_time))),' events/s).'];
disp(EndMessage);

% Close a parpool for different versions of MATLAB
if strcmp(labversion(2),'R2014a') || strcmp(labversion(2),'R2014b')
    delete(gcp)
elseif strcmp(labversion(2),'R2013a') || strcmp(labversion(2),'R2013b')
    delete(DistanceMeasurers)
end

end