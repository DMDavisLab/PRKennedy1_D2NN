function [ppFFD_tmp,ppAFD_tmp,ppFFA_tmp,ppAFA_tmp,ppDistancesToN,ppAnglesToN] = DTF2ParaFunc3(data_query,ProcSettings,data_ref)

% data_query = coordinates table to which distance measurements will originate
% data_ref = coordinates table to which distance measurements will terminate

    if nargin < 3
        data_ref = data_query;
        InfoMessage = ['Examining ',num2str(size(data_query,1)),' events in 1 channel:'];
    else
        InfoMessage = ['Examining ',num2str(size(data_query,1)),' (subject) and ',num2str(size(data_ref,1)),' (reference) events in 2 channels:'];
    end
    disp(InfoMessage);
        
    % Start a parpool for different versions of MATLAB
    labversion = strsplit(version,{'(',')'},'CollapseDelimiters',true); 
    if strcmp(labversion(2),'R2014a') || strcmp(labversion(2),'R2014b')
        parpool;
    elseif strcmp(labversion(2),'R2013a') || strcmp(labversion(2),'R2013b')
        DistanceMeasurers = matlabpool;
    end

    ppDistancesToN = NaN(ProcSettings.FurthestFriendID,size(data_query,1));
    ProcInit = tic;
% % %
% % % Check for (and break up) large data sets.
% % %
        
    if size(data_query,1) > 1000000
        
        %add a UID to data for each row/event so we can return the data to its rightful row.
        data_query(:,end+1) = 1:size(data_query);
        
        % find the number of segments to split up the image into for approximately this many events (500k, assumes even distribution)
        number_of_chunks = ceil(sqrt(size(data_query,1)/500000))^2; 
        chunk_width = (ceil((ProcSettings.PlotWidth / sqrt(number_of_chunks))/100))*100; % adjust the width of the chunks to the nearest 100.
        
        InfoMessage = ['Large data set will be processed as ',num2str(number_of_chunks),' chunks.'];
        disp(InfoMessage);
        
        %meshgrid of xy starting coords for each sub-region
        [X_query, Y_query] = meshgrid(ProcSettings.AxisLimits(1):chunk_width:ProcSettings.AxisLimits(2)-chunk_width,ProcSettings.AxisLimits(3):chunk_width:ProcSettings.AxisLimits(4)-chunk_width);
        subregion_query = zeros(number_of_chunks,4);
        for chunkID = 1:number_of_chunks
            subregion_query(chunkID,1:2) = [X_query(chunkID),Y_query(chunkID)];
            subregion_query(chunkID,3:4) = subregion_query(chunkID,1:2) + chunk_width;
        end
        
        subregion_ref = subregion_query;
        subregion_ref(:,1) = subregion_ref(:,1) - (0.1*chunk_width);
        subregion_ref(:,2) = subregion_ref(:,2) - (0.1*chunk_width);
        subregion_ref(:,3) = subregion_ref(:,3) + (0.1*chunk_width);
        subregion_ref(:,4) = subregion_ref(:,4) + (0.1*chunk_width);
        
        % initialise the output tables with NaNs to receive measurements as they are processed
        ppFFD_tmp = NaN(size(data_query,1),1);
        ppAFD_tmp = NaN(size(data_query,1),1);
        for chunkID = 1:number_of_chunks
            
            InfoMessage = ['Processing chunk ',num2str(chunkID),' of ',num2str(number_of_chunks),'.'];
            disp(InfoMessage);
            
            Query_CropBounds = [subregion_query(chunkID,1),subregion_query(chunkID,3),subregion_query(chunkID,2),subregion_query(chunkID,4)];
            Chunk_query = RegionCropper(data_query,Query_CropBounds,[1,2]);        
            
            if size(Chunk_query,1) > 0
                Ref_CropBounds = [subregion_ref(chunkID,1),subregion_ref(chunkID,3),subregion_ref(chunkID,2),subregion_ref(chunkID,4)];
                Chunk_ref = RegionCropper(data_ref,Ref_CropBounds,[1,2]);
                
            % Adaptive the starting local-region size to stand a chance at getting a goldilocks number of events on the first attempt
                ref_data_density = size(Chunk_ref,1) / (ProcSettings.ImageSize*ProcSettings.ImageSize);
                ProcSettings.TestRegionSize_for_minFF = 1.2 * ceil(sqrt(ProcSettings.FurthestFriendID / ref_data_density));
                
                [ppFFD_chunk, ppAFD_chunk, DistsToN_chunk] = ppGetDistances(Chunk_query,Chunk_ref,ProcSettings);
                
            % return results to approriate indices
                ppFFD_tmp(Chunk_query(:,3)) = ppFFD_chunk;
                ppAFD_tmp(Chunk_query(:,3)) = ppAFD_chunk;
                ppDistancesToN(:,Chunk_query(:,3)) = DistsToN_chunk;
            end
            
            clear ppFFD_chunk ppAFD_chunk DistsToN_chunk Chunk_query Chunk_ref ref_data_density
            
        end

        % find nan cells that were missed - generally these are exactly on the sub-region boundaries
        missed_idx = data_query(isnan(ppAFD_tmp),3);
        if ~isempty(missed_idx)
        %process the missing events.
            InfoMessage = ['Processing residual events (',num2str(numel(missed_idx)),').'];
            disp(InfoMessage);
            Chunk_query = data_query(missed_idx,:);
            [ppFFD_chunk, ppAFD_chunk, DistsToN_chunk] = ppGetDistances(Chunk_query,data_ref,ProcSettings);
        % return results to approriate indices
            ppFFD_tmp(Chunk_query(:,3)) = ppFFD_chunk;
            ppAFD_tmp(Chunk_query(:,3)) = ppAFD_chunk;
            ppDistancesToN(:,Chunk_query(:,3)) = DistsToN_chunk;
            clear ppFFD_chunk ppAFD_chunk DistsToN_chunk Chunk_query missed_events
        end
        
%         data_query(:,end) = []; % delete the UID column ... or don't ... this isn't sent back
        clear missed_idx subregion_query subregion_ref number_of_chunks chunkID
        
    else
        % Do it the normal way
        ref_data_density = size(data_query,1) / (ProcSettings.ImageSize*ProcSettings.ImageSize);
        ProcSettings.TestRegionSize_for_minFF = ceil(sqrt(ProcSettings.FurthestFriendID / ref_data_density));

        [ppFFD_tmp,ppAFD_tmp,ppFFA_tmp,ppAFA_tmp,ppDistancesToN,ppAnglesToN] = ppGetDistances(data_query,data_ref,ProcSettings);     
        
    end

    ProcTerm = toc(ProcInit);
    proc_hours = floor(ProcTerm / 3600);
    ProcTerm = ProcTerm -  proc_hours * 3600;
    proc_mins = floor(ProcTerm / 60);
    proc_secs = floor(10 * (ProcTerm -  proc_mins * 60)) / 10;
    EndMessage = ['Table examined in ',num2str(proc_hours),' hours, ',num2str(proc_mins),' minutes, and ',num2str(proc_secs),' seconds (',num2str(floor((size(data_query,1) / ProcTerm))),' events/s).'];
    disp(EndMessage);

    
    % Close a parpool for different versions of MATLAB
    if strcmp(labversion(2),'R2014a') || strcmp(labversion(2),'R2014b')
        delete(gcp)
    elseif strcmp(labversion(2),'R2013a') || strcmp(labversion(2),'R2013b')
        delete(DistanceMeasurers)
    end

end




function [FFD, AFD, FFA, AFA, DTN, ATN] = ppGetDistances(dq,dr,ps)
    fprintf([9,repmat('.',1,50),'\n',9,'\n']);
    ProgressStepping = floor(size(dq,1) / 50);
    xCol = ps.xCol_procd;
    yCol = ps.yCol_procd;
    xData = dq(:,xCol); % splitting the data into separate arrays helps to
    yData = dq(:,yCol); %   reduce communications overhead in parfor loop.
    parfor eventID1 = 1:size(dq,1);
        x1 = xData(eventID1);
        y1 = yData(eventID1);
        if isfield(ps,'TestRegionSize_for_minFF')
            ppTestRegionSize = ps.TestRegionSize_for_minFF;
        else
            ppTestRegionSize = 15;
        end
        UseableRegion = false; % TestRegion is considered "unsuitable" until it contains enough neighbours to be useful to us.
        while ~UseableRegion 
            TestRegionBounds = [x1-ppTestRegionSize x1+ppTestRegionSize y1-ppTestRegionSize y1+ppTestRegionSize];
            TestRegion = RegionCropper(dr,TestRegionBounds,[xCol,yCol]);
            if size(TestRegion,1) >= ps.FurthestFriendID + 1
                UseableRegion = true; % TestRegion contains self + N-neighbours and is useable
            else
                ppTestRegionSize = 2*ppTestRegionSize; % Insufficient events: enlarge the region and try again
            end
        end
        
        % measure all the distances
        ppGeometryTemp = zeros(size(TestRegion,1),4);  
        for eventID2 = 1:size(TestRegion,1)
            x2 = TestRegion(eventID2,xCol);
            y2 = TestRegion(eventID2,yCol);
            deltaX = x1-x2;
            deltaY = y1-y2;
            
            ppEucDist = sqrt((deltaX*deltaX) + (deltaY*deltaY));
            ppGeometryTemp(eventID2,1) = ppEucDist;
            
            ppAngle = atan2(deltaY, deltaX); % angle in radians relative to horizontal axis
            ppGeometryTemp(eventID2,2) = ppAngle;
            
            ppGeometryTemp(eventID2,3) = x2;
            ppGeometryTemp(eventID2,4) = y2;
        end
        
        ppGeometryTemp = sortrows(ppGeometryTemp,1);
%         ppGeometryTemp(1,:) = []; %delete the distance-to-self
        ppGeometryTemp(ppGeometryTemp(:,1) == 0,:) = []; %delete the distance-to-self

        % collect distances and angles for this event to NNs
        DTN(:,eventID1) = ppGeometryTemp(1:ps.FurthestFriendID,1); % Distances
        ATN(:,eventID1) = ppGeometryTemp(1:ps.FurthestFriendID,2); % Angles
        
        FFD(eventID1,1) = ppGeometryTemp(ps.FurthestFriendID,1);        % Furthest NN distance
        AFD(eventID1,1) = sum(ppGeometryTemp(1:ps.FurthestFriendID,1)); % Sum of distances
                
        FFA(eventID1,1) = ppGeometryTemp(ps.FurthestFriendID,2);        % Furthest NN angle
%               
%         % weighted mean of neighbours, including self
%         dist_weights = ones(1,ps.FurthestFriendID+1)./(1:ps.FurthestFriendID+1);
%         avg_x2 = sum(([0;(ppGeometryTemp(1:ps.FurthestFriendID,3)-x1)].*dist_weights'))/ps.FurthestFriendID;
%         avg_y2 = sum(([0;(ppGeometryTemp(1:ps.FurthestFriendID,4)-y1)].*dist_weights'))/ps.FurthestFriendID; % computer average of the unit vectors of NNs and calculate its angle
% 
%         % weighted mean of neighbours, excluding self
%         dist_weights = ones(1,ps.FurthestFriendID)./(1:ps.FurthestFriendID);
%         avg_x2 = sum(((ppGeometryTemp(1:ps.FurthestFriendID,3)-x1).*dist_weights'))/ps.FurthestFriendID;
%         avg_y2 = sum(((ppGeometryTemp(1:ps.FurthestFriendID,4)-y1).*dist_weights'))/ps.FurthestFriendID; % computer average of the unit vectors of NNs and calculate its angle

        % unweighted mean angle, excluding self
        avg_x2 = mean(ppGeometryTemp(1:ps.FurthestFriendID,3)-x1);
        avg_y2 = mean(ppGeometryTemp(1:ps.FurthestFriendID,4)-y1);

        AFA(eventID1,1) = atan2(avg_y2,avg_x2);         % Average NN angle, for degrees just use (atan2(avg_y2,avg_x2)*180)/pi;

%                 figure;
%                 heatlines = [linspace(1,0,10)',zeros(10,1),linspace(0,1,10)'];
%                 for n=1:ps.FurthestFriendID
%                     line([x1,ppGeometryTemp((ps.FurthestFriendID-n)+1,3)],[y1,ppGeometryTemp((ps.FurthestFriendID-n)+1,4)],'Color',heatlines((ps.FurthestFriendID-n)+1,:),'LineWidth',1);
%                 end
%                 hold on
%                 scatter(TestRegion(:,1),TestRegion(:,2),10,'ok','filled');
%                 scatter(x1,y1,10,'ok','filled');
%                 scatter(x1,y1,100,'or');
%                 scatter(ppGeometryTemp(1:10,3),ppGeometryTemp(1:10,4),100,'og')
%                 scatter(x1+avg_x2,y1+avg_y2,10,'dm','filled')
%                 line([x1,x1+avg_x2],[y1,y1+avg_y2],'LineStyle',':','Color','m','LineWidth',1);
%                 axis  equal



        % Do randomisation of the data within a certain distance of the key event
        % local_rand_region = randomise(within large_r)
        % distances for the smaller randomised data
%         [ppFFD_tmp,ppAFD_tmp,ppFFA_tmp,ppAFA_tmp,ppDistancesToN,ppAnglesToN] = ppGetDistances(data_query,data_ref,ProcSettings);




        if rem(eventID1,ProgressStepping) == 0
            fprintf('\b|\n');  %update tracking bar
        end      
    end % end parfor
end


        
