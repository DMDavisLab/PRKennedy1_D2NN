% Similarity Clustering

    % Add UID to each point -- move this to the main Distance2Friends
    eventUID = (1:size(data_proc_corrLR,1))';
    data_proc_corrLR = horzcat(eventUID,data_proc_corrLR);
    clear eventUID
    xCol = 2;
    yCol = 3;
    ClustValsCol = 4;
       
    %find the point with the highest ilsd (inverse log sum distances, aka invlogAFD) value
    seedpoint = find(data_proc_corrLR(:,ClustValsCol) == max(data_proc_corrLR(:,ClustValsCol)));

    % round all clustering values to nearest 0.01
    RoundedClustValsCol =  size(data_proc_corrLR,2) + 1;
    data_proc_corrLR(:,RoundedClustValsCol) = ceil(data_proc_corrLR(:,ClustValsCol).*100)./100;
    
    % go through table and match similarity of xy and clusval for rows with a cluster ID == -1
    % most similar get matched to this clusterID. Unmatchable points get clusterID == 0

    % Set all events to have a non-clustered member ID of -1
    ClustIDCol = size(data_proc_corrLR,2) + 1;
    data_proc_corrLR(:,ClustIDCol) = -1;


    % Start at the seedpoint, initialise the cluster ID, and give the seedpoint the first ID
    eventID1 = seedpoint;
    ClusterID = 1;
    data_proc_corrLR(seedpoint,ClustIDCol) = ClusterID;

    InitialLocalRegionSize = 50;




% find NN

% assess similarity

% if similar give ClusterID or adopt lowest ClusterID

% next event



% for each event, if without a clusterID
while min(data_proc_corrLR(:,ClustIDCol)) == -1
    
    %Get coords of the 'focus' event (begins by using the event with the
    %highest cluster value.
    x1 = data_proc_corrLR(eventID1,xCol);
    y1 = data_proc_corrLR(eventID1,yCol);

    % crop to region around this focus event. %TODO - add check for sufficient neighbours?
    LocalRegionBounds = [x1-(0.55 * InitialLocalRegionSize) x1+(0.55 * InitialLocalRegionSize) y1-(0.55 * InitialLocalRegionSize) y1+(0.55 * InitialLocalRegionSize)];
    LocalRegion = RegionCropper(data_proc_corrLR,LocalRegionBounds,[xCol yCol]);
    
    
        if size(LocalRegion,1) == 1
            %everything is already too far away and this is a lonely event :( SAD FACE
            data_proc_corrLR(LocalRegion(1,1),ClustIDCol) = 0;
            dregs = find(data_proc_corrLR(:,ClustIDCol) == -1);
            if ~isempty(dregs)
                eventID1 = data_proc_corrLR(dregs(randi(size(dregs,1),1)),1); % New eventID is a random event from the list of unassigned events.
                ClusterID = ClusterID + 1;
            end
        else
            % we have neighbours, let's check them out
            for eventID2 = 1:size(LocalRegion,1)
                x2 = LocalRegion(eventID2,xCol);
                y2 = LocalRegion(eventID2,yCol);

                deltaX = x1-x2;
                deltaY = y1-y2;

                EucDistSqrd = (deltaX*deltaX) + (deltaY*deltaY);
                LocalRegion(eventID2,ClustIDCol+1) = EucDistSqrd;
            end

    %         % Delete points beyond the Initial Region Size reach (makes search region circular)
    %         LocalRegion(LocalRegion(:,ClustIDCol+1)>InitialLocalRegionSize,:) = [];

            % Sort the remaining points by distance
            LocalRegion = sortrows(LocalRegion,ClustIDCol+1);

    %         mean_n10 = mean(RegionTemp(1:MaxNeighbours,3)); % Mean Including self
    %         twostdev_n10 = 2 * std(RegionTemp(1:MaxNeighbours,3)); % StdDev Including self

            RegionTemp = LocalRegion(1:11,:); % delete self
            figure
            scatter(LocalRegion(:,7),LocalRegion(:,4));
            axis([0 20 0 5]);
            
            if (LocalRegion(1,4) - LocalRegion(2,4)) / LocalRegion(1,4) < 0.05
                % the NN event is 95% similar to the query event
            end
        
            
            if size(find(LocalRegion(2:end,ClustIDCol) == -1),1) == size(LocalRegion,1)-1
                % all neighbours (excluding self) are unassigned
                data_proc_corrLR(LocalRegion(:,1),ClustIDCol) = ClusterID; % Give all the points the same cluster ID
                eventID1 = LocalRegion(end,1); % Set the new eventID to the one furthest away
            elseif size(find(LocalRegion(:,ClustIDCol) == -1),1) == 0
                % all neighbours assigned to clusters already
                dregs = find(data_proc_corrLR(:,ClustIDCol) == -1);
                if ~isempty(dregs)
                    eventID1 = data_proc_corrLR(dregs(randi(size(dregs,1),1)),1); % New eventID is a random event from the list of unassigned events.
                    ClusterID = ClusterID + 1;
                end
            else
                % only some neighbours are already in clusters so we should
                % inherit their IDs
                CluIDtmp = mode(LocalRegion(:,ClustIDCol));
                % give all others this ID
                data_proc_corrLR(LocalRegion(:,1),ClustIDCol) = CluIDtmp;
            end
        end

%         EventCount = EventCount + 1;
end

