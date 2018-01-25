function CroppedRegion = RegionCropper(FullRegion,RegionBounds,CoordColumns)
    %% cropping
    % FullRegion = table with XY data in columns specified by CoordColumns
    % RegionBounds = [minX maxX minY maxY]
    % CoordColumns = [xColumn yColumn]

    % FullRegion = data2;
    % RegionBounds = [region_xmin region_xmax region_ymin region_ymax];
    % CoordColumns = [xCoordsColumn yCoordsColumn];

%     if nargin < 4
%         keep_on_boundary_events = false;
%     end

    % if length(CoordColumns) == 2 && size(RegionBounds,2) ~= 4
    %     error(message('You must supply min/max values for both dimensions'));
    % end

    
    if length(CoordColumns) == 2 % 2D data

        CroppedRegion=FullRegion(FullRegion(:,CoordColumns(1)) > RegionBounds(1) & FullRegion(:,CoordColumns(1)) < RegionBounds(2) & ...
                                 FullRegion(:,CoordColumns(2)) > RegionBounds(3) & FullRegion(:,CoordColumns(2)) < RegionBounds(4),:);

    elseif length(CoordColumns) == 3 % 3D data

        CroppedRegion=FullRegion(FullRegion(:,CoordColumns(1)) > RegionBounds(1) & FullRegion(:,CoordColumns(1)) < RegionBounds(2) & ...
                                 FullRegion(:,CoordColumns(2)) > RegionBounds(3) & FullRegion(:,CoordColumns(2)) < RegionBounds(4) & ...
                                 FullRegion(:,CoordColumns(3)) > RegionBounds(5) & FullRegion(:,CoordColumns(3)) < RegionBounds(6),:);

    elseif length(CoordColumns) == 1 % 1D data

        CroppedRegion=FullRegion(FullRegion(:,CoordColumns(1)) > RegionBounds(1) & FullRegion(:,CoordColumns(1)) < RegionBounds(2),:);

    end 
    
end
