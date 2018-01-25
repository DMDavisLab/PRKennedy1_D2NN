ROI_stats = zeros(size(FileList,1),4);

for FileID = 1:size(FileList,1)

    Phase1Source_filename = [FileList{FileID,2},' - Phase1Files.mat'];
    load(fullfile(DataDirName,'Z - MAT Files',Phase1Source_filename));
    clear Phase1Source_filename
    
    ROIarea = polyarea(pos(:,1),pos(:,2))/1e6;
    
    ROI_stats(FileID,1) = ROIarea;
    ROI_stats(FileID,2) = size(datatable.data,1);
    
    [InROI_idx,~] = inpolygon(datatable.data(:,ProcSettings.xCoordsColumn),datatable.data(:,ProcSettings.yCoordsColumn),pos(:,1),pos(:,2));

    ROI_stats(FileID,3) = sum(InROI_idx);

    
end