data = data_proc;                           % if your data is already in 'data' array then don't run this line

data_length = length(data(:,1));            % Get the length of the input vectors.
data_rand_idx = randperm(data_length);      % Create a vector of random integers from 1 to len

% make new tables containing a range of randomly-picked events from 100% to 1% of the original data
data_proc_100 = data;           
data_proc_90 = data(data_rand_idx(1:floor(0.90*data_length)),:);
data_proc_75 = data(data_rand_idx(1:floor(0.75*data_length)),:);
data_proc_50 = data(data_rand_idx(1:floor(0.50*data_length)),:);
data_proc_25 = data(data_rand_idx(1:floor(0.25*data_length)),:);
data_proc_10 = data(data_rand_idx(1:floor(0.10*data_length)),:);
data_proc_05 = data(data_rand_idx(1:floor(0.05*data_length)),:);
data_proc_01 = data(data_rand_idx(1:floor(0.01*data_length)),:);

% Get density within ROI
ROIarea = polyarea(pos(:,1),pos(:,2));
[in_roi,on_roi] = inpolygon(data_proc(:,1),data_proc(:,2),pos(:,1),pos(:,2));
density = sum(in_roi) / ROIarea;

% seed events randomly within whole image area to match density (or 2x, 3x density etc)

% delete events outside ROI


