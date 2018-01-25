
% Useage
%
% Beforehand:
%    SaveToFolder = 'D - Threshold','Mean plus StDev';
%    SaveFileName = [FileList{FileID,2},' - events over threshold=',num2str(Cluster_thr_meanplus1stdev),' (Mean - StDev).png'];
%    Cluster_thr_meanplus1stdev = mean(data_proc(:,ProcSettings.invlogAFDCol)) + std(data_proc(:,ProcSettings.invlogAFDCol));
%    RefDataCol = ProcSettings.invlogAFDCol;
%
% Calling the function:
%   [num_events_over, num_events_under] = thresholder(data_proc,Cluster_thr_meanplus1stdev,SaveToFolder,SaveFileName,ProcSettings)
% 
% Afterwards:
% 
%         percent_in_meanminusstdev = (num_events_over/total_events_table) * 100;
%         
%         % update the summary table
%         stats_mmsd(1,FileID) = Cluster_thr_meanplus1stdev;
%         stats_mmsd(2,FileID) = num_events_over
%         stats_mmsd(3,FileID) = total_events_table;
%         stats_mmsd(4,FileID) = percent_in_meanminusstdev;

function [events_over, events_under] = thresholder(inputdata,thresholdvalue,RefDataCol,outputfolder,outputname,ProcSettings)

    if ~isdir(fullfile(pwd,outputfolder))
        mkdir(fullfile(pwd,outputfolder));
    end
    
    % Check ProcSettings color settings
    if any(horzcat(ProcSettings.ColourEventsOver > 1,ProcSettings.ColourEventsUnder > 1,ProcSettings.ColourBackground > 1))
        % correct the 8 bit RGB colours to fractional
        ProcSettings.ColourEventsOver = ProcSettings.ColourEventsOver ./ 255;
        ProcSettings.ColourEventsUnder = ProcSettings.ColourEventsUnder ./ 255;
        ProcSettings.ColourBackground = ProcSettings.ColourBackground ./ 255;
    end


    % Cluster threshold (mean minus 1 stdev)
    events_over = inputdata(inputdata(:,RefDataCol) > thresholdvalue,:);
    events_under = inputdata(inputdata(:,RefDataCol) <= thresholdvalue,:);
    
    % num_events_over = size(events_over,1);
%     num_events_under = size(plot_pts_under,1); 
    
    if ProcSettings.SaveThresholdImages
        
        %Scramble over threshold indicies
        data_over_len = length(events_over(:,1));     %Get the length of the input vectors.
        data_over_rand_ind=randperm(data_over_len);     %Create a vector of random integers from 1 to len
        data_over_randx=events_over(data_over_rand_ind,ProcSettings.xCol_procd);  %Randomize the input vectors.
        data_over_randy=events_over(data_over_rand_ind,ProcSettings.yCol_procd);  %Randomize the input vectors.

        %Scramble under threshold indices
        data_under_len = length(events_under(:,1));           %Get the length of the input vectors.
        data_under_rand_ind=randperm(data_under_len);           %Create a vector of random integers from 1 to len
        data_under_randx=events_under(data_under_rand_ind,ProcSettings.xCol_procd);   %Randomize the input vectors.
        data_under_randy=events_under(data_under_rand_ind,ProcSettings.yCol_procd);   %Randomize the input vectors.

        % Plot the under-points then the over-points
        [fig_plot_pts, axes_plot_pts] = DoMeAFigure(ProcSettings.AxisLimits,ProcSettings.ColourBackground);
                
        scatter(data_under_randx,data_under_randy,1,ProcSettings.ColourEventsUnder,'.');
        scatter(data_over_randx,data_over_randy,1,ProcSettings.ColourEventsOver,'.');

        cbar = colormap;
        if isempty(events_over(:,RefDataCol))
            thold_posn = size(cbar,1);
        else
            thold_posn = floor(thresholdvalue / max(events_over(:,RefDataCol)) * size(cbar,1));
        end
        cbar(1:thold_posn,:) = repmat(ProcSettings.ColourEventsUnder, thold_posn,1);
        cbar(thold_posn+1:end,:) = repmat(ProcSettings.ColourEventsOver, size(cbar,1) - thold_posn,1);
        colormap(cbar);
        caxis([0 max(inputdata(:,RefDataCol))]);
        
        set(axes_plot_pts,'Color',ProcSettings.ColourBackground);
        set(fig_plot_pts,'InvertHardcopy','off')
        
        SaveColourBar(cbar,caxis,fullfile(pwd,outputfolder,[outputname,'_colourbar']));
        
        print(fig_plot_pts,'-dpng',ProcSettings.SaveHighDPI,fullfile(pwd,outputfolder,[outputname,'.png']));
        close(fig_plot_pts);
               
        clear axes_plot_pts
        
    end
    
end