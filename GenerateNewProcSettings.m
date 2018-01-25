% Make a ProcSettings File
    prompt = {...
                'xCoordsColumn',...
                'yCoordsColumn',...
                'ChannelIDColumn',...
                'FurthestFriendID',...
                'DataTableFileExt',...
                'DataDelimiter',...
                'DataTableScale',...
                'ImageSize',...
                'SaveClusterColourImages',...
                'DoStrictdNN',...
                'DoSumdNN',...
                'DoInvLogSumdNN',...
                'NormaliseToDensity',...
                'DoThresholds',...
                'SaveThresholdImages',...
                'ColourEventsOver',...
                'ColourEventsUnder',...
                'ColourBackground',...
                'FixedThreshold',...
                'SaveMATFiles',...
              };
      
    dlg_title = 'Check settings to apply...';
    num_lines = 1;
    
    defaults = {...
        '3',...
        '4',...
        'None',...
        '10',...
        'csv',...
        'comma',...
        '1',...
        '26932',...
        'True',...
        'True',...
        'True',...
        'True',...
        'True',...
        'True',...
        'True',...
        '[241 215 0]',...
        '[46 146 166]',...
        '[0 0 0]',...
        '0.5',...
        'True',...
       };

    answer = inputdlg(prompt,dlg_title,num_lines,defaults);

    if ~isempty(answer)

        fileID = fopen('ProcSettings.txt','w');
        fprintf(fileID,'%s\n',['# ProcSettings for D2NN - Generated ',datestr(date)]);
        for p = 1:length(prompt)
            varlength = floor(length(prompt{p})/4);
            
            if varlength == 5
                formatspec = '%s\t%s\t%s\n';
            elseif varlength == 4
                formatspec = '%s\t\t%s\t%s\n';
            elseif varlength == 3
                formatspec = '%s\t\t\t%s\t%s\n';
            else
                formatspec = '%s\t\t\t\t%s\t%s\n';
            end
            fprintf(fileID,formatspec,prompt{p},':',answer{p});
        end
        fclose(fileID);
    else
        error('Cancelled?! So rude.');
    end
