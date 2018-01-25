function ProcSet = LoadProcSettings
    ProcSet = struct();
    
    if exist(fullfile(cd, 'ProcSettings.txt'), 'file') == 0
        error('Cannot find ''ProcSettings.txt'' file for your data. You''ll need one to proceed.');
    else
        fileID = fopen(fullfile(cd,'ProcSettings.txt'));
        PS_RawInput = textscan(fileID,'%s %s %s','delimiter',':');  %'delimiter','/[#|:]');
        fclose(fileID);
    end
    
    PS_Temp = cell(1,2);

    % Remove comment lines (begining with # chars)
    ProcSettIdx = 1;
    for s=1:length(PS_RawInput{1,1})
        CommentTest = PS_RawInput{1,1}(s,1);  
        if ~strcmpi(CommentTest{1,1}(1,1),'#')
            % if the whole line isn't a comment line, process it...
            PSVarName = PS_RawInput{1,1}(s,:);
            PSVarValue = PS_RawInput{1,2}(s,:);
            
            if ~isempty(strfind(char(PSVarValue),'#')) % check if variable value contains a # comment within it
                tmpPSVarValue = strsplit(char(PSVarValue),{'#','\t'},'CollapseDelimiters',true);
                PSVarValue = cellstr(tmpPSVarValue{1,1});
                clear tmpPSVarValue
            end
                        
            PS_Temp{ProcSettIdx,1} = strtrim(PSVarName);
            PS_Temp{ProcSettIdx,2} = strtrim(PSVarValue);
            ProcSettIdx = ProcSettIdx + 1;
        end      
    end
    
    % Check if what remains makes sense
    if strcmp(PS_Temp{1,1},'')
        error('I think you are using an old version of PS_Temp.txt so I can''t continue.')
    end

    % Parse the PS_Temp fields into variables.
    
    for s=1:size(PS_Temp,1)
        
        %is the setting value a number?
        SettingType = isstrprop(PS_Temp{s,2},'digit');
        if all(SettingType{1,1}) % the whole thing is a number
            valtmp = str2double(PS_Temp{s,2});
        
        % if not an integer and only one non-integer value, it could be a
        % decimal point
        elseif nnz(SettingType{1,1}) == length(SettingType{1,1}) - 1 
            valtmp = PS_Temp{s,2};
            idxtmp = find(SettingType{1,1} == 0);
            if strcmp(valtmp{1,1}(1,idxtmp),'.')
                valtmp = str2double(PS_Temp{s,2});
            end
            clear idxtmp
            
        % if it's not a number, check if we have an array bounded by []
        elseif ~isempty(regexp(PS_Temp{s,2}{1,1},'[+\w','match')) && ~isempty(regexp(PS_Temp{s,2}{1,1},'\w+]','match'))
            cell_tmp = strrep(PS_Temp{s,2}{1,1},'[','');    % delete first bracket
            cell_tmp = strrep(cell_tmp,']','');             % delete last bracket
            valtmp = str2num(cell_tmp);
            
        %if it's not a number or an array, that leaves either boolean chars
        %or simple chars
        else
            if strcmpi(PS_Temp{s,2},'True')
                valtmp = true;
            elseif strcmpi(PS_Temp{s,2},'False')
                valtmp = false;
            else
                valtmp = cell2mat(PS_Temp{s,2});
            end
            
        end

        if strcmpi(PS_Temp{s,2},'None')
        	valtmp = 0;
        end 

        vartmp = cell2mat(PS_Temp{s,1});
        ProcSet.(vartmp) = valtmp;
        
    end

    %fix data delimiter input
    if isfield(ProcSet,'DataDelimiter')
        if strcmp(ProcSet.DataDelimiter,'tab')
            ProcSet.DataDelimiter = '\t';
        elseif strcmp(ProcSet.DataDelimiter,'comma')
            ProcSet.DataDelimiter = ',';
        elseif strcmp(ProcSet.DataDelimiter,'space')
            ProcSet.DataDelimiter = ' ';
        elseif strcmp(ProcSet.DataDelimiter,'semicolon')
            ProcSet.DataDelimiter = ';';
        elseif strcmp(ProcSet.DataDelimiter,'pipe')
            ProcSet.DataDelimiter = '|';
        end
    end

    %Remangle ExcelHeaders into a cell array
    if isfield(ProcSet,'GetHeadersFromTable')
        if ~ProcSet.GetHeadersFromTable
            try
                ProcSet.ExcelHeaders = textscan(char(ProcSet.ExcelHeaders),'%s','delimiter',',');
                ProcSet.ExcelHeaders = ProcSet.ExcelHeaders{1,1}';
            catch
                errordlg('You have not specified any data headers in ProcSettings.txt file and you have not requested headers to be read from the table. Please update your ProcSettings.txt file.');
            end
        end
    end
    
disp('---------------------------------------------------------');
SettingsLoadedStr = ['Loaded ',num2str(size(PS_Temp,1)),' variables from ProcSettings file.'];
    disp(SettingsLoadedStr);
end