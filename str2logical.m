function bool_out = str2logical(str_in)
    if any(strcmp(str_in,{'1','True','true','Yes','yes'}))
        bool_out = true;
    elseif any(strcmp(str_in,{'0','False','false','No','no'}))
        bool_out = false;
    end
end