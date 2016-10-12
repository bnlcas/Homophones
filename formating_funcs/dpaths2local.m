function evnt_out = dpaths2local(evnt, root_dir)
%% This function takes evnt structures generated on a different machine
% and ammends the event structures to the local directory
evnt_out = evnt;
for i = 1:length(evnt)
    % Split dpath (check for mac or window format:
    dpath_mac = strsplit(evnt(i).dpath,'/');
    wname_mac = strsplit(evnt(i).wname,'/');
    
    dpath_win = strsplit(evnt(i).dpath,'\');
    wname_win = strsplit(evnt(i).wname,'\');
    if length(dpath_win) > length(dpath_mac)
        dpath_stem_dir = dpath_win((end-1):end);
        wname_stem_dir = wname_win((end-1):end);
    else
        dpath_stem_dir = dpath_mac((end-1):end);
        wname_stem_dir = wname_mac((end-1):end);
    end
    dpath_stem_dir = strjoin(dpath_stem_dir, filesep); % combine end of directory
    dpath_out = [root_dir filesep dpath_stem_dir];
    
    wname_stem_dir = strjoin(wname_stem_dir, filesep);
    wname_out = [root_dir filesep wname_stem_dir];
    
    evnt_out(i).dpath = dpath_out;
    evnt_out(i).wname = wname_out;
end



end
