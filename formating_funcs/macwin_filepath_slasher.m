function dpath_out = macwin_filepath_slasher(dpath_in)
% converts file directories in mac or windows
% format to file directories in the accepted local format

split_mac = strsplit(dpath_in, '/');
split_win = strsplit(dpath_in, '\');

if length(split_mac) > length(split_win)
    dpath_out = strjoin(split_mac, filesep);
else
    dpath_out = strjoin(split_win, filesep);
end


end