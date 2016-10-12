function [] = localize_evnt_files(root_dir, evnt_dir)
%% this function loads a directory of event files that were generated
% on a different machine and changes the file directories to include
% the localize directory locations specified in root_dir


%% get list of evnt files:
event_files = dir(evnt_dir); % list of all possible event structures
is_file = true(size(event_files));
for i = 1:length(event_files)
    file_name = event_files(i).name;
    if ~strcmpi(file_name(1),'B');
        is_file(i) = false;
    end
end
event_files(~is_file) = [];


%% Loop through and change:
for k = 1:length(event_files)
    load([evnt_dir filesep event_files(k).name]);
    evnt = dpaths2local(evnt, root_dir); % Localize dpathes and wnames
    save([evnt_dir filesep event_files(k).name], 'evnt');
end

end