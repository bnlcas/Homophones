file_names_tmp = dir(behavior_dat_dir);
file_names_tmp(1:2) = [];

for i = 1:length(file_names_tmp)
    file_names{i} = file_names_tmp(i).name;
end
file_names(strcmpi(file_names,'')) = [];
block_tags = cell(1,length(file_names));
for i = 1:length(file_names)
    filename_split = strsplit(file_names{i},'_');
    if strcmpi(filename_split{1}, 'Homophone')
        if strcmpi(filename_split{2}, 'Exp')
            end_tag = filename_split{end};
            if ~strcmpi(end_tag(1),'B')
                block_tags{i} = filename_split{end-1};
            end
            
        end
    else    
        block_tags{i} = '';
    end   
end
unique_blocks = unique(block_tags); % Find unique blocks...
unique_blocks(strcmpi(unique_blocks,'')) = [];

