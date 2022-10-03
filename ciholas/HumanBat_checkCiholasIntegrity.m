function [] = HumanBat_checkCiholasIntegrity(ciholas_dir_path, allocated_lines)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   ciholas_dir_path: Path to folder containing *_cdp_*.txt file

ciholas_file = dir(fullfile(ciholas_dir_path,'*cdp*.txt'));
ciholas_filepath = fullfile(ciholas_file.folder, ciholas_file.name);
[s,w] = system(sprintf('wc -l %s', ciholas_filepath));

disp(w);
numLines = strsplit(w, ' ');
numLines = str2num(numLines{1});
assert(numLines == allocated_lines, sprintf('Number of lines in %s does not match expected %d', w, allocated_lines));

end