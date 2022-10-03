% Stop-gap solution for if you need to re-cluster ciholas flights

% 1. After making the new ciholas_flight_struct in HumanBat_Align2Ephys,
% 2. relabel it as "ciholas_flight_struct_new".
% 3. Load in the old ciholas_flight_struct.

% Make vectors of the fclus labels
new_struct_cluster_labels=[]; for i=1:length(ciholas_flight_struct); new_struct_cluster_labels(i) = ciholas_flight_struct_new{i}.fclus; end
old_struct_cluster_labels=[]; for i=1:length(ciholas_flight_struct); old_struct_cluster_labels(i) = ciholas_flight_struct{i}.fclus; end
cluster_labels_to_change = find(new_struct_cluster_labels~=old_struct_cluster_labels);

% Check that the start times are the same!
old_struct_time_labels=[]; for i=1:length(ciholas_flight_struct); old_struct_time_labels(i) = ciholas_flight_struct{i}.fstart; end
new_struct_time_labels=[]; for i=1:length(ciholas_flight_struct_new); new_struct_time_labels(i) = ciholas_flight_struct_new{i}.fstart; end
if sum(old_struct_time_labels == new_struct_time_labels) == length(ciholas_flight_struct)
    disp("all good!");
else
    disp("PROBLEM");
end

%% Replace the old labels with the new labels
for i=1:length(cluster_labels_to_change)
    ciholas_flight_struct{cluster_labels_to_change(i)}.fclus = ciholas_flight_struct_new{cluster_labels_to_change(i)}.fclus;
end

% Resave the struct so it is correct
save(strcat(exp_data_path,'ciholas/ciholas_bat_final_flight_structure.mat'),'ciholas_flight_struct','-v7.3');

% Repeat this process for the trimmed, pruned, and tripodded matfiles,
% bearing in mind that the ordering of the flights is different
other_matfiles = dir(fullfile(strcat(exp_data_path,'ciholas/*flight_structure*')));
for i=1:length(other_matfiles)
    if strcmp(other_matfiles(i).name,'ciholas_bat_final_flight_structure.mat')
        continue
    else
        to_alter = load(strcat(exp_data_path,'ciholas/',other_matfiles(i).name));
        disp(strcat("Loading ",other_matfiles(i).name));

        if isfield(to_alter,'ciholas_flight_struct')
            for j=1:length(to_alter.ciholas_flight_struct)
                for k=1:length(ciholas_flight_struct)
                    if ~isempty(to_alter.ciholas_flight_struct{j})
                        if to_alter.ciholas_flight_struct{j}.fstart == ciholas_flight_struct{k}.fstart
                            to_alter.ciholas_flight_struct{j}.fclus = ciholas_flight_struct{k}.fclus;
                        end
                    end
                end
            end
            ciholas_flight_struct_old = ciholas_flight_struct;
            ciholas_flight_struct = to_alter.ciholas_flight_struct;
            save(strcat(exp_data_path,'ciholas/',other_matfiles(i).name),'ciholas_flight_struct','-v7.3');

        elseif isfield(to_alter,'ciholas_flight_struct_resort')
            for j=1:length(to_alter.ciholas_flight_struct_resort)
                for k=1:length(ciholas_flight_struct)
                    if to_alter.ciholas_flight_struct_resort{j}.fstart == ciholas_flight_struct_old{k}.fstart
                        to_alter.ciholas_flight_struct_resort{j}.fclus = ciholas_flight_struct_old{k}.fclus;
                    end
                end
            end
            ciholas_flight_struct_resort = to_alter.ciholas_flight_struct_resort;

            % Resort according to cluster now.
            cv_k = [];
            for k=1:length(ciholas_flight_struct_resort)
                cv_k(k) = ciholas_flight_struct_resort{k}.fclus;
            end
            [im,imk] = sort(cv_k); ciholas_flight_struct_resort_resorted = {};
            for kk=1:length(imk)
                ciholas_flight_struct_resort_resorted{kk} = ciholas_flight_struct_resort{imk(kk)};
            end

            clear ciholas_flight_struct_resort
            ciholas_flight_struct_resort = ciholas_flight_struct_resort_resorted;

            save(strcat(exp_data_path,'ciholas/',other_matfiles(i).name),'ciholas_flight_struct_resort','-v7.3');
        end
    end
end
        

