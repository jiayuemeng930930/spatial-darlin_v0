function results = subsetCARLIN(dir, number_select, num_mismatches, min_reads, num_perturbation,outdir) 
    % Load Data
    fprintf('Starting loading samples\n');
    load(append(dir,'Intermediate/FQ'));
    fprintf('Finish loading FQ\n');
    load(append(dir,'Intermediate/matching_results'));
    fprintf('Finish loading matching results\n');
    load(append(dir,'Intermediate/tag_called_allele'));
    fprintf('Finish loading tag_called_alleles\n');
    load(append(dir,'Intermediate/spatial_collection'));
    fprintf('Finish loading spatial_collection\n');
    results = cell(num_perturbation,1);
    for i=1:num_perturbation
        
        for N = 1:length(number_select)
            index = randperm(FQ.Nreads, number_select(N))';
            index = intersect(index, FQ.masks.valid_lines);
            CB_index = unique(FQ.read_CB(index));
            CBs = FQ.CB(CB_index);
            [~,matched_CB_index] = ismember(CBs,matching_results(:,1));
            matched_CB = matching_results(matched_CB_index,:);
            matched_CB(strcmp(matched_CB(:, 2), 'No Match'), :) = [];
            matched_CB(strcmp(matched_CB(:, 2), 'Duplicate Match'), :) = [];
            [~, tag_idx] = ismember(vertcat(matched_CB{:,2}), {tag_called_allele.CB}');
            tag_idx(tag_idx == 0) = [];
            tag_idx = unique(tag_idx);
            matched_tag = tag_called_allele(tag_idx,:);
            matched_tag_filered = matched_tag([matched_tag.min_mismatches]' <= num_mismatches,:);
            [~, spatial_idx] = ismember(vertcat(matched_tag_filered(:).CB),{spatial_collection.CBs.CB}');
            matched_tag_filtered = {spatial_collection.CBs(spatial_idx).CB}';
                CB_number = 0;
                for a = 1:size(spatial_idx,1)
                    line_number = vertcat(spatial_collection.CBs(spatial_idx(a)).UMIs.SEQ_ind);
                    line_number = intersect(line_number, index);
	                if size(line_number,1) >= min_reads	
                        CB_number = CB_number + 1;
	                end
                end

            results{i}.number_CBs(N) = CB_number;
            
        end
    fprintf(append('Perturbation ', string(i), ' finished\n'));
    end

    if (~exist(outdir, 'dir'))
        mkdir(outdir);        
    end
    fid = fopen(append(outdir,'/number_CB_called_by_perterbation.txt'), 'wt'); 
    s = join(repmat("%d\t ",1,length(number_select)-1),'');
    s = join([s, "%d\n"],'');
    t = join(repmat("%s\t ",1,length(number_select)-1),'');
    t = join([t, "%s\n"],'');
    fprintf(fid, t, strcat('Num_reads_', string(number_select)));
    for i = 1:length(results)
        fprintf(fid, s, results{i}.number_CBs);
    end
    fclose(fid);    
end
