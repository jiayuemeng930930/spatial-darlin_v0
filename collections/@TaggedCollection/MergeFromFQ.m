function out = MergeFromFQ(FQ, CB, mapped_CBs, called_CBs, num_mismatches, ref_list, x_coord, y_coord)

    idx = find(strcmp(mapped_CBs, CB));
    called_CB = called_CBs(idx);
    [UMIs, filter] = get_UMIs_by_CB(FQ, called_CB);
    uniq_UMIs = unique(vertcat(UMIs{:}));
    UMI = cell(uniq_UMIs);
    filter = vertcat(filter{:});
    filter = intersect(filter, FQ.masks.valid_lines);
    [is, which_UMI] = ismember(UMI, FQ.UMI); 
    if (any(is))
        [SEQ_ind, ~, SEQ_weight] = unique([FQ.read_UMI(filter) filter], 'rows');
        SEQ_weight = accumarray(SEQ_weight,1);
        [SEQ_ind, SEQ_weight] = ...
            arrayfun(@(i) deal(SEQ_ind(SEQ_ind(:,1)==i,2), SEQ_weight(SEQ_ind(:,1)==i)), which_UMI(is), 'un', false);
        SEQ_ind = cellfun(@(x,y) sortrows([x, y], 2, 'descend'), SEQ_ind, SEQ_weight, 'un', false);
        [SEQ_ind, SEQ_weight] = cellfun(@(x) deal(x(:,1), x(:,2)), SEQ_ind, 'un', false);
        SEQ_ind(is) = SEQ_ind;
        SEQ_weight(is) = SEQ_weight;
    end
    SEQ_ind(~is) = {[]};
    SEQ_weight(~is) = {[]};
    UMI_data = cellfun(@(u,i,w) UMIData(u,i,w), UMI, SEQ_ind, SEQ_weight, 'un', false);

    % get other information
    mismatches = vertcat(num_mismatches{idx});
    index = find(strcmp(ref_list, CB));
    position = {x_coord{index}, y_coord{index}};
    
    % output to a struct
    out = SpatialData(CB, called_CB, mismatches, vertcat(UMI_data{:}), position);

end

