function [matched_CB, num_mismatches] = match_CB(input_CB, ref_list)
    
    matched_CB = ref_list(strcmp(input_CB, ref_list));
    if (isempty(matched_CB))
        [matched_CB, num_mismatches] = fuzzy_match(input_CB, ref_list, 3);
    else
        num_mismatches = 0;
    end
   
end