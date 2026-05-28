% fuzzy match a single string to a list of strings using Hamming Score
function [matching_CB, num_mismatches] = fuzzy_match(input_string, ref_list, cutoff)

    hamming_score = cellfun(@(x) sum(input_string ~= x), ref_list);
    min_dist = min(hamming_score);

    if (le(min_dist,cutoff))
        matching_CB = {ref_list(hamming_score == min_dist)};
        num_mismatches = cellfun(@(x) min_dist-count(x(input_string ~= x),'N'), matching_CB{1,1});
        num_match = min(num_mismatches);
        if (sum(num_mismatches(:) == num_match) == 1)
            matching_CB = matching_CB{1,1}(num_mismatches == num_match);
            num_mismatches = num_match;
        else
            matching_CB = 'Duplicate Match';
            num_mismatches = 'NaN';
        end
    else
        matching_CB = 'No Match';
        num_mismatches = 'NaN';
    end
   
end