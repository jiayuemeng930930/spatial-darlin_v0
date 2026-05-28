function call_result = spatial_call_alleles(obj, FQ, CARLIN_def, aligned, depth)

    assert(isscalar(depth) && depth > 0);
    call_result.allele = [];
    call_result.constituents = [];

    if (sum(obj.SEQ_weight) >= depth)
        assert(isa(aligned, 'AlignedSEQDepot'));  
        aligned_SEQ = arrayfun(@(x) aligned.get_alignment_for_SEQ_ind(x), FQ.read_SEQ_valid(obj.SEQ_ind), 'un', false);
        [call_result.allele, call_result.constituents, weight_contribution] ...
                    = CallableCollection.call_alleles_coarse_grain(CARLIN_def, aligned_SEQ, obj.SEQ_weight, true);

        if (~isempty(call_result.allele))
            assert(sum(obj.SEQ_weight(call_result.constituents)) == sum(weight_contribution));
            event = call_result.allele.get_event_structure;

            %if (startsWith(event, 'E') || endsWith(event, 'E'))
            %    call_result.allele = [];
            %    call_result.constituents = [];
            %end
        end
    end

end
