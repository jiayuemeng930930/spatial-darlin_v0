classdef SpatialData < TaggedCollection & CallableCollection
 
    properties(SetAccess = immutable)
        CB
        UMIs
        called_CB
        num_mismatches
        position
    end        

    methods (Access = public)
        
        function obj = SpatialData(CB, called_CB, mismatches, UMI_data, position)
            obj.CB = CB;
            obj.called_CB = called_CB;
            obj.num_mismatches = mismatches;
            obj.position = position;
            obj.UMIs = UMI_data;
        end

        function call_result = call_multiple_alleles(obj, FQ, CARLIN_def, aligned, depth)
          assert(isequal(size(depth), [1, 2]) && all(depth > 0) && (depth(1) >= depth(2)));
          call_result = []; 
          if (~isempty(obj.UMIs))
                umi_weights = vertcat(obj.UMIs.SEQ_weight);
                if (sum(umi_weights) >= depth(1))
                    assert(isa(aligned, 'AlignedSEQDepot'));
                    call_result.umi_call_result = cell(size(obj.UMIs,1),1);
                    for i = 1:size(obj.UMIs,1)
                        call_result.umi_call_result{i} = ...
                            CallableCollection.spatial_call_alleles(obj.UMIs(i), FQ, CARLIN_def, aligned, depth(2));
                    end
                    call_result.umi_call_result = vertcat(call_result.umi_call_result{:});
                    [call_result.umi_call_result.UMI] = obj.UMIs.UMI;
                    [call_result.allele, call_result.constituents, call_result.read_number, call_result.allele_contribution] = ...
                        CallableCollection.call_alleles_multiple_grain(CARLIN_def, {call_result.umi_call_result(:).allele}', ...
                        {call_result.umi_call_result(:).constituents}');
                end
            end
        end

       
        function call_result = call_alleles(obj, FQ, CARLIN_def, aligned, depth)
            
            assert(isequal(size(depth), [1, 2]) && all(depth > 0) && (depth(1) >= depth(2)));
            call_result = [];
            
            if (~isempty(obj.UMIs))
                umi_weights = vertcat(obj.UMIs.SEQ_weight);
                if (sum(umi_weights) >= depth(1))
                    assert(isa(aligned, 'AlignedSEQDepot'));
                    call_result.umi_call_result = cell(size(obj.UMIs,1),1);
                    for i = 1:size(obj.UMIs,1)
                        call_result.umi_call_result{i} = ...
                            CallableCollection.spatial_call_alleles(obj.UMIs(i), FQ, CARLIN_def, aligned, depth(2));
                    end
                    call_result.umi_call_result = vertcat(call_result.umi_call_result{:});
                    [call_result.umi_call_result.UMI] = obj.UMIs.UMI;
                    [call_result.allele, call_result.constituents] = ...
                        CallableCollection.call_alleles_coarse_grain(CARLIN_def, {call_result.umi_call_result(:).allele}');
                end
            end
        end

          function [out, collapse_map] = denoise(obj, ref_list)
            
            fprintf('Denoising CB collection\n');

            [is, which_ref] = ismember({obj.CBs.CB}', ref_list);
            CB_weight = cellfun(@(x) sum(vertcat(x.SEQ_weight)), {obj.CBs.UMIs}');
            CB_MiSEQ_internal = TaggedCollection.directional_adjacency_top_down_denoiser({obj.CBs.CB}', CB_weight, is);
            
            self_same = CB_MiSEQ_internal==[1:length(CB_MiSEQ_internal)]';

            matched = ismember(CB_MiSEQ_internal,find(self_same&is));
            CB_cleaned = cell(size(CB_weight));
            CB_cleaned(matched) = ref_list(which_ref(CB_MiSEQ_internal(matched)));

            matched = find(matched);
            k = {obj.CBs(matched).CB};
            v = CB_cleaned(matched);
            collapse_map.CB = containers.Map(k, v);
            [v, ~, ind] = unique(v);
            N = size(v,1);
            
            fprintf('...(%d/%d) MiSEQ CBs matched to (%d/%d) reference CBs\n', ...
                length(matched), length(obj.CBs), N, length(ref_list));
            
            fprintf('Denoising UMIs\n')
            CB_data = cell(N,1);
            UMI_map = cell(N,1); 
            merged_data = arrayfun(@(i) vertcat(obj.CBs(matched(ind==i))), 1:N, 'un', false);            
            
            parfor i = 1:N
                temp = CBData.FromMerge(v{i}, merged_data{i});
                [CB_data{i}, UMI_map{i}] = temp.denoise();
            end
            collapse_map.UMI = containers.Map(v, UMI_map);
            out = CBCollection(vertcat(CB_data{:}));
            
        end
        
        function CB_freq = get_CB_freq(obj)
            CB_freq = sort(cellfun(@(x) sum(vertcat(x.SEQ_weight)), {obj.CBs.UMIs}'), 'descend');
        end
        
        function UMI_freq = get_UMI_freq(obj)
            UMI_freq = cellfun(@(x) cellfun(@(y) sum(y), {x.SEQ_weight}'), {obj.CBs.UMIs}', 'un', false);
            UMI_freq = sort(vertcat(UMI_freq{:}), 'descend');
        end
        
        function thresholds = compute_thresholds(obj, params, FQ, N_ref_CBs)
            
            fprintf('Computing thresholds for CB collection\n');
            p = cellfun(@(x) 10.^((double(x)-33)/-10), FQ.QC(FQ.masks.valid_lines), 'un', false);
            p = sort(horzcat(p{:}), 'descend');
            p = p(ceil(length(p)/20));
            CB_freq = obj.get_CB_freq();
            thresholds.CB = TaggedCollection.threshold_function(CB_freq, min(params.Results.max_cells, N_ref_CBs), ...
                                                                length(FQ.masks.valid_lines), p, ...
                                                                params.Results.read_cutoff_CB_denoised, ...
                                                                params.Results.read_override_CB_denoised, 'CB');
            fprintf('...%d common CBs\n', sum(CB_freq >= thresholds.CB.chosen));
            
            UMI_freq = obj.get_UMI_freq();
            thresholds.UMI = TaggedCollection.threshold_function(UMI_freq, inf, ...
                                                                 length(FQ.masks.valid_lines), p, ...
                                                                 params.Results.read_cutoff_UMI_denoised, ...
                                                                 params.Results.read_override_UMI_denoised, 'UMI');
            fprintf('...%d common UMIs\n', sum(UMI_freq >= thresholds.UMI.chosen));
                                           
        end
    end
end