classdef (Sealed=true) SpatialCollection < TaggedCollection & CallableCollection
    
    properties(SetAccess = immutable)
        CBs
    end
    
    methods (Static)
        
        function out = FromFQ(FQ, matching_results, ref_CBs)
            assert(isa(FQ, 'SpatialFastQData'));
            called_CBs = matching_results(:,1);
            mapped_CBs = vertcat(matching_results{:,2});
            num_mismatches = matching_results(:,3);
            ref_list = ref_CBs{1};
            x_coord = ref_CBs{2};
            y_coord = ref_CBs{3};
            unique_CBs = unique(mapped_CBs);
            unique_CBs(ismember(unique_CBs,'Duplicate Match')) = [];
            unique_CBs(ismember(unique_CBs,'No Match')) = [];
            N = length(unique_CBs);
            Spatial_data = cell(N,1);

            fprintf('Building CB collection from FastQ with %d unique CBs\n', N);

            for i = 1:N
                Spatial_data{i} = TaggedCollection.MergeFromFQ(FQ, unique_CBs{i}, mapped_CBs, ...
                    called_CBs, num_mismatches, ref_list, x_coord, y_coord);
            end
            out = SpatialCollection(vertcat(Spatial_data{:}));
        end


        function matching_results = match_CB(obj, ref_CBs)

            CBs = obj.get_CBs();
            N = size(CBs,1);
            fprintf('Matching %d CBs to reference\n', N);
            ref_list = ref_CBs{1}; 
            matching_results = cell(N,3); 
            
            matching_CB = cell(N,1);
	    mismatch = cell(N,1);
	    input_CB = cell(N,1);
	    	
            parfor i=1:N
                [matching_CB{i},mismatch{i}] = fuzzy_match(CBs{i}, ...
                    ref_list, 3);
                input_CB{i} = CBs{i};
            end 
            
            matching_results = table(input_CB, matching_CB, mismatch);
            matching_results = table2cell(matching_results);
	    
            index = cellfun(@(x) isnumeric(x), matching_results(:,3), 'UniformOutput', true);
            temp = cell2mat(matching_results(index,3));

            mismatch0 = sum(temp(:) == 0);
            mismatch1 = sum(temp(:) == 1);
            mismatch2 = sum(temp(:) == 2);
            mismatch3 = sum(temp(:) == 3);

            fprintf('From %d called CBs, (%d,%d,%d,%d) matched to reference with 0,1,2,3 mismatches\n', ...
                N, mismatch0, mismatch1, mismatch2, mismatch3);
            
        end        
   
    end

        %[SEQ_ind, SEQ_weight] = get_SEQ_ind_by_UMI(obj, UMI, filter)

    methods (Access = public)

        function obj = SpatialCollection(CB_data)
            obj.CBs = CB_data;
        end
        
        function call_result = call_alleles(obj, FQ, CARLIN_def, aligned, depth)
            
            fprintf('Calling alleles for CB collection\n');
            assert(isequal(size(depth), [1, 2]) && all(depth > 0) && (depth(1) >= depth(2)));                        
            assert(isa(aligned, 'AlignedSEQDepot'));
            
            N = size(obj.CBs,1);
            call_result = cell(N,1);
            temp = obj.CBs;
            
            % This can be parfor on a bigmem machine
            for i = 1:N
                call_result{i} = temp(i).call_multiple_alleles(FQ, CARLIN_def, aligned, depth);
            end
            
            % Call result has a per-CB allele (which may or may not be empty) and a 
            % per-UMI call result (which may or may not be
            % empty, even when per-CB allele is empty).
            
            % Get cell array elements where all fields are empty.
            uncalled = cellfun(@isempty, call_result);
            
            % This concatenation collapses these all-empty elements
            call_result = vertcat(call_result{:});
            
            if (isempty(call_result))
                temp = {'CB', 'allele'};
                temp{2,1} = {};
                call_result = struct(temp{:});
            else
                % Get mapped CBs of uncollapsed elements and tack them on as a field
                temp = {obj.CBs(~uncalled).CB};
                [call_result.CB] = temp{:};

                % Get called CBs
                temp = {obj.CBs(~uncalled).called_CB};
                [call_result.called_CB] = temp{:};

                % Get x and y coordinates
                temp = vertcat(obj.CBs(~uncalled).position);
                [call_result.x_coord] = temp{:,1};
                [call_result.y_coord] = temp{:,2};

                % get minimum number of mismatches
                temp = {obj.CBs(~uncalled).num_mismatches};
                temp = cellfun(@(x) min(x), temp, 'un', true);
                temp = num2cell(temp);
                [call_result.min_mismatches] = temp{:};

                % Remove CBs where there is no callable per-CB allele
                call_result(cellfun(@isempty, {call_result.allele}')) = [];

                % Corner case, if nothing is callable, make me a 0x1 array
                % instead of a 1x0 array
                if (isempty(call_result))
                    call_result = call_result';
                end
            end
            fprintf('...%d CBs callable\n', size(call_result,1));

            temp = vertcat(call_result.min_mismatches);
            fprintf('...(%d,%d,%d,%d) callable CBs matched to reference with (0,1,2,3) mismatches \n', ...
                sum(temp(:) == 0), sum(temp(:) == 1), sum(temp(:) == 2), sum(temp(:) == 3));

        end
        
        function [out, collapse_map] = denoise(obj)
            
            fprintf('Denoising CB collection\n');
            
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
