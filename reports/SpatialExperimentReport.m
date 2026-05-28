classdef SpatialExperimentReport < ExperimentReport
    
    methods (Access = public)
        
        function obj = SpatialExperimentReport(CARLIN_def, N, reads, nUMI)
            
            obj = obj@ExperimentReport(CARLIN_def, N, reads, nUMI);
        end
        
    end
    
    methods (Static)
        
        function obj = create(CARLIN_def, CB_collection, CB_called_allele, FQ, thresholds, ref_CBs)
            
            fprintf('Generating report for Spatial experiment\n');

            % output mismatches, positions, called_CB and UMIs
            
            temp = {CB_collection.CBs.UMIs}';
            cb_weights = zeros(size(temp));
            cb_weights(~cellfun(@isempty, temp)) = cellfun(@(x) sum(vertcat(x.SEQ_weight)), temp(~cellfun(@isempty, temp)));
            
            N.reference_tags = length(ref_CBs);
            N.uncleaned_tags = length(FQ.get_CBs());
            N.matched_tags   = length(CB_collection.CBs);            
            N.common_tags    = sum(cb_weights >= thresholds.CB.chosen);
            N.common_tags_3  = sum(cb_weights >= 3);  % for test; will be removed
            N.common_tags_10  = sum(cb_weights >= 10);
            N.called_tags    = size(CB_called_allele,1);
            
            fprintf('Finish generating report for N\n');
            
            reads.in_fastq            = FQ.Nreads;
            for fn = fieldnames(FQ.masks)'
                reads.(fn{1}) = length(FQ.masks.(fn{1}));
            end
            [is, where] = ismember(vertcat(CB_collection.CBs(:).called_CB), FQ.CB);
            assert(all(is));
            reads.matched_tags         = sum(ismember(FQ.read_CB(FQ.masks.valid_lines), where));
            reads.common_tags          = sum(cb_weights(cb_weights >= thresholds.CB.chosen));
            
            [is, where] = ismember({CB_called_allele.CB}',{CB_collection.CBs.CB}');
            assert(all(is));
            reads.called_tags_total  = cb_weights(where);
            reads.called_tags_allele = zeros(length(where),1);

            temp = cell(length(where),1);
	    parfor i = 1:length(where)
                %reads.called_tags_allele(i) = sum(arrayfun(@(j) sum(CB_collection.CBs(where(i)).UMIs(j).SEQ_weight ...
                %                                  (CB_called_allele(i).umi_call_result(j).constituents)), vertcat(CB_called_allele(i).constituents{:})));
            
                temp{i} = sum(arrayfun(@(j) sum(CB_collection.CBs(where(i)).UMIs(j).SEQ_weight ...
                                                  (CB_called_allele(i).umi_call_result(j).constituents)), vertcat(CB_called_allele(i).constituents{:})));
	    end

            reads.called_tags_allele = cell2mat(temp);

            %for i = 1:length(where)
            %    template_index = strcmp(CB_called_allele(i).allele, CARLIN_def.seq.CARLIN);
            %    temp = CB_called_allele(i).constituents(~template_index);
            %    reads.eventful_tags_allele(i) = sum(arrayfun(@(j) sum(CB_collection.CBs(where(i)).UMIs(j).SEQ_weight ...
            %                                      (CB_called_allele(i).umi_call_result(j).constituents)), vertcat(temp{:})));  
            %end

            reads.called_tags_total = sum(reads.called_tags_total);
            reads.called_tags_allele = sum(reads.called_tags_allele);
            %freads.eventful_tags_allele = sum(reads.eventful_tags_allele);
            
            fprintf('Finish generating report for reads\n');
            % generate UMI data
            nUMI.nUMI_uncleaned_unique  = length(unique(FQ.read_UMI(FQ.masks.valid_lines)));
            nUMI.nUMI_matched           = sum(cellfun(@length, {CB_collection.CBs.UMIs}));
            nUMI.nUMI_common            = sum(cellfun(@(x) sum(cellfun(@sum, {x.SEQ_weight})...
                                           >=thresholds.UMI.chosen), {CB_collection.CBs.UMIs}));
            nUMI.nUMI_called            = sum(cellfun(@(x) size(vertcat(x.allele),1), {CB_called_allele.umi_call_result}));
            
            fprintf('Finish generating report for nUMI\n');
            obj = SpatialExperimentReport(CARLIN_def, N, reads,nUMI);
            
        end        
        
        function obj = loadobj(s)
            
            if (~isempty(s.CARLIN_def))
                obj = s;
                return;
            end
            
            assert(isa(s, 'SpatialExperimentReport'), 'loadobj only defined for existing SCExperimentReport objects');
            fprintf('Reformatting legacy SpatialExperimentReport object\n');
            CARLIN_def = ExperimentSummary.get_legacy_amplicon(s);            
            obj = SpatialExperimentReport(CARLIN_def, s.alleles, s.allele_colony, s.N, s.reads);
            
        end
        
    end
    
end
