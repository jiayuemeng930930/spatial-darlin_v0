classdef ExperimentSummary
    
    properties (SetAccess = immutable)
        CARLIN_def
        alleles
        allele_freqs
    end
    
    methods (Access = public)
        
        function obj = ExperimentSummary(CARLIN_def)
                        
            % Always realign when we get to this stage. If we're merging datasets, 
            % makes sure results are updated, and makes sure that the same
            % sequence always provides the same allele regardless of
            % alignment.
            %
            % NB. Technically don't have to do this anymore since
            % CARLIN_def now stored with summary object, and FromMerge
            % checks that all samples to merge have same CARLIN_def.
            
            obj.CARLIN_def = CARLIN_def;
            
        end
        
    end
    
    methods (Static)
        
        function con = get_legacy_amplicon(s)
            
            assert(isa(s, 'ExperimentSummary'));

            if (~isempty(s.alleles))
                ref_seq = degap(s.alleles{1}.get_ref);
                OC = CARLIN_amplicon(parse_amplicon_file('OriginalCARLIN'));
                TC = CARLIN_amplicon(parse_amplicon_file('TigreCARLIN'));
                if (isequal(ref_seq, OC.seq.CARLIN))
                    con = OC;
                elseif (isequal(ref_seq, TC.seq.CARLIN))                    
                    con = TC;
                else
                    error('Unable to determine CARLIN amplicon used to create legacy file');
                end
            else
                error('Unable to determine CARLIN amplicon used to create legacy file');
            end
        end        
        
        function obj = loadobj(s)
            if (~isempty(s.CARLIN_def))
                obj = s;
                return;
            end
            assert(isa(s, 'ExperimentSummary'), 'loadobj only defined for existing ExperimentSummary objects');
            fprintf('Reformatting legacy ExperimentSummary object\n');
            CARLIN_def = ExperimentSummary.get_legacy_amplicon(s);
            obj = ExperimentSummary(CARLIN_def, s.alleles, s.allele_freqs);
        end
        
        function [obj, sample_map, allele_breakdown_by_sample] = FromMerge(samples)
            
            assert(isa(samples, 'ExperimentSummary'));
            if (size(samples,1)==1)                
                obj = ExperimentSummary(samples.CARLIN_def, samples.alleles, samples.allele_freqs);
                sample_map = {[1:length(obj.alleles)]'};
                allele_breakdown_by_sample = samples.allele_freqs;
                return;
            end
            
            N_samples = length(samples);
            for i = 2:N_samples
                assert(isequal(samples(1).CARLIN_def, samples(i).CARLIN_def), 'Cannot merge ExperimentSummaries from different CARLIN amplicons');
            end
            
            joint_frequency = vertcat(samples.allele_freqs);
            joint_alleles = {samples.alleles}';
            N_alleles = cellfun(@length, joint_alleles);            
            [joint_alleles, ~, sample_map] = ...
                unique_by_freq(cellfun(@(x) degap(x.get_seq), vertcat(joint_alleles{:}), 'un', false), joint_frequency);
            
            obj = ExperimentSummary(samples(1).CARLIN_def, joint_alleles, accumarray(sample_map, joint_frequency));
            
            sample_map = mat2cell(sample_map, N_alleles);        
            
            allele_breakdown_by_sample = zeros(length(joint_alleles), N_samples);
            for i = 1:N_samples
                allele_breakdown_by_sample(sample_map{i},i) = samples(i).allele_freqs;
            end
        end
        
    end
end