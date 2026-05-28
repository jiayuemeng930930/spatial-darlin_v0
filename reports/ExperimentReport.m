classdef (Abstract=true) ExperimentReport < ExperimentSummary
    
    properties (SetAccess = immutable)
        reads
        N
        nUMI
        CB_called_allele
    end
    
    methods (Access = public)
        
        function obj = ExperimentReport(CARLIN_def, N, reads, nUMI)
            obj = obj@ExperimentSummary(CARLIN_def);
            obj.N = N;
            obj.reads = reads;
            obj.nUMI = nUMI;
        end
        
        function conv = ExperimentSummary(obj)
            conv = ExperimentSummary(obj.alleles, obj.allele_freqs, true);
        end
        
    end
    
end