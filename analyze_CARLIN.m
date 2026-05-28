function analyze_CARLIN(fastq_file, cfg_type, outdir, varargin)

% For slideSeq, sample script:
% analyze_CARLIN({'C:/Users/ymj5a/Documents/Spatial_tracing/raw_data/TC_right_brain_R1.fastq.gz',...
%    'C:/Users/ymj5a/Documents/Spatial_tracing/raw_data/TC_right_brain_R2.fastq.gz'},...
%    'slideSeqV2','C:/Users/ymj5a/Documents/Spatial_tracing/results_matching/TC_right_brain',...
%    'CARLIN_amplicon','TigreCARLIN',...
%    'ref_CB_file','C:/Users/ymj5a/Documents/Spatial_tracing/barcode_file/right_brain_barcode_matching.txt',...
%    'read_override_CB_denoised', 1, 'read_override_UMI_denoised', 1)
%

    tic;

    % Start logging
    diary([tempname '.txt']);
    diary on;
    fprintf('Writing log to temporary location: %s\n', get(0,'DiaryFile'));
    
    % Parse parameters
    cfg = parse_config_file(cfg_type);
    params = get_parameters(cfg);    
    parse(params, fastq_file, cfg_type, outdir, varargin{:});
    clear fastq_file cfg_type outdir varargin
        
    % Setup directory
    if (~exist(params.Results.outdir, 'dir'))
        mkdir(params.Results.outdir);        
    end

    % Setup CARLIN amplicon
    CARLIN_def = CARLIN_amplicon(parse_amplicon_file(params.Results.CARLIN_amplicon));
    
    % Make FQ representation
    if (strcmp(cfg.type, 'Bulk'))
        FQ = BulkFastQData(params.Results.fastq_file, cfg, CARLIN_def);
    end
    if (strcmp(cfg.type, 'Spatial'))
        ref_CBs = get_Spatial_ref_BCs(params.Results.ref_CB_file);
        FQ = SpatialFastQData(params.Results.fastq_file, cfg, CARLIN_def);                         
    end
    if (strcmp(cfg.type, 'SC'))
        ref_CBs = get_SC_ref_BCs(params.Results.ref_CB_file);
        FQ = SCFastQData(params.Results.fastq_file, cfg, CARLIN_def);                         
    end

    if (isempty(FQ.get_SEQs()))
        fprintf('ERROR: No reads survive filtering. Ensure that your CFG_TYPE and CARLIN_Amplicon settings are correct.\n');
        fprintf('Saving intermediate results...\n');
        try
            save(sprintf('%s/Analysis.mat', params.Results.outdir));
        catch
            save(sprintf('%s/Analysis.mat', params.Results.outdir), '-v7.3', '-nocompression');
        end        
        gzip(sprintf('%s/Analysis.mat', params.Results.outdir));
        delete(sprintf('%s/Analysis.mat', params.Results.outdir));        
        fprintf('Pipeline ABORTED after %g seconds\n', toc);    
        diary off;
        copyfile(get(0,'DiaryFile'), [params.Results.outdir '/Log.txt']);
        return;
    end
    
    % Align unique sequences
    aligned = AlignedSEQDepot(FQ.get_SEQs(), CARLIN_def);
    aligned.sanitize_prefix_postfix();
    aligned.sanitize_conserved_regions(CARLIN_def);
    
    % Make CB collection and call alleles
    if (strcmp(cfg.type, 'Bulk'))
        tag_collection = BulkUMICollection.FromFQ(FQ);
        [tag_collection_denoised, tag_denoise_map] = tag_collection.denoise(CARLIN_def, aligned);
        thresholds = tag_collection_denoised.compute_thresholds(params, FQ);
        tag_called_allele = tag_collection_denoised.call_alleles(CARLIN_def, aligned, thresholds.chosen);
        summary = BulkExperimentReport.create(CARLIN_def, tag_collection_denoised, tag_denoise_map, tag_called_allele, FQ, thresholds);
    end
    if (strcmp(cfg.type, 'SC'))       
        tag_collection = CBCollection.FromFQ(FQ);
        [tag_collection_denoised, tag_denoise_map] = tag_collection.denoise(ref_CBs);
        thresholds = tag_collection_denoised.compute_thresholds(params, FQ, length(ref_CBs));
        tag_called_allele = tag_collection_denoised.call_alleles(CARLIN_def, aligned, [thresholds.CB.chosen, thresholds.UMI.chosen]);
        summary = SCExperimentReport.create(CARLIN_def, tag_collection_denoised, tag_collection, tag_denoise_map, ...
                                            tag_called_allele, FQ, thresholds, ref_CBs);
    end
    if (strcmp(cfg.type, 'Spatial'))   
        matching_results = SpatialCollection.match_CB(FQ, ref_CBs);
        spatial_collection = SpatialCollection.FromFQ(FQ, matching_results, ref_CBs);
        thresholds = spatial_collection.compute_thresholds(params, FQ, length(ref_CBs{1}));
        tag_called_allele = spatial_collection.call_alleles(FQ, CARLIN_def, aligned, [thresholds.CB.chosen, thresholds.UMI.chosen]);
        summary = SpatialExperimentReport.create(CARLIN_def, spatial_collection, tag_called_allele, FQ, thresholds, ref_CBs{1});
    end
   
    % Save just summary values needed for further analysis separately, so
    % they can be opened quickly. Full results saved later, can be big and
    % clunky to reopen for one-off analysis
    
    fprintf('Saving results...');
    
    if (strcmp(cfg.type, 'Bulk'))
        save(sprintf('%s/Summary.mat', params.Results.outdir), 'summary', 'thresholds', 'params');
    end
    if (strcmp(cfg.type, 'SC'))
        save(sprintf('%s/Summary.mat', params.Results.outdir), 'summary', 'thresholds', 'params', 'ref_CBs');
        save(sprintf('%s/Analysis.mat', params.Results.outdir), '-v7.3', '-nocompression');
    end
    if (strcmp(cfg.type, 'Spatial'))
        if (~exist(append(params.Results.outdir,'/Intermediate'), 'dir'))
        mkdir(append(params.Results.outdir,'/Intermediate'));        
        end
        save(sprintf('%s/Summary.mat', params.Results.outdir), 'summary', 'thresholds', 'params', 'ref_CBs');
        save(sprintf('%s/Intermediate/matching_results.mat', params.Results.outdir), 'matching_results');
        save(sprintf('%s/Intermediate/spatial_collection.mat', params.Results.outdir), '-v7.3','-nocompression','spatial_collection');
        save(sprintf('%s/Intermediate/FQ.mat', params.Results.outdir), '-v7.3','-nocompression','FQ');
        save(sprintf('%s/Intermediate/aligned.mat', params.Results.outdir), '-v7.3','-nocompression','aligned');
        save(sprintf('%s/Intermediate/tag_called_allele.mat', params.Results.outdir), '-v7.3','-nocompression','tag_called_allele');
    end
      
    
    % Generate outputs
    
    if (strcmp(cfg.type, 'Bulk'))
        generate_text_output(summary, params, thresholds, params.Results.outdir);
    end
    if (strcmp(cfg.type, 'SC'))
        generate_text_output(summary, params, thresholds, params.Results.outdir, ...
                             tag_collection, tag_collection_denoised, tag_denoise_map, tag_called_allele);
    end
    if (strcmp(cfg.type, 'Spatial'))
        generate_spatial_text_output(summary, params, thresholds, params.Results.outdir, ...
            matching_results, spatial_collection, tag_called_allele);
    end

    %fprintf('Generating allele plot\n');    
    %warning('off', 'MATLAB:hg:AutoSoftwareOpenGL');
    %plot_summary(summary, params.Results.outdir);
    
    %fprintf('Generating diagnostic plot\n');
    %if (strcmp(cfg.type, 'Bulk'))    
    %    suspect_alleles = plot_diagnostic(cfg, FQ, aligned, tag_collection_denoised, tag_denoise_map, tag_called_allele, ...
    %                                      summary, thresholds, params.Results.outdir);
    %else
    %    suspect_alleles = plot_diagnostic(cfg, FQ, aligned, tag_collection_denoised, tag_denoise_map, tag_called_allele, ...
    %                                      summary, thresholds, ref_CBs, params.Results.outdir);
    %end
    %warning('on', 'MATLAB:hg:AutoSoftwareOpenGL');
                                  
    %generate_warnings(summary, params, suspect_alleles, params.Results.outdir);
    
    fprintf('Pipeline COMPLETED in %g seconds\n', toc);
    
    diary off;
    copyfile(get(0,'DiaryFile'), [params.Results.outdir '/Log.txt']);
    
end
