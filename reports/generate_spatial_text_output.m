function generate_spatial_text_output(summary, params, thresholds, output_path, matching_results, spatial_collection, tag_called_allele, varargin)

    fprintf('Generating results file\n');
    
    assert(isa(summary, 'ExperimentReport'));
    
    outdir = params.Results.outdir;
    filename = append(outdir,'/Results.txt'); 
    fid = fopen(filename, 'wt');    
    
    fprintf(fid, 'INPUT\n\n');
    
    fprintf(fid, '%-55s %s\n', 'CARLIN Amplicon:', params.Results.CARLIN_amplicon);
    fprintf(fid, '%-55s %s\n', 'CfgType:', params.Results.cfg_type);    
    if (iscell(params.Results.fastq_file))
        fastq_files = {params.Results.fastq_file{:}};
        fprintf(fid, '%-55s %s\n', 'Input FastQ File(s):', fastq_files{1});
        for i = 2:length(fastq_files)
            fprintf(fid, '%-55s %s\n', '', fastq_files{i});
        end
    else
        fprintf(fid, '%-55s %s\n', 'Input FastQ File(s):', params.Results.fastq_file);
    end
    fprintf(fid, '%-55s %s\n', 'Reference barcode file:', params.Results.ref_CB_file);
    fprintf(fid, '%-55s %s\n', 'Output directory:', params.Results.outdir);
    
    fprintf(fid, '\n%-55s %+10s %+6s\n\n', 'READ BREAKDOWN', 'Reads', '%');
    
    for x = fieldnames(summary.reads)'
        fprintf(fid, '%-55s %10d %6.0f\n', [ x{1} ':'], ...
            summary.reads.(x{1}), max(round(double(summary.reads.(x{1}))/double(summary.reads.in_fastq)*100),0));
    end
    
    fprintf(fid, '\nTHRESHOLDS\n\n');
    
    for x = fieldnames(thresholds.CB)'
        fprintf(fid, '%-55s %10d\n', ['CB_' x{1} ':'], thresholds.CB.(x{1}));
    end
    fprintf(fid, '\n');
    for x = fieldnames(thresholds.UMI)'
        fprintf(fid, '%-55s %10d\n', ['UMI_' x{1} ':'], thresholds.UMI.(x{1}));
    end
 
    fprintf(fid, '\nCELL BARCODE BREAKDOWN\n\n');    
    for x = fieldnames(summary.N)'
        fprintf(fid, '%-55s %10d\n', [erase(x{1}, '_tags') ':'], summary.N.(x{1}));
    end

    fprintf(fid, '\nUMI BREAKDOWN\n\n');
    for x = fieldnames(summary.nUMI)'
        fprintf(fid, '%-55s %10d\n', [x{1} ':'], summary.nUMI.(x{1}));
    end
          
    fprintf(fid, '\nMAPPING MISMATCHES of CALLED/UNIQUE CBs\n\n'); 
    
    called_CB = vertcat(tag_called_allele.called_CB);
    temp = matching_results(cellfun(@(x) ismember(x, called_CB),matching_results(:,1), 'UniformOutput', true),:);
    N = length(temp);
    fprintf(fid, '%-55s %10d\n', 'Exact Match:', sum(vertcat(temp{:,3}) == 0));
    fprintf(fid, '%-55s %10d\n', '1 Mismatch:', sum(vertcat(temp{:,3}) == 1));
    fprintf(fid, '%-55s %10d\n', '2 Mismatches:', sum(vertcat(temp{:,3}) == 2));
    fprintf(fid, '%-55s %10d\n', '3 Mismatches:', sum(vertcat(temp{:,3}) == 3));
    fprintf(fid, '%-55s %10d\n', 'Total Called/Unique CBs:', N);

    fprintf(fid, '\nPosition of CALLED/UNIQUE CB MAPPING MISMATCHES\n\n'); 
    % Get the position of mismatch 
    position = cell(N,1);
    nucleotide = cell(N,1);
    for i = 1:N
        input = temp{i,2}{:};
        position{i} = find(temp{i,1} ~= temp{i,2}{:});
        nucleotide{i} = input(position{i});
    end
    position = horzcat(position{:})';
    num_mismatch = length(position)-sum(horzcat(nucleotide{:}) =='N');
    assert(num_mismatch == sum(vertcat(temp{:,3})));

    % Count number of mismatches at each position (1-14)
    N1  = sum(position(:) == 1);
    N2  = sum(position(:) == 2);
    N3  = sum(position(:) == 3);
    N4  = sum(position(:) == 4);
    N5  = sum(position(:) == 5);
    N6  = sum(position(:) == 6);
    N7  = sum(position(:) == 7);
    N8  = sum(position(:) == 8);
    N9  = sum(position(:) == 9);
    N10 = sum(position(:) == 10);
    N11 = sum(position(:) == 11);
    N12 = sum(position(:) == 12);
    N13 = sum(position(:) == 13);
    N14 = sum(position(:) == 14);

    fprintf(fid, '%-55s %10d %6.0f\n', '1st Nucleotide:', N1, round(N1/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '2nd Nucleotide:', N2, round(N2/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '3rd Nucleotide:', N3, round(N3/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '4th Nucleotide:', N4, round(N4/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '5th Nucleotide:', N5, round(N5/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '6th Nucleotide:', N6, round(N6/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '7th Nucleotide:', N7, round(N7/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '8th Nucleotide:', N8, round(N8/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '9th Nucleotide:', N9, round(N9/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '10th Nucleotide:', N10, round(N10/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '11th Nucleotide:', N11, round(N11/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '12th Nucleotide:', N12, round(N12/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '13th Nucleotide:', N13, round(N13/num_mismatch*100));
    fprintf(fid, '%-55s %10d %6.0f\n', '14th Nucleotide:', N14, round(N14/num_mismatch*100));

    fprintf(fid, '\nMAPPING MISMATCHES of CALLED/DENOISED CBs\n\n'); 

    temp = vertcat(tag_called_allele.min_mismatches);
    fprintf(fid, '%-55s %10d\n', 'Exact Match:', sum(temp(:) == 0));
    fprintf(fid, '%-55s %10d\n', '1 Mismatch:', sum(temp(:) == 1));
    fprintf(fid, '%-55s %10d\n', '2 Mismatches:', sum(temp(:) == 2));
    fprintf(fid, '%-55s %10d\n', '3 Mismatches:', sum(temp(:) == 3));
    fprintf(fid, '%-55s %10d\n', 'Total Called/Denoised CBs:', length(temp));

    fclose(fid);

    %% Output CB mapped results

    fprintf('Generating list of called and mapped CBs\n');

    fid = fopen(append(outdir,'/CB_Mapping.txt'), 'wt'); 
    fprintf(fid, 'called_CB\tMapped_CB\tNum_mismatches\n');
    for i = 1:length(matching_results)
        if isnumeric(matching_results{i,3})
            fprintf(fid, '%s\t %s\t %d\n', matching_results{i,1},matching_results{i,2}{:}, matching_results{i,3});
        else
            fprintf(fid, '%s\t %s\t %s\n', matching_results{i,1},matching_results{i,2}, matching_results{i,3});
        end
    end
    fclose(fid);

    
    %% output UMIs and number of reads file

    fprintf('Generating number of reads for CBs and UMIs\n');
    
    idx_called = find(ismember({spatial_collection.CBs.CB}',{tag_called_allele.CB}'));
    CB = cellfun(@(x) {spatial_collection.CBs(x).CB}, num2cell(idx_called), 'UniformOutput', false);
    CB = vertcat(CB{:});
    weight_CB = cellfun(@(x) sum(vertcat(spatial_collection.CBs(x).UMIs.SEQ_weight)), num2cell(idx_called), 'UniformOutput', false);
    weight_CB = vertcat(weight_CB{:});
    UMI = cellfun(@(x) {spatial_collection.CBs(x).UMIs.UMI}, num2cell(idx_called), 'UniformOutput', false);
    num_UMI = cellfun(@(x) length(x), UMI);

    weight_UMI = cellfun(@(x) {spatial_collection.CBs(x).UMIs.SEQ_weight}, num2cell(idx_called), 'UniformOutput', false);
    UMI = horzcat(UMI{:})';
    weight_UMI = horzcat(weight_UMI{:})';
    weight_UMI = cellfun(@(x) length(x), weight_UMI);

    % output number of UMIs and number of reads associated with each CB
    assert(length(CB) == length(weight_CB));
    assert(length(CB) == length(num_UMI));
    fid = fopen(append(outdir,'/CB_nUMI_nReads.txt'), 'wt'); 
    fprintf(fid, 'Mapped_CB\tnUMIs\tnReads\n');
    for i = 1:length(CB)
        fprintf(fid, '%s\t %d\t %d\n', CB{i}, num_UMI(i), weight_CB(i));
    end
    fclose(fid);

    % output number of reads associated with each UMI
    assert(length(UMI) == length(weight_UMI));
    assert(length(CB) == length(num_UMI));
    fid = fopen(append(outdir,'/UMI_nReads.txt'), 'wt'); 
    fprintf(fid, 'UMIs\tnReads\n');
    for i = 1:length(CB)
        fprintf(fid, '%s\t %d\n', UMI{i}, weight_UMI(i));
    end
    fclose(fid);
    

    %%% output AlleleAnnotations.txt
    fprintf('Generating allele files\n');
    fid = fopen(append(outdir,'/Alleles.txt'), 'wt');
    N = length(tag_called_allele);
    for i = 1:N
        if (isa(tag_called_allele(i).allele{1}, 'AlignedSEQ'))
            event = tag_called_allele(i).allele{1}.get_event_structure;
            if (startsWith(event, 'E') || endsWith(event, 'E'))
                allele_printout = 'collapsed_allele';
            else
                try
                    allele_printout = Mutation.ToAnnotate( ...
                        summary.CARLIN_def, tag_called_allele(i).allele(1));
                    allele_printout_str = strjoin(allele_printout, ';');
                catch
                allele_printout = 'collapsed_allele';
                end
            end
        end
        fprintf(fid, '%s\t%s\n', tag_called_allele(i).CB, ...
            allele_printout_str);
    end
    fclose(fid);


    %%% output MinMismatches.txt
    fprintf('Generating mismatch files\n');
    fid = fopen(append(outdir,'/nMismatches_by_CB.txt'), 'wt'); 
    fprintf(fid, 'Mapped CBs\tnMismatches\n');
    for i = 1:length(tag_called_allele)
        fprintf(fid, '%s\t %d\n', tag_called_allele(i).CB, tag_called_allele(i).min_mismatches);
    end
    fclose(fid);


end
