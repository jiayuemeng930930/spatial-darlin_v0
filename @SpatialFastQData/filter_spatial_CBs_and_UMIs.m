function masks = filter_spatial_CBs_and_UMIs(cfg, CB, read_CB, UMI, read_UMI, QC, linker, polyA)

    fprintf('Analyzing Spatial CBs and UMIs\n');

    assert(strcmp(cfg.type, 'Spatial'), 'Invalid CFG passed function');
    assert(size(read_CB,1) == size(QC,1), 'Dimension mismatch QC and CB');
    assert(size(read_UMI,1) == size(QC,1), 'Dimension mismatch QC and UMI');
    
    % exact match of the up_primer_sequence
    %masks.primer_site_match = cellfun(@(x) eq(x,"TCTTCAGCGTTCCCGAGA"), linker);

    % fuzzy match of the up_primer_sequence
    hamming_score = cellfun(@(x) sum('TCTTCAGCGTTCCCGAGA' ~= x), linker);
    masks.primer_site_match = le(hamming_score, 3);
    
    % CB_UMI with QC > 18
    %masks.good_CB_UMI_QC = cellfun(@(x) all(double(x)-33 >= 18), QC);

    masks.CB_no_N = ~cellfun(@(x) any(x=='N'), CB);
    masks.CB_no_N = masks.CB_no_N(read_CB);

    masks.UMI_no_N = ~cellfun(@(x) any(x=='N'), UMI);
    masks.UMI_no_N = masks.UMI_no_N(read_UMI);

    % filter 3' poly-T based on criteria
    if (cfg.CB_trim.PolyA == "exact")
        % exact match on the 3' poly-T (VVTTTTTTTT)
        masks.three_prime_match = cellfun(@(x) (size(regexp(x,'[ACG][ACG]TTTTTTTT'),1)>0), polyA);
    elseif (cfg.CB_trim.PolyA == "ignore")
        masks.three_prime_match = ~cellfun(@isempty, polyA);
    end

    %masks.valid_provenance_structure = ( masks.three_prime_match    & ...
        %masks.primer_site_match & masks.good_CB_UMI_QC & masks.CB_no_N & ...
        %masks.UMI_no_N);

    masks.valid_provenance_structure = ( masks.three_prime_match    & ...
        masks.primer_site_match & masks.CB_no_N & ...
        masks.UMI_no_N);

    masks.primer_site_match          = uint32(find(masks.primer_site_match));
    masks.three_prime_match          = uint32(find(masks.three_prime_match));
    %masks.good_CB_UMI_QC             = uint32(find(masks.good_CB_UMI_QC));
    masks.CB_no_N                    = uint32(find(masks.CB_no_N));
    masks.UMI_no_N                   = uint32(find(masks.UMI_no_N));
    masks.valid_provenance_structure = uint32(find(masks.valid_provenance_structure));

    N = length(read_CB);

    %fprintf('From %d reads, found good (primer_site_match,three_primer_match,QC,no N CB,no N UMI,all) (%d,%d,%d,%d,%d,%d) times\n', ...
        %N, length(masks.primer_site_match), length(masks.three_prime_match), ... 
        %length(masks.good_CB_UMI_QC), length(masks.CB_no_N), ...
        %length(masks.UMI_no_N), length(masks.valid_provenance_structure));

    fprintf('From %d reads, found good (primer_site_match,three_primer_match,no N CB,no N UMI,all) (%d,%d,%d,%d,%d) times\n', ...
        N, length(masks.primer_site_match), length(masks.three_prime_match), ... 
        length(masks.CB_no_N), ...
        length(masks.UMI_no_N), length(masks.valid_provenance_structure));

end
    