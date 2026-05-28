function [alleles, which_seqs, weight_contribution, allele_contribution] = call_alleles_multiple_grain(CARLIN_def, aligned_seqs, constituents)
    
    total_size = size(aligned_seqs, 1);
    aligned_seq_weights = ones(size(aligned_seqs));
    
    valid_mask = find(~cellfun(@isempty, aligned_seqs));
    aligned_seqs = aligned_seqs(valid_mask);
    aligned_seq_weights = aligned_seq_weights(valid_mask); 
    aligned_seq_umis = cellfun('size',constituents, 1);
    aligned_seq_umis = aligned_seq_umis(valid_mask);
    
    alleles = [];
    which_seqs = [];
    weight_contribution = [];
    allele_contribution = [];
    
    if (isempty(aligned_seqs))
        return;
    end

    events = cellfun(@(x) x.get_event_structure(), aligned_seqs, 'un', false);
    [events, ~, which_event] = unique_by_freq(events, aligned_seq_weights);            
    event_weight = accumarray(which_event, aligned_seq_weights);
    umi_weights = accumarray(which_event, aligned_seq_umis);
    assert(issorted(event_weight, 'descend'));
    allele_weight = event_weight / total_size;
    assert(isequal(size(allele_weight),size(umi_weights)));    
    seq = 1:size(allele_weight,1);
    umi_call_allele = horzcat(allele_weight, umi_weights,seq');
    umi_call_allele = sortrows(umi_call_allele,[1 2 3], 'descend');
 
    mask = find(umi_call_allele(:,1,:) > 0.1);
    umi_call_allele = umi_call_allele(mask,:,:);

    if (isempty(umi_call_allele))
        return;
    end
    
    N = size(mask,1);
    alleles = cell(N,1);
    which_seqs = cell(N,1);
    weight_contribution = cell(N,1);
    allele_contribution = cell(N,1);
    allele_idx = umi_call_allele(:,3,:);
    idx = [];

    for i = 1:N
        event_mask = find(which_event == allele_idx(i));
        seqs_for_event = cellfun(@(x) x.get_seq(), aligned_seqs(event_mask), 'un', false);
        umi_weights_for_event = aligned_seq_umis(event_mask);
        seq_weights_for_event = aligned_seq_weights(event_mask);
        seq_lengths = cellfun(@length, seqs_for_event);

        % Empirically we find almost all reads grouped by UMI/CB are of the same length
        % so insertions/deletions due to sequencing and RT are negligible. Although we 
        % could do a more sophisticated length-dependent multialign, in practice it's
        % not worth the trouble just to include a few more reads.
        if (length(unique(seq_lengths)) == 1)           
            temp = repelem(seqs_for_event, umi_weights_for_event);
            temp = vertcat(temp{:});
            temp_seq = repelem(seqs_for_event, seq_weights_for_event);
            if (size(temp_seq,1) > 1)
                dominant_seq = mode(temp);
            else
                dominant_seq = temp;
            end
            ref_ind = find(strcmp(temp_seq, dominant_seq));
            if (isempty(ref_ind))
                idx = [idx, i];
            else
                ref_ind = ref_ind(1);
                alleles{i} = CARLIN_def.desemble_sequence(mode(temp,1), aligned_seqs{event_mask(ref_ind)}.get_ref());
                which_seqs{i} = valid_mask(event_mask);
                weight_contribution{i} = sum(umi_weights_for_event); 
                allele_contribution{i} = umi_call_allele(i,1,:);
            end
        else
            mode_length = mode(repelem(seq_lengths, umi_weights_for_event));
            length_mask = find(seq_lengths == mode_length);
            temp = repelem(seqs_for_event(length_mask), umi_weights_for_event(length_mask));
            temp = vertcat(temp{:});
            temp_seq = repelem(seqs_for_event, seq_weights_for_event);
            if (size(temp,1) > 1)
                dominant_seq = mode(temp);
            else
                dominant_seq = temp;
            end
            ref_ind = find(strcmp(temp_seq, dominant_seq));
            if (isempty(ref_ind))
                idx = [idx, i];
            else
                ref_ind = ref_ind(1);                    
                alleles{i} = CARLIN_def.desemble_sequence(mode(temp,1), aligned_seqs{event_mask(ref_ind)}.get_ref());
                which_seqs{i} = valid_mask(event_mask(length_mask));
                weight_contribution{i} = sum(umi_weights_for_event(length_mask));
                allele_contribution{i} = umi_call_allele(i,1,:);
            end
        end
    end
    alleles(idx) =[];
    which_seqs(idx)=[];
    weight_contribution(idx)=[];
    allele_contribution(idx)=[];
end