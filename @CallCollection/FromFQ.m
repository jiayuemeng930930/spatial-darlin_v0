 function out = FromFQ(FQ)
      assert(isa(FQ, 'SpatialFastQData'));
      CBs = FQ.get_CBs();
      N = size(CBs,1);
      fprintf('Building CB collection from FastQ with %d CBs\n', N);
            
      [UMIs, filter] = FQ.get_UMIs_by_CB(CBs);
      CB_data = cell(N,1);

      for i = 1:N
        if (size(UMIs{i},1) > 1)                
            assert(iscell(UMIs{i}));
        else
            UMI = cell(UMIs{i});
        end
        
        if (nargin == 3)
            filter = intersect(obj.masks.valid_lines, filter);
        else
            filter = obj.masks.valid_lines;
        end

        [is, which_UMI] = ismember(UMIs{i}, obj.UMI);            
        if (any(is))
            filter = filter(ismember(FQ.read_UMI(filter), which_UMI(is)));
            [SEQ_ind, ~, SEQ_weight] = unique([obj.read_UMI(filter) filter], 'rows');
            SEQ_weight = accumarray(SEQ_weight,1);
            [SEQ_ind, SEQ_weight] = ...
                arrayfun(@(i) deal(SEQ_ind(SEQ_ind(:,1)==i,2), SEQ_weight(SEQ_ind(:,1)==i)), which_UMI(is), 'un', false);
            SEQ_ind = cellfun(@(x,y) sortrows([x, y], 2, 'descend'), SEQ_ind, SEQ_weight, 'un', false);
            [SEQ_ind, SEQ_weight] = cellfun(@(x) deal(x(:,1), x(:,2)), SEQ_ind, 'un', false);
            SEQ_ind(is) = SEQ_ind;
            SEQ_weight(is) = SEQ_weight;
        end
            
        SEQ_ind(~is) = {[]};
        SEQ_weight(~is) = {[]};
        UMI_data = cellfun(@(u,i,w) UMIData(u,i,w), UMIs{i}, SEQ_ind, SEQ_weight, 'un', false);
        CB_data{i} = CBData(CBs{i}, vertcat(UMI_data{:}));
      end
      out = CBCollection(vertcat(CB_data{:}));
 end

     