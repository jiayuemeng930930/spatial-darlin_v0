 function [SEQ_ind, SEQ_weight] = get_SEQ_ind_by_UMI(obj, UMI, filter)
            
      if (size(UMI,1) > 1)                
           %assert(iscell(UMI));
      else
           UMI = cell(UMI);
      end
            
      if (nargin == 3)
           filter = intersect(obj.masks.valid_lines, filter);
      else
           filter = obj.masks.valid_lines;
      end
            
            [is, which_UMI] = ismember(UMI, obj.UMI);            
            if (any(is))
                filter = filter(ismember(obj.read_UMI(filter), which_UMI(is)));
                [SEQ_ind, ~, SEQ_weight] = unique([obj.read_UMI(filter) obj.filter], 'rows');
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
            
 end
        
    