function [CB, UMI, linker, polyA, QC] = parse_slideSeq_provenance(CBs, QCs)
    
    N = size(CBs, 1);
    CB = cell(N,1);
    UMI = cell(N,1);
    linker = cell(N,1);
    polyA = cell(N,1);
    QC = cell(N,1);
    length = cellfun(@(x) strlength(x), CBs);
    for i = 1:N
        if (length(i) > 41)
            CB{i} = [CBs{i}(1:8) CBs{i}(27:32)];
            UMI{i} = CBs{i}(33:40);
            linker{i} = CBs{i}(9:26);
            polyA{i} = CBs{i}(41:end);
            QC{i} = [QCs{i}(1:8) QCs{i}(27:40)];
        else
            CB{i} = ['N'];
            UMI{i} = ['N'];
            linker{i} = ['NNNNNNNNNNNNNNNNNN'];
            polyA{i} = ['N'];
            QC{i} = ['N'];
        end
    %[CB{idx}, UMI{idx}, linker{idx}, polyA{idx}] = cellfun(@(x) deal([x(1:8) x(27:32)],...
    %    x(33:40),x(9:26),x(41:end)), CBs, 'un', false);
    end

end