function ref_BCs = get_Spatial_ref_BCs(barcode_file)
    
    fprintf('Parsing CB reference list: %s\n', barcode_file);

    assert(exist(barcode_file, 'file') == 2, 'Missing CB reference file: %s', barcode_file);
    
    [~, tempfile, ext] = fileparts(barcode_file);
    if (strcmp(ext, '.gz'))
        unzip_path = tempname;
        fprintf('Unzipping CB reference file to: %s/%s\n', unzip_path, tempfile);
        gunzip(barcode_file, unzip_path);
        barcode_file = sprintf('%s/%s', unzip_path, tempfile);
        assert(exist(barcode_file, 'file') == 2, sprintf('Unzip unsuccessful'));
    end
    
    fid = fopen(barcode_file);
    assert(fid ~= -1, 'Failed to open CB reference file: %s', barcode_file);    
    ref_BCs = textscan(fid, '%s%s%s','Delimiter', ',');
    fclose(fid);
    
    if (strcmp(ext, '.gz'))
        fprintf('Removing unzipped CB reference file: %s\n', barcode_file);
        delete(barcode_file);
    end

    assert(length(unique(ref_BCs{1})) == length(ref_BCs{1}), 'Barcodes should be unique');
    fprintf('%d CBs in reference list\n', size(ref_BCs{1},1));
    
end