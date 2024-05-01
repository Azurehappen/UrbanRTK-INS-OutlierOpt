function vtec_dict = parseSsrVtec(file_path)
% This function read the VTEC parameters from a file obtained
% from real-time CNES stream
% Para:
%   file_path: the file path of the correction file
% Return: VTEC dict.

disp("Reading SSR VTEC parameters......");
vtec_dict = dictionary;

fid = fopen(file_path, 'r');

% Loop over lines in the file
while ~feof(fid)
    % Read line 
    line = fgetl(fid);
    line = strtrim(line);
    parts = strsplit(line);

    if strcmp(parts{2}, 'VTEC')
        % Parse GPS time string to datetime
        dt = datetime([parts{3} ' ' parts{4} ' ' parts{5} ' ' parts{6}...
            ' ' parts{7} ' ' parts{8}], 'InputFormat',...
            'yyyy MM dd HH mm ss.S');
        dt_posix = posixtime(dt);
        % Read next line 
        line = fgetl(fid);
        line = strtrim(line);
        parts = strsplit(line);
        data.datetime = dt;
        data.num_deg = str2double(parts{2});
        data.num_order = str2double(parts{3});
        data.height_m = str2double(parts{4});
        % 0 is one of the degrees and orders.
        data.sin_coeffs = zeros(data.num_deg+1, data.num_order+1);
        data.cos_coeffs = zeros(data.num_deg+1, data.num_order+1);
        for i = 1:data.num_deg+1
            line = fgetl(fid);
            line = strtrim(line);
            parts = strsplit(line);
            for j = 1:data.num_order+1
                data.cos_coeffs(i, j) = str2double(parts{j});
            end
        end
        for i = 1:data.num_deg+1
            line = fgetl(fid);
            line = strtrim(line);
            parts = strsplit(line);
            for j = 1:data.num_order+1
                data.sin_coeffs(i, j) = str2double(parts{j});
            end
        end

        % Store data in inner dictionary
        if ~isConfigured(vtec_dict) || ~isKey(vtec_dict, dt_posix)
            vtec_dict(dt_posix) = data;
        end
    end
end

% Close the file
fclose(fid);
disp("Done.");
