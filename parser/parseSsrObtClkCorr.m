function [orbit_dict, clock_dict] = parseSsrObtClkCorr(file_path)
% This function read the orbit and clock correction from a file obtained
% from real-time WHU stream
% Para:
%   file_path: the file path of the correction file
% Return: orbit dict and clock dict.

disp("Reading SSR Orbit and Clock parameters......");
orbit_dict = dictionary;
clock_dict = dictionary;

% Mapping of system designators to field
sys_map = dictionary('G', 'gps', 'R', 'glo', 'E', 'gal', 'C', 'bds');

fid = fopen(file_path, 'r');

% Initialize a variable to keep track of the current data type
data_type = '';

% Loop over lines in the file
while ~feof(fid)
    % Read line
    line = fgetl(fid);
    parts = strsplit(line);

    % Check for ORBIT or CLOCK line
    if strcmp(parts{2}, 'ORBIT') || strcmp(parts{2}, 'CLOCK')
        % Store data type
        data_type = parts{2};
        % Parse GPS time string to datetime
        dt = datetime([parts{3} ' ' parts{4} ' ' parts{5} ' ' parts{6}...
            ' ' parts{7} ' ' parts{8}], 'InputFormat',...
            'yyyy MM dd HH mm ss.S');
        dt_posix = posixtime(dt);
        num_of_sat = str2double(parts{10});
        for i = 1:num_of_sat
            line = fgetl(fid);
            line = strtrim(line);
            parts = strsplit(line);
            % Read system, PRN, and data
            sys = parts{1}(1);
            if ~isKey(sys_map, sys)
                continue;
            end
            prn = str2double(parts{1}(2:end));
            data = struct();
            data.datetime = dt;

            % Depending on the dataType fill the struct
            if strcmp(data_type, 'ORBIT')
                data.iod = str2double(parts{2});
                data.x = str2double(parts{3});
                data.y = str2double(parts{4});
                data.z = str2double(parts{5});
                data.dx = str2double(parts{6});
                data.dy = str2double(parts{7});
                data.dz = str2double(parts{8});

                % Store data in inner dictionary
                if ~isConfigured(orbit_dict) || ~isKey(orbit_dict, dt_posix)
                    orbit_dict(dt_posix) = struct;
                end
                if ~isfield(orbit_dict(dt_posix), sys_map(sys))
                    orbit_dict(dt_posix).(sys_map(sys)) = dictionary;
                end
                orbit_dict(dt_posix).(sys_map(sys))(prn) = data;

            elseif strcmp(data_type, 'CLOCK')
                data.iod = str2double(parts{2});
                data.corr = str2double(parts{3});
                data.vel = str2double(parts{4});
                data.acc = str2double(parts{5});

                % Store data in inner dictionary
                if ~isConfigured(clock_dict) || ~isKey(clock_dict, dt_posix)
                    clock_dict(dt_posix) = struct;
                end
                if ~isfield(clock_dict(dt_posix), sys_map(sys))
                    clock_dict(dt_posix).(sys_map(sys)) = dictionary;
                end
                clock_dict(dt_posix).(sys_map(sys))(prn) = data;
            end
        end
    end
end

% Close the file
fclose(fid);
disp("Done.");
