function obs = parserGnssObs(p, obspath)
% Parse observation data from rinex file to matlab data file.
% Supportting RINEX version 3
%
%%%%%-----Reference
% ftp://igs.org/pub/data/format/rinex303.pdf
%
%%%%%-----Input
% obs file path
%
%%%%%-----Output
% Struct of constellation observation
%
% Specification
%           Reference: Section 5 in RINEX document
%           Output obs will include GPS, GAL, GLO and BDS
%           In this statement, attribute 'sys' denote GPS/GLO/GAL/BDS
%           obs.sys(i).type : The frequency type and measurement code
%           obs.sys(i).data : P:psedorange (meter), 
%                             C: carrier phase (cycle),
%                             D: doppler (Hz), 
%                             S: signal strength
%           Default:
%                   GPS(1)->type/data: (C1C) L1 frequency, C/A code;
%                   GPS(2)->type/data: (C2L/C2S/C2X) L2 frequency;
%                   GPS(3)->type/data: (C1W) L1 frequency, P code;
%                   GPS(4)->type/data: (C2W) l2 frequency, P code;
%
%                   GLO(1)->type/data: (C1C) L1 frequency, C/A code;
%                   GLO(2)->type/data: (C2C) L2 frequency, C/A code;
%                   GLO(3)->type/data: (C1P) L1 frequency, P code;
%                   GLO(4)->type/data: (C2P) l2 frequency, P code;
%
%                   GAL(1)->type/data: (C1X/C1C) E1 frequency;
%                   GAL(2)->type/data: (C7I/C7Q/C7X) E5b frequency;
%
%                   BDS(1)->type/data: (C2I/C2Q/C2X) B1(B1-2 in RINEX 3.04) frequency;
%                   BDS(2)->type/data: (C7I/C7Q/C7X) B2b frequency;
%           see RINEX 3.03 table A2
%
% Author: Azurehappen
fprintf ('Loading observations...\n');
obsfile = fopen(obspath);
%-----------------------------------%
% read header
fprintf ('Reading header...\n');
while (true)
    line = fgetl(obsfile);                                                   %get line
    if contains(line,'SYS / # / OBS TYPES')
        constellation = line(1);
        switch(constellation)
            case 'G'
                fprintf('File contains GPS observations \n')
                num_gpstype  = str2double(line(5:6)); % number of obs type
                gpstype = strsplit(line(8:58));
                gpstype(cellfun(@isempty,gpstype))=[]; % Delete empty cite
                i = length(gpstype);
                while (i<num_gpstype)
                    line = fgetl(obsfile);
                    insert = strsplit(line(8:58));
                    insert(cellfun(@isempty,insert))=[]; % Delete empty cite
                    gpstype = [gpstype insert];
                    i = length(gpstype);
                end
            case 'R'
                fprintf('File contains GLO observations \n')
                num_glotype  = str2double(line(5:6)); % number of obs type
                glotype = strsplit(line(8:58));
                glotype(cellfun(@isempty,glotype))=[]; % Delete empty cite
                i = length(glotype);
                while (i<num_glotype)
                    line = fgetl(obsfile);
                    insert = strsplit(line(8:58));
                    insert(cellfun(@isempty,insert))=[]; % Delete empty cite
                    glotype = [glotype insert];
                    i = length(glotype);
                end
            case 'E'
                fprintf('File contains GAL observations \n')
                num_galtype  = str2double(line(5:6)); % number of obs type
                galtype = strsplit(line(8:58));
                galtype(cellfun(@isempty,galtype))=[]; % Delete empty cite
                i = length(galtype);
                while (i<num_galtype)
                    line = fgetl(obsfile);
                    insert = strsplit(line(8:58));
                    insert(cellfun(@isempty,insert))=[]; % Delete empty cite
                    galtype = [galtype insert];
                    i = length(galtype);
                end
            case 'C'
                fprintf('File contains BDS observations \n')
                num_bdstype  = str2double(line(5:6)); % number of obs type
                bdstype = strsplit(line(8:58));
                bdstype(cellfun(@isempty,bdstype))=[]; % Delete empty cite
                i = length(bdstype);
                while (i<num_bdstype)
                    line = fgetl(obsfile);
                    insert = strsplit(line(8:58));
                    insert(cellfun(@isempty,insert))=[]; % Delete empty cite
                    bdstype = [bdstype insert];
                    i = length(bdstype);
                end
        end
    elseif contains(line,"APPROX POSITION XYZ")
        aprox_pos_str = strsplit(line(1:58));
        aprox_pos_str(cellfun(@isempty,aprox_pos_str))=[]; % Delete empty cite
        aprox_pos = str2double(aprox_pos_str');

    elseif contains(line, 'TIME OF FIRST OBS')
        start_t = sscanf(line,'%f');
    elseif contains(line, 'TIME OF LAST OBS')
        end_t = sscanf(line,'%f');
    elseif contains(line,'END OF HEADER')
        break;                                                              % End of header loop
    end
end
No_obs = round(seconds(datetime(end_t')-datetime(start_t')))+10;
% Initialization
obs = initObsStruct(p, No_obs);
obs.aprox_pos = aprox_pos;
%-----------------------------------%
% read observables
count = 0;
% Step size for each data and its gap
gap = 4; len = 12;
fprintf ('Parsing observables \n');
while ~feof(obsfile)
    line = fgetl(obsfile);
    if strcmp(line(1),'>') % new observables
        count = count+1;
        lineSplit = strsplit(line);
        % [year;month;date;hour;minute;second]
        obs.tr_prime(:,count) = str2double(lineSplit(2:7))';
        obs.datetime(1, count) = datetime(obs.tr_prime(:,count)');
        obs.tr_posix(1, count) = posixtime(obs.datetime(1, count));
        [obs.tr_week(1,count),~,obs.tr_sow(1,count)] = date2gnsst(str2double(lineSplit(2:7))); %  GPS seconds
    else
        switch line(1)
            case{'G'}
                prn = str2double(line(2:3));
                idx = 1+gap;
                for i = 1:length(gpstype)
                    if idx + len <=length(line)
                        % Avioding short line that don't have other types of data
                        switch(gpstype{i}(1))
                            case 'C' % Psedorange
                                if strcmp(gpstype{i},'C1C')
                                    val=str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.gps(1).data.P(prn,count)=val;
                                    end
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'C2X')
                                    val = str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.gps(2).data.P(prn,count)=str2double(line(idx:idx+len));
                                    end
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'C1W')
                                    obs.gps(3).data.P(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'C2W')
                                    obs.gps(4).data.P(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'L' % Carrier phase
                                if strcmp(gpstype{i},'L1C')
                                    obs.gps(1).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'L2X')
                                    obs.gps(2).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'L1W')
                                    obs.gps(3).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'L2W')
                                    obs.gps(4).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'D' % Doppler
                                if strcmp(gpstype{i},'D1C')
                                    obs.gps(1).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'D2X')
                                    obs.gps(2).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'D1W')
                                    obs.gps(3).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'D2W')
                                    obs.gps(4).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'S' % Signal strength
                                if strcmp(gpstype{i},'S1C')
                                    obs.gps(1).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'S2X')
                                    obs.gps(2).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'S1W')
                                    obs.gps(3).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(gpstype{i},'S2W')
                                    obs.gps(4).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                        end
                        %[obs,idx]=readgps(gpstype{i},prn,count,obs,idx,line,len);
                        idx = idx + min(gap,length(line)-idx);
                    end
                end
            case{'R'}
                prn = str2double(line(2:3));
                idx = 1+gap;
                for i = 1:length(glotype)
                    if idx + len <=length(line)
                        switch(glotype{i}(1))
                            case 'C' % Psedorange
                                if strcmp(glotype{i},'C1C')
                                    val=str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.glo(1).data.P(prn,count)=val;
                                    end
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'C2C')
                                    val = str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.glo(2).data.P(prn,count)=str2double(line(idx:idx+len));
                                    end
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'C1P')
                                    obs.glo(3).data.P(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'C2P')
                                    obs.glo(4).data.P(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'L' % Carrier phase
                                if strcmp(glotype{i},'L1C')
                                    obs.glo(1).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'L2C')
                                    obs.glo(2).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'L1P')
                                    obs.glo(3).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'L2P')
                                    obs.glo(4).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'D' % Doppler
                                if strcmp(glotype{i},'D1C')
                                    obs.glo(1).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'D2C')
                                    obs.glo(2).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'D1P')
                                    obs.glo(3).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'D2P')
                                    obs.glo(4).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'S' % Signal strength
                                if strcmp(glotype{i},'S1C')
                                    obs.glo(1).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'S2C')
                                    obs.glo(2).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'S1P')
                                    obs.glo(3).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(glotype{i},'S2P')
                                    obs.glo(4).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                        end
                        %                         [obs,idx]=readglo(glotype{i},prn,count,obs,idx,line,len);
                        idx = idx + min(gap,length(line)-idx);
                    end
                end
            case{'E'}
                prn = str2double(line(2:3));
                idx = 1+gap;
                for i = 1:length(galtype)
                    if idx + len <=length(line)
                        switch(galtype{i}(1))
                            case 'C' % Psedorange
                                if strcmp(galtype{i},'C1C')||strcmp(galtype{i},'C1X')
                                    val=str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.gal(1).data.P(prn,count)=val;
                                    end
                                    idx = idx + len;
                                elseif strcmp(galtype{i},'C7I')||strcmp(galtype{i},'C7Q')||strcmp(galtype{i},'C7X')
                                    val = str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.gal(2).data.P(prn,count)=val;
                                    end
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'L' % Carrier phase
                                if strcmp(galtype{i},'L1C')||strcmp(galtype{i},'L1X')
                                    obs.gal(1).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(galtype{i},'L7I')||strcmp(galtype{i},'L7Q')||strcmp(galtype{i},'L7X')
                                    obs.gal(2).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'D' % Doppler
                                if strcmp(galtype{i},'D1C')||strcmp(galtype{i},'D1X')
                                    obs.gal(1).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(galtype{i},'D7I')||strcmp(galtype{i},'D7Q')||strcmp(galtype{i},'D7X')
                                    obs.gal(2).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'S' % Signal strength
                                if strcmp(galtype{i},'S1C')||strcmp(galtype{i},'S1X')
                                    obs.gal(1).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(galtype{i},'S7I')||strcmp(galtype{i},'S7Q')||strcmp(galtype{i},'S7X')
                                    obs.gal(2).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                        end
                        idx = idx + min(gap,length(line)-idx);
                    end
                end
            case{'C'}
                prn = str2double(line(2:3));
                idx = 1+gap;
                for i = 1:length(bdstype)
                    if idx + len <=length(line)
                        switch(bdstype{i}(1))
                            case 'C' % Psedorange
                                if strcmp(bdstype{i},'C2I')||strcmp(bdstype{i},'C2Q')||strcmp(bdstype{i},'C1I')
                                    val=str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.bds(1).data.P(prn,count)=val;
                                    end
                                    idx = idx + len;
                                elseif strcmp(bdstype{i},'C7I')||strcmp(bdstype{i},'C7Q')||strcmp(bdstype{i},'C7X')
                                    val = str2double(line(idx:idx+len));
                                    if ~isnan(val)
                                        obs.bds(2).data.P(prn,count)=str2double(line(idx:idx+len));
                                    end
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'L' % Carrier phase
                                if strcmp(bdstype{i},'L2I')||strcmp(bdstype{i},'L2Q')||strcmp(bdstype{i},'L1I')
                                    obs.bds(1).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(bdstype{i},'L7I')||strcmp(bdstype{i},'L7Q')||strcmp(bdstype{i},'L7X')
                                    obs.bds(2).data.C(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'D' % Doppler
                                if strcmp(bdstype{i},'D2I')||strcmp(bdstype{i},'D2Q')||strcmp(bdstype{i},'D1I')
                                    obs.bds(1).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(bdstype{i},'D7I')||strcmp(bdstype{i},'D7Q')||strcmp(bdstype{i},'D7X')
                                    obs.bds(2).data.D(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                            case 'S' % Signal strength
                                if strcmp(bdstype{i},'S2I')||strcmp(bdstype{i},'S2Q')||strcmp(bdstype{i},'S1I')
                                    obs.bds(1).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                elseif strcmp(bdstype{i},'S7I')||strcmp(bdstype{i},'S7Q')||strcmp(bdstype{i},'S7X')
                                    obs.bds(2).data.S(prn,count)=str2double(line(idx:idx+len));
                                    idx = idx + len;
                                else
                                    idx = idx + len;
                                end
                        end
                        idx = idx + min(gap,length(line)-idx);
                    end
                end
        end
    end
end
obs = simplyobs(obs,count);
fclose(obsfile);
%----------------------------------------------------%
fprintf ('\nObservables loaded correctly\n \n');

end

