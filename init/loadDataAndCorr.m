function [p, obs] = loadDataAndCorr(p, files, eph, obs)

% PPP data
if isfield(files, 'ssr') && ~isempty(files.ssr) % If PPP, parse the iono correction
    [p.orbit_dict, p.clock_dict] = parseSsrObtClkCorr(files.ssr);
    p.vtec_dict = parseSsrVtec(files.vtec);
    p.code_bia = parser_bia(p,files.code_bias);
    p.bia_type = 1;

    matname = files.ustec_data + "USTEC.mat";
    if exist(matname,'file')==2
        load(matname);
        p.USTEC = USTEC;
    else
        load_USTEC(p.t(1),p.t(end),files.ustec_data);
        USTEC = parser_USTEC(files.ustec_data);  
        p.USTEC = USTEC;
        save(matname,'USTEC')
    end
end

if p.L2enable == true && p.bia_type == 1
    if p.enableGPS == 1
        Factor = p.L2freq^2/(p.L1freq^2-p.L2freq^2);
        [n,r] = size(obs.gps(1).data.P);
        obs.Iono_GPS = zeros(n,r);
        for i = 1:n
            bia_L1 = p.code_bia.GPS.bia_C1C(i);
            bia_L2 = p.code_bia.GPS.bia_C2L(i);
            if ~isnan(bia_L2) && ~isnan(bia_L1)
                for j = 1:r-400
                    data_L1 = obs.gps(1).data.P(i,j:j+400);
                    data_L2 = obs.gps(2).data.P(i,j:j+400);
                    if isempty(find(data_L2 == 0, 1)) && isempty(find(data_L1 == 0, 1))
                        obs.Iono_GPS(i,j) = Factor * ...
                            mean(data_L2 - p.c*bia_L2*1e-9 - data_L1 + p.c*bia_L1*1e-9);
                    end
                end
            end
        end
    end
    if p.enableGAL == 1
        Factor = p.E5bfreq^2/(p.E1freq^2-p.E5bfreq^2);
        [n,r] = size(obs.gal(1).data.P);
        obs.Iono_GAL = zeros(n,r);
        for i = 1:n
            bia_L1 = p.code_bia.GAL.bia_C1C(i);
            bia_L2 = p.code_bia.GAL.bia_C7Q(i);
            if ~isnan(bia_L2) && ~isnan(bia_L1)
                for j = 1:r-400
                    data_L1 = obs.gal(1).data.P(i,j:j+400);
                    data_L2 = obs.gal(2).data.P(i,j:j+400);
                    if isempty(find(data_L2 == 0, 1)) && isempty(find(data_L1 == 0, 1))
                        obs.Iono_GAL(i,j) = Factor * ...
                            mean(data_L2 - p.c*bia_L2*1e-9 - data_L1 + p.c*bia_L1*1e-9);
                    end
                end
            end
        end
    end
    if p.enableBDS == 1
        Factor = p.B2afreq^2/(p.B1freq^2-p.B2afreq^2);
        [n,r] = size(obs.bds(1).data.P);
        obs.Iono_BDS = zeros(n,r);
        for i = 1:n
            bia_L1 = p.code_bia.BDS.bia_C2I(i);
            bia_L2 = p.code_bia.BDS.bia_C7I(i);
            if ~isnan(bia_L2) && ~isnan(bia_L1)
                for j = 1:r-400
                    data_L1 = obs.bds(1).data.P(i,j:j+400);
                    data_L2 = obs.bds(2).data.P(i,j:j+400);
                    if isempty(find(data_L2 == 0, 1)) && isempty(find(data_L1 == 0, 1))
                        obs.Iono_BDS(i,j) = Factor * ...
                            mean(data_L2 - p.c*bia_L2*1e-9 - data_L1 + p.c*bia_L1*1e-9);
                    end
                end
            end
        end
    end
end
p.eph_b = eph; p.obs_b = [];
% Base station data
if isfield(files, 'data_base')
    % Get observables data (.obs file, RINEX verion 3.03)
    p.obs_b = parserGnssObs(p, files.data_base);
end


end