function cpt = trop_iono_compute(p,eph,cpt,obs,re_pos,tdoy,user_gpst,rovert_posix)
% Compute tropospheric delay and iono delay for the measurements

lla = ecef2lla(re_pos', 'WGS84');
p.lat_deg = lla(1);
p.lon_deg = lla(2);
p.h_r = lla(3);
%%%%% convert geodetic height to orthometric height
%%%%% website: https://www.mathworks.com/matlabcentral/answers/97079-how-can-i-extract-the-orthometric-height-from-the-ellpsoiidal-height-in-the-mapping-toolbox-2-5-r20
% load geoid;
%%%%% geodetic ellipsoidal separation
% N = ltln2val(geoid, geoidrefvec, p.lat*180/pi, p.lon*180/pi);

% ecef2llh longitude range is [-pi,+pi] but IGGTrop takes longitude input
% range is [0,2pi]. This mapping is done here. 
lon_deg = p.lon_deg;
if p.lon_deg <= 0
    lon_deg = mod(p.lon_deg,360);
end
N = geoidheight(p.lat_deg,lon_deg);
%%%%% orthometric height of the receiver (Groves 2.121)
H_r=p.h_r-N;
% computing tropospheric delay for the reciever (s) using various
len = length(cpt.corr_range);
ind_prn = find(cpt.svprn_mark~=0);
for i = 1:len    
    % tropo delay (meter) computation using UNB3M model
%     [cpt.trop_delay(i), ~, ~, ~, cpt.IoFac(i)]=UNB3M(p.lat,H_r,tdoy,cpt.elev(i));
%-----------------------------------%        
    % tropo delay (meter) computation using IGGTrop model
    % Reference Paper: IGGtrop_SH & IGGtrop_rH: Two Improved Empirical
    % Tropospheric Delay Models Based on Vertical Reduction Functions 
    IGGtrop_ZenithTropDelay = IGGtropSH_bl(lon_deg,p.lat_deg,H_r/1000,tdoy);
    cpt.trop_delay(i) = (1.001/sqrt(0.002001 + sin(cpt.elev(i))^2))*IGGtrop_ZenithTropDelay;
%-----------------------------------%    
    % Iono data from USTEC: https://www.ngdc.noaa.gov/stp/iono/ustec/products/    
    %---------------------------%
    % Select the frequncy
    sysi = cpt.svprn_mark(ind_prn(i));
    switch sysi
        case 1 % GPS
            freq = p.L1freq;
            freq2 = p.L2freq;
%             switch obs.GPS.f1
%                 case 'L1'
%                     freq = p.L1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
        case 2 % GLO
            freq = p.L1freq;
%             switch obs.GLO.f1
%                 case 'G1'
%                     freq = p.L1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
        case 3 % GAL
            freq = p.E1freq;
            freq2 = p.E5bfreq;
%             switch obs.GAL.f1
%                 case 'E1'
%                     freq = p.E1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
        case 4 % BDS
            freq = p.B1freq;
            freq2 = p.B2afreq;
%             switch obs.BDS.f1
%                 case 'B1'
%                     freq = p.B1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
    end
    %---------------------------%
    if p.L2enable == 1 && p.bia_type == 1
        prn = cpt.prn_record(ind_prn(i));
        switch sysi
            case 1
                iono = obs.Iono_GPS(prn,p.i);
                if iono ~= 0
                    cpt.iono_delay(i) = iono;
                end
            case 3
                iono = obs.Iono_GAL(prn,p.i);
                if iono ~= 0
                    cpt.iono_delay(i) = iono;
                end
            case 4
                iono = obs.Iono_BDS(prn,p.i);
                if iono ~= 0
                    cpt.iono_delay(i) = iono;
                end
        end
        if isnan(cpt.iono_delay(i))
            cpt.num_sv(sysi) = cpt.num_sv(sysi) - 1;
            cpt.prn_record(ind_prn(i)) = 0;
            cpt.svprn_mark(ind_prn(i)) = 0;
        end
    elseif p.post_mode == 1 && ~isempty(p.vtec_dict) && p.enable_vtec
        % Computing Iono delay using SSR VTEC
        [cpt.iono_delay(i),cpt.iono_map_m(i)] = ssrVtecComputation(p,p.vtec_dict,re_pos,cpt.elev(i),cpt.az(i),user_gpst,rovert_posix,freq);
        % iono_delay(i) = ustec_iono_delay_computation(p,p.USTEC,cpt.elev(i),cpt.az(i),user_gpst,freq);
    elseif p.post_mode == 1 && ~isempty(p.USTEC) && ~p.enable_vtec
        % Computing Iono delay
        [cpt.iono_delay(i),cpt.iono_map_m(i)] = ustec_iono_delay_computation(p,p.USTEC,cpt.elev(i),cpt.az(i),user_gpst,freq);
    else
        if ~isempty(eph.ionoParameters)
            [cpt.iono_delay(i)] = klobuchar_model(p,eph.ionoParameters,cpt.elev(i),cpt.az(i),user_gpst.sow);
        else
            cpt.iono_delay(i) = 0;
        end
    end
end

ind0 = find(isnan(cpt.iono_delay));
cpt.corr_range(ind0) = [];
cpt.s_pos_ecef(:,ind0) = [];
if p.post_mode == 1
    cpt.s_pos_prc(:,ind0) = [];
    cpt.sat_posprc_Rcorr(:,ind0) = [];
end
cpt.s_v_ecef(:,ind0) = [];
cpt.sat_pos_Rcorr(:,ind0) = [];
cpt.sat_v_Rcorr(:,ind0) = [];
cpt.tp(ind0) = [];
cpt.elev(ind0) = [];
cpt.az(ind0) = [];
cpt.trop_delay(ind0) = [];
cpt.iono_delay(ind0) = [];
cpt.IoFac(ind0) = [];