function [satlog] = satpost_corrpsedR_singlefreq(p,eph,obs,ind,len_prn,sys_type)
% This function is to compute the satellite positions and correct 
% the psedoranges by sat clock bias
% Input: 
%       p --parameters
%       eph --ephemeris data
%       obs --observables
%       ind --index of observables data
%       sys_type --The system that be computed, 'gps', 'gal', 'glo' and 'bds'
% Outpu:
%       satlog.svprn_mark -- Mark the sat prn that be computed
%       satlog.s_pos_ecef -- Satellite position in ECEF frame
%       satlog.corr_range -- corrected pseudorange
%----------------------------%
% Initialize
satlog.svprn_mark = zeros(len_prn,1);
satlog.prn_record = zeros(len_prn,1);
satlog.s_pos_ecef = zeros(3,len_prn);
satlog.s_pos_prc = zeros(3,len_prn);
satlog.s_v_ecef = zeros(3,len_prn);
satlog.corr_range = zeros(len_prn,1);
satlog.phase_m = zeros(len_prn,1);
satlog.wavelength = zeros(len_prn,1);
satlog.doppler = zeros(len_prn,1);
satlog.tp = zeros(len_prn,1);
obs_tr_gps = obs.tr_posix(ind);
obs_tr = obs_tr_gps;
if p.post_mode == p.mode_ppp && ~isempty(p.orbit_dict) && ~isempty(p.clock_dict)
    [~, orbit_single] = closestKeyAndValue(p.orbit_dict, obs_tr_gps, 0, 60);
    if isempty(orbit_single)
        satlog.num_sv = 0;
        return;
    end
end

obs_range = obs.(sys_type)(1).data.P(:,ind);
doppler = obs.(sys_type)(1).data.D(:,ind);
phase_cyc = obs.(sys_type)(1).data.C(:,ind);
Strength = obs.(sys_type)(1).data.S(:,ind);
eph_info = eph.(sys_type);

switch(sys_type)
    case 'gps'
        message_duration = p.gps.message_duration;
    case 'glo'
        message_duration = p.glo.message_duration;
        obs_tr = obs_tr_gps - p.glo.lps_gps; % Correct time diff from GPS time to GLO time
    case 'gal'
        message_duration = p.gal.message_duration;
        obs_tr = obs_tr_gps - p.gal.lps_gps; % Correct time diff from GPS time to GAL time
    case 'bds'
        message_duration = p.bds.message_duration;
        obs_tr = obs_tr_gps - p.bds.lps_gps; % Correct time diff from GPS time to BDS time
end

for prn = 1 :len_prn
    if p.post_mode == 1 && (~isfield(orbit_single, sys_type)...
            || ~isKey(orbit_single.(sys_type), prn))
        continue;
    end
    orbit_corr = [];
    clock_corr_m = 0;
    if p.post_mode == 1
        orbit_corr = orbit_single.(sys_type)(prn);
        clock_corr_m = obtainSatClockCorr(p.clock_dict, obs_tr_gps, 0, 60, prn, ...
            orbit_corr.iod, sys_type);
        if isnan(clock_corr_m)
            continue;
        end
    end

    if obs_range(prn) == 0 || doppler(prn) == 0 ...
        || isnan(obs_range(prn)) || prn>size(eph_info.SV_health,1)
        continue;
    end

    if p.post_mode == p.mode_rtkfloat && (phase_cyc(prn) == 0 || isnan(phase_cyc(prn)))
        continue;
    end
    doppler_i = doppler(prn);
    phase_m = 0;
    wavelength = 0;
    tp_prime = obs_range(prn)/p.c;
    t_sv = obs_tr-tp_prime;
    tidx = ephtidx(eph_info,t_sv,prn,message_duration,sys_type,...
        p.post_mode == p.mode_ppp);
    % Check the signal strength and sv health (health message 0 means ok)
    if isempty(tidx) || Strength(prn) < p.sig_strg
        continue;
    end

    switch(sys_type)
        % compute ephemeris to satellite position and clock bias
        case 'gps'
            [sat, dt_sv, ddt_sv] = eph2pos(p,eph_info,prn,tidx,t_sv,'gps',orbit_corr);
            wavelength = p.gps.l1_wavelength;
            doppler_i = -doppler_i*wavelength + p.c*ddt_sv;
            phase_m = phase_cyc(prn)*wavelength + p.c*dt_sv;
            if ~isnan(sat.pos_ecef(1))
                satlog.svprn_mark(prn) = p.gps.sys_num;satlog.prn_record(prn) = prn;
            end
        case 'glo'
            [sat, dt_sv, ddt_sv] = geph2pos(p,eph_info,prn,tidx,t_sv);
            fcn = eph_info.freq(prn,tidx(end));
            wavelength = p.c/(p.L1glo+p.glo.df1*fcn);
            doppler_i = -doppler_i*wavelength + p.c*ddt_sv;
            phase_m = phase_cyc(prn)*wavelength + p.c*dt_sv;
            if ~isnan(sat.pos_ecef(1))
                satlog.svprn_mark(prn) = p.glo.sys_num;satlog.prn_record(prn) = prn;
            end
         case 'gal'
            [sat, dt_sv, ddt_sv] = eph2pos(p,eph_info,prn,tidx,t_sv,'gal',orbit_corr);
            wavelength = p.gal.e1_wavelength;
            doppler_i = -doppler_i*wavelength + p.c*ddt_sv;
            phase_m = phase_cyc(prn)*wavelength + p.c*dt_sv;
            if ~isnan(sat.pos_ecef(1))
                satlog.svprn_mark(prn) = p.gal.sys_num;satlog.prn_record(prn) = prn;
            end
        case 'bds'
            [sat, dt_sv, ddt_sv] = eph2pos(p,eph_info,prn,tidx,t_sv,'bds',orbit_corr);
            wavelength = p.bds.e1_wavelength;
            doppler_i = -doppler_i*wavelength + p.c*ddt_sv;
            phase_m = phase_cyc(prn)*wavelength + p.c*dt_sv;
            if ~isnan(sat.pos_ecef(1))
                satlog.svprn_mark(prn) = p.bds.sys_num;satlog.prn_record(prn) = prn;
            end
    end
    if ~isnan(sat.pos_ecef(1))
        satlog.tp(prn) = tp_prime+dt_sv;
        satlog.s_pos_ecef(:,prn) = sat.pos_ecef;
        satlog.s_v_ecef(:,prn) = sat.v_ecef;
        if p.post_mode == 1
            satlog.s_pos_prc(:,prn) = sat.pos_prc;
        end
        % ts = tr - rho/c - dt_sv - dt_corr + code_bias;
        satlog.corr_range(prn) = obs_range(prn)+p.c*dt_sv+clock_corr_m;
        satlog.phase_m(prn) = phase_m;
        satlog.wavelength(prn) = wavelength;
        satlog.doppler(prn) = doppler_i;
    end
end
satlog.num_sv = sum(satlog.svprn_mark~=0);
end