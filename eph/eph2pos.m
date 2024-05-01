function [sat, dt_sv, ddt_sv] = eph2pos(p, eph, prn, tidx, t_sv,sys_type,orbit_corr)
% compute ephemeris to satellite position and clock bias
% for GPS, GAL, BDS system
%%%%%% Inputs
% p - parameters
% prn - svid
% tidx - time index in eph data
% eph  - ephemeris data structure
% obs  - oberservable data structure
% t_sv - signal transmit time
%%%%%% Outputs
% sat_pos_ecef - satellite position in ECEF
% sat_v_ecef - satellite velocity in ECEF
% dt_sv - satellite clock bias correction
% Set the system parameters and compute group delay for single frequency user
%
% Author: Wang Hu
%
% Notes from RTKLIB
% satellite position and clock are values at signal transmission time
% satellite position is referenced to antenna phase center
% satellite clock does not include code bias correction (tgd or bgd)
%------------------%Define the frequency of system and getting the GROUP
%DELAY
sat.pos_ecef = NaN(3,1);
sat.v_ecef = NaN(3,1);
if p.post_mode == p.mode_ppp
    % Save the precise pos and velocity in PPP mode
    sat.pos_prc = NaN(3,1);
    %     sat.v_prc = NaN(3,1);
end
dt_sv = 0;
ddt_sv = 0;

% Find if there is PPP Orbit correction
if p.post_mode == 1 && p.IGS_enable == 1
    if strcmp(sys_type,'bds')&&(prn<=5 || prn>=59)
        return;
    end
    if length(tidx)>1
        ei = find(eph.IODE(prn,tidx) == orbit_corr.iod);
        if ~isempty(ei)
            tidx = tidx(ei(end));
        else
            return;
        end
    else
        if ~orbit_corr.iod == eph.IODE(prn,tidx)
            return;
        end
    end
else
    tidx=tidx(end);
end

switch sys_type
    case 'gps'
        sysp = p.gps;
        if p.post_mode == 1 && p.IGS_enable == 1
            code_bias = p.code_bia.GPS.bia_C1C;
            if isnan(code_bias(prn))
                return;
            end
            group_delay = code_bias(prn)*1e-9;
        else
            group_delay = eph.TGD(prn,tidx);
        end
    case 'gal'
        sysp = p.gal;
        if p.post_mode == 1 && p.IGS_enable == 1
            code_bias = p.code_bia.GAL.bia_C1C;
            if isnan(code_bias(prn))
                return;
            end
            group_delay = code_bias(prn)*1e-9;
        else
            group_delay = eph.BGD_E5a(prn,tidx);
        end
        %         if strcmp(obs.GAL.f1,'E1')
        %             group_delay = eph.BGD_E5a(prn,tidx);
        %         elseif strcmp(obs.GAL.f1,'E5a')
        %             group_delay = (p.E1freq/p.E5afreq)^2*eph.BGD_E5a(prn,tidx);
        %         elseif strcmp(obs.GAL.f1,'E5b')
        %             group_delay = (p.E1freq/p.E5bfreq)^2*eph.BGD_E5b(prn,tidx);
        %         end
    case 'bds'
        sysp = p.bds;
        % Currently RINEX 3.03 only provide info for B1, no B2a
        if p.post_mode == 1 && p.IGS_enable == 1
            code_bias = p.code_bia.BDS.bia_C2I;
            if isnan(code_bias(prn))
                return;
            end
            group_delay = code_bias(prn)*1e-9;
        else
            group_delay = eph.TGD1(prn,tidx);
        end
end
%------------------------------------------------%
% Initialize
Ek = p.Ek0;
% compute clock correction estimate
dt_sv = sat_clock(sysp, prn, tidx, eph, Ek, t_sv);
tm = t_sv - dt_sv;         % Corr. mess. trans. time
% estimate satellite position and velocity (m) & (m/s) in ECEF at corrected GPS SV transmit time
[sat.pos_ecef, sat.v_ecef, Ek] = sat_posvel(sysp, eph, tm, prn, tidx,sys_type);
[pos_s, ~, ~] = sat_posvel(sysp, eph, tm+0.001, prn, tidx,sys_type);
sat.v_ecef = (pos_s-sat.pos_ecef)/0.001;
dt_sv2 = sat_clock(sysp, prn, tidx, eph, Ek, t_sv+0.001);
ddt_sv = (dt_sv2 - dt_sv) / 0.001;
% relativistic correction (s)
dt_sv = dt_sv + sysp.F*eph.e(prn,tidx)*eph.sqrtA(prn,tidx)*sin(Ek)- group_delay;
% check for convergence
%     sat.pos_prc = sat.pos_ecef;
%     R = norm(sat_pos_old - sat.pos_ecef);
%     if R < p.satdelta
%         break;
%     end
% end
if p.post_mode == 1 && p.IGS_enable == 1
    sat.pos_prc = sat_position_precise(orbit_corr,sat.pos_ecef,sat.v_ecef,tm);
end
% p.Ek0 = Ek; % To be a initial value at next obs.
% if and(iter>p.NsatEk-1,  R > p.satdelta)
%     sat.pos_ecef = NaN(3,1);
%     sat.v_ecef = NaN(3,1);
%     warning('GNSS path length iteration failed in satellite position computation');
% end
%------------------------------------------------%
end