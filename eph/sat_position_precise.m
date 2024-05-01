function sat_prec = sat_position_precise(orbit_corr,sat_pos_ecef,sat_v_ecef,t_tsm)
% Compute precise satellit postion
% Input:
%       p: parameters
%       prn: SVID
%       igs_idx: the index in IGS data
%       t_tsm: time of transmit


% orbit correction in radial along cross track direction
dP_rac = [orbit_corr.x; orbit_corr.y; orbit_corr.z];
% orbit velocity correction in radial along cross track direction
dV_rac = [orbit_corr.dx; orbit_corr.dy; orbit_corr.dz];

% theta = p.omge * limit_tgps(t_tsm - p.IGS.orbit_iTOW(igs_idx));
% R = [ cos(theta)  sin(theta)  0;
%      -sin(theta)  cos(theta)  0;
%       0               0       1];
% dP_rac = R*dP_rac;
% dV_rac = R*dV_rac;
% Compute position error
deph = limit_tgps(t_tsm - posixtime(orbit_corr.datetime));
dP_rac = dP_rac + dV_rac*deph;

dP_ecef=RAC2ECEF(dP_rac,sat_pos_ecef,sat_v_ecef);

sat_prec = sat_pos_ecef-dP_ecef;

end