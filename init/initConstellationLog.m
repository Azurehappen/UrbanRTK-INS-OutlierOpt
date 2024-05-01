function satlog = initConstellationLog(p,log)
% Mark the sat prn that be computed
gpslog.svprn_mark = zeros(log.num_obs_gps,1);glolog.svprn_mark = zeros(log.num_obs_glo,1);
gallog.svprn_mark = zeros(log.num_obs_gal,1);bdslog.svprn_mark = zeros(log.num_obs_bds,1);
% Record the prn of each system for whose satllite been used
gpslog.prn_record = zeros(log.num_obs_gps,1);glolog.prn_record = zeros(log.num_obs_glo,1);
gallog.prn_record = zeros(log.num_obs_gal,1);bdslog.prn_record = zeros(log.num_obs_bds,1);
% Satellite position in ECEF frame
gpslog.s_pos_ecef = [];glolog.s_pos_ecef = [];gallog.s_pos_ecef = [];bdslog.s_pos_ecef = [];
if p.post_mode == 1
    % precise satellite position in ECEF frame
    gpslog.s_pos_prc = [];glolog.s_pos_prc = [];gallog.s_pos_prc = [];bdslog.s_pos_prc = [];
end
% Satellite velocity in ECEF frame
gpslog.s_v_ecef = [];glolog.s_v_ecef = [];gallog.s_v_ecef = [];bdslog.s_v_ecef = [];
% Signal propagation time
gpslog.tp = [];glolog.tp = [];gallog.tp = [];bdslog.tp = [];
% corrected pseudorange
gpslog.corr_range = [];glolog.corr_range = [];gallog.corr_range = [];bdslog.corr_range = [];
% Phase
gpslog.phase_m = [];glolog.phase_m = [];gallog.phase_m = [];bdslog.phase_m = [];
% Wavelength
gpslog.wavelength = [];glolog.wavelength = [];gallog.wavelength = [];bdslog.wavelength = [];
% Doppler
gpslog.doppler = [];glolog.doppler = [];gallog.doppler = [];bdslog.doppler = [];
% The number of satellite be computed
gpslog.num_sv = 0;glolog.num_sv = 0;gallog.num_sv = 0;bdslog.num_sv = 0;

satlog.gpslog = gpslog;
satlog.glolog = glolog;
satlog.gallog = gallog;
satlog.bdslog = bdslog;

end