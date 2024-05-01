function satInfo = obtainSatInfoStruct(p,satlog,eph,obs,ind,log)
% Initialize a satellite information struct
satInfo = struct;
if ~isempty(obs.gps)&& p.freq==1 && p.enableGPS == 1
    % GPS satellite position computation, Single frenquency receiver mode
    satlog.gpslog = satpost_corrpsedR_singlefreq(p,eph,obs,ind,log.num_obs_gps,'gps');
end
if ~isempty(obs.glo)&& p.freq==1 && p.enableGLO == 1
    % GLO satellite position computation, Single frenquency receiver mode
    satlog.glolog = satpost_corrpsedR_singlefreq(p,eph,obs,ind,log.num_obs_glo,'glo');
end
if ~isempty(obs.gal)&& p.freq==1 && p.enableGAL == 1
    % GAL satellite position computation, Single frenquency receiver mode
    satlog.gallog = satpost_corrpsedR_singlefreq(p,eph,obs,ind,log.num_obs_gal,'gal');
end
if ~isempty(obs.bds)&& p.freq==1 && p.enableBDS == 1
    % BDS satellite position computation, Single frenquency receiver mode
    satlog.bdslog = satpost_corrpsedR_singlefreq(p,eph,obs,ind,log.num_obs_bds,'bds');
end
gpslog = satlog.gpslog;
glolog = satlog.glolog;
gallog = satlog.gallog;
bdslog = satlog.bdslog;
satInfo.prn_record = [gpslog.prn_record;glolog.prn_record;gallog.prn_record;bdslog.prn_record];
satInfo.svprn_mark = [gpslog.svprn_mark;glolog.svprn_mark;gallog.svprn_mark;bdslog.svprn_mark];
satInfo.corr_range = [gpslog.corr_range;glolog.corr_range;gallog.corr_range;bdslog.corr_range];
satInfo.phase_m = [gpslog.phase_m;glolog.phase_m;gallog.phase_m;bdslog.phase_m];
satInfo.wavelength = [gpslog.wavelength;glolog.wavelength;gallog.wavelength;bdslog.wavelength];
satInfo.doppler = [gpslog.doppler;glolog.doppler;gallog.doppler;bdslog.doppler];
ind = find(satInfo.prn_record==0);
satInfo.prn_record(ind) = [];
satInfo.svprn_mark(ind) = [];
ind = find(satInfo.corr_range==0);
satInfo.corr_range(ind) = [];
satInfo.doppler(ind) = [];
satInfo.phase_m(ind) = [];
satInfo.wavelength(ind) = [];
satInfo.s_pos_ecef = [gpslog.s_pos_ecef,glolog.s_pos_ecef,gallog.s_pos_ecef,bdslog.s_pos_ecef];
satInfo.s_pos_ecef(:,ind) = [];
if p.post_mode == 1
    satInfo.s_pos_prc = [gpslog.s_pos_prc,glolog.s_pos_prc,gallog.s_pos_prc,bdslog.s_pos_prc];
    satInfo.s_pos_prc(:,ind) = [];
end
satInfo.s_v_ecef = [gpslog.s_v_ecef,glolog.s_v_ecef,gallog.s_v_ecef,bdslog.s_v_ecef];
satInfo.s_v_ecef(:,ind) = [];
satInfo.tp = [gpslog.tp;glolog.tp;gallog.tp;bdslog.tp];
satInfo.tp(ind) = [];
satInfo.num_sv = [gpslog.num_sv,glolog.num_sv,gallog.num_sv,bdslog.num_sv];

% elevation & azimuth
satInfo.elev = NaN(length(satInfo.corr_range),1); satInfo.az = NaN(length(satInfo.corr_range),1);
% trop delay and iono delay
satInfo.trop_delay = NaN(length(satInfo.corr_range),1); satInfo.iono_delay = NaN(length(satInfo.corr_range),1);
satInfo.iono_map_m = NaN(length(satInfo.corr_range),1);
end