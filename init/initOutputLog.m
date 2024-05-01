function log = initOutputLog(p, obs)

N = length(obs.tr_sow); % The number of positioning points
% Initialize output
log.epoch_t = datetime.empty;
log.gps_sec = [];
log.gpst = obs.tr_sow-obs.tr_sow(1);
log.err = NaN(1,N); % The position (Norm) error between estimated pos and true pos
log.hor_err = NaN(1,N); % The horizontal position (Norm) error between estimated pos and true pos
log.ned_err_norm = NaN(1,N); % NED frame norm error
log.ned_err = NaN(3,N); % NED frame error
log.cost = NaN(1,N); % Cost in MAP estimate
log.cost_bcd = NaN(1,N); % Cost in MAP estimate for BCD
log.pos_ecef = NaN(3,N); % Estimated position in ECEF
log.rover_clk = NaN(1,N); % Receiver clock bias (meter)
log.open_sky = NaN(1,N); % Open sky indicator
log.comp_time = NaN(1,N); % Computation time for each epoch
log.comp_time_bcd = NaN(1,N); % Computation time for BCD
log.sv_num_GPS = NaN(1,N); % The amount of GPS satellite be used
log.sv_num_GLO = NaN(1,N); % The amount of GLO satellit  used
log.sv_num_GAL = NaN(1,N); % The amount of GAL satellite be used
log.sv_num_BDS = NaN(1,N); % The amount of BDS satellite be used
log.num_meas_used = NaN(1,N); % The amount of measurements be used (with Doppler)
log.clk_gps = NaN(1,N);
log.clk_glo = NaN(1,N);
log.clk_gal = NaN(1,N);
log.clk_bds = NaN(1,N);
log.GDOP = NaN(1,N);
if p.state_mode ~= p.pva_mode
    log.vel_ned = NaN(3,N);
    log.vel_ned_err = NaN(3,N);
    log.vel_ned_cov = NaN(3,N);
end

if p.double_diff == false
    numOfState = p.modeToNumUserErrStates(p.state_mode) + 2; % user_states,clk,clk_drift
else
    numOfState = p.modeToNumUserErrStates(p.state_mode);
end
if p.enableGPS == 1 && p.double_diff == false
    % Error state size;
    log.state_cov = NaN(numOfState+p.enableGLO+p.enableGAL+p.enableBDS, N);
    log.state_info = log.state_cov;
else
    % Only used to test a single constellation.
    log.state_cov = NaN(numOfState, N);
    log.state_info = log.state_cov;
end
log.ned_cov = NaN(3, N);
if p.est_mode == p.raps_est || p.est_mode == p.raps_ned_est
    log.pos_risk = NaN(1,N); % Pos risk in MAP estimate
    log.pos_risk_bcd = NaN(1,N); % Pos risk in MAP estimate for BCD
    log.raps_num_iter = NaN(1,N); % Number of iteration of BCD
    log.raps_constraint = NaN(6,N);
    log.raps_num_sat = NaN(1,N);
    log.raps_flag = NaN(1,N);
    log.raps_penalty = NaN(1,N);
    % log.raps_num_nodes = NaN(1,N);
    log.raps_spec_xyz = NaN(3, N);
    log.pos_info_ned = NaN(3, N);
end
if ~isempty(obs.gps)
    log.num_obs_gps = size(obs.gps(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gps = 0;
end
if ~isempty(obs.glo)
    log.num_obs_glo = size(obs.glo(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_glo = 0;
end
if ~isempty(obs.gal)
    log.num_obs_gal = size(obs.gal(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gal = 0;
end
if ~isempty(obs.bds)
    log.num_obs_bds = size(obs.bds(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_bds = 0;
end
log.res_GPS = NaN(log.num_obs_gps,N); % The residual at the end
log.res_GLO = NaN(log.num_obs_glo,N); % The residual at the end
log.res_GAL = NaN(log.num_obs_gal,N); % The residual at the end
log.res_BDS = NaN(log.num_obs_bds,N); % The residual at the end
log.elev_GPS = NaN(log.num_obs_gps,N); % The elevation of satellites
log.elev_GLO = NaN(log.num_obs_glo,N);
log.elev_GAL = NaN(log.num_obs_gal,N);
log.elev_BDS = NaN(log.num_obs_bds,N);
log.res = [log.res_GPS;log.res_GAL;log.res_GLO;log.res_BDS];
log.elev = [log.elev_GPS;log.elev_GAL;log.elev_GLO;log.elev_BDS];

end