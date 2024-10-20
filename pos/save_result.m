function log = save_result(p,cpt,log,i,estState,res,grd,epoch_t,reset_flag)

%------------------------%
% [pos_llh,~,~]=ecef2llh_iter(estState.pos);
% R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
% err_pos = grdpos - estState.pos;
% ned_err=R_e2g*err_pos;
lla_gt_deg = ecef2lla(grd.pos', 'WGS84');
wgs84 = wgs84Ellipsoid('meter');
pos = estState.pos;
log.pos_ecef(:,i) = pos;
[xNorth,yEast,zDown] = ecef2ned(pos(1),pos(2),pos(3),lla_gt_deg(1),lla_gt_deg(2),lla_gt_deg(3),wgs84);
ned_err = [xNorth;yEast;zDown];
log.ned_err(:,i) = ned_err;
%------------------------%
log.open_sky(i) = cpt.is_open_sky;
log.ned_err_norm(i) = norm(ned_err);
log.hor_err(i) = norm(ned_err(1:2));
log.err(i) = norm(grd.pos - estState.pos);
log.rover_clk(i) = estState.clock_bias;
log.sv_num_GPS(i) = cpt.num_sv(1);log.sv_num_GLO(i) = cpt.num_sv(2);
log.sv_num_GAL(i) = cpt.num_sv(3);log.sv_num_BDS(i) = cpt.num_sv(4);
log.num_meas_used(i) = p.num_meas_used;
ind_mark = cpt.svprn_mark ~= 0;
% log.res(ind_mark,i) = res;
log.elev(ind_mark,i) = cpt.elev;
start = 1; endi = log.num_obs_gps;
log.res_GPS(:,i) = log.res(start:endi,i);
log.elev_GPS(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_gps; endi = endi + log.num_obs_glo;
log.res_GLO(:,i) = log.res(start:endi,i);
log.elev_GLO(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_glo; endi = endi + log.num_obs_gal;
log.res_GAL(:,i) = log.res(start:endi,i);
log.elev_GAL(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_gal; endi = endi + log.num_obs_bds;
log.res_BDS(:,i) = log.res(start:endi,i);
log.elev_BDS(:,i) = log.elev(start:endi,i);
log.clk_gps(i) = estState.clock_sys(p.gps.sys_num);
log.clk_glo(i) = estState.clock_sys(p.glo.sys_num);
log.clk_gal(i) = estState.clock_sys(p.gal.sys_num);
log.clk_bds(i) = estState.clock_sys(p.bds.sys_num);
log.GDOP(i) = p.GDOP;
log.state_cov(:,i) = diag(p.state_cov)';
log.comp_time(i) = p.comp_t;
log.comp_time_bcd(i) = p.comp_t_bcd;
R_e2g=computeRotForEcefToNed(lla_gt_deg');
if p.state_mode ~= p.pos_mode
    log.vel_ned(:,i) = R_e2g*estState.vel;
    vel_ned_cov = R_e2g * p.state_cov(4:6, 4:6) * R_e2g';
    log.vel_ned_cov(:,i) = diag(vel_ned_cov);
    if isfield(grd, 'vel')
        log.vel_ned_err(:,i) = grd.vel - log.vel_ned(:,i);
    end
end
ned_cov = R_e2g * p.state_cov(1:3, 1:3) * R_e2g';
log.ned_cov(:,i) = diag(ned_cov);
if p.est_mode ~= p.ekf_est
    log.cost(i) = p.augcost;
    log.cost_bcd(i) = p.augcost_bcd;
    log.state_info(:,i) = diag(p.infor_ned)';
end
if (p.est_mode == p.raps_est || p.est_mode == p.raps_ned_est) && reset_flag == false
    log.pos_risk(i) = p.pos_risk;
    log.pos_risk_bcd(i) = p.pos_risk_bcd;
    log.raps_penalty(i) = p.raps_penalty;
    log.raps_flag(i) = p.raps_flag;
    %log.raps_constraint(1:6,i) = p.constraint;
    log.raps_num_sat(i) = p.raps_num_sat;
    log.raps_num_iter(i) = p.raps_num_iter;
    % log.raps_spec_xyz(:,i) = p.raps_spec_xyz;
    % J_ned = R_e2g' * p.raps_J(1:3, 1:3) * R_e2g;
    % log.pos_info_ned(:,i) = diag(p.raps_J(1:3, 1:3));
end
end