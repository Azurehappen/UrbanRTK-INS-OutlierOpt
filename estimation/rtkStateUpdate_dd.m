function [p,estState,res_rho] = rtkStateUpdate_dd(p, cpt, dt)
% RTK Instantaneous
%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;
x_minus = p.state0;
%------------------%
y_rho = cpt.corr_range;
p.num_sats_window = [p.num_sats_window(2:length(p.num_sats_window)), length(y_rho)];
y_dop = [];
num_rho = length(y_rho); % The number of measurement
num_dop = 0;
if p.state_mode ~= p.pos_mode
    y_dop = cpt.doppler;
    num_dop = length(y_dop);
end

y_phi = [];
num_phi = 0;
T_dd = [];
num_sv_rtk = cpt.num_sv;
if p.post_mode == p.mode_rtkfloat
    y_phi = cpt.phase_m;
    num_phi = length(y_phi);
    % Find the pivot satellite
    sysnumToPivotPrn = dictionary;
    T_dd = eye(num_phi); % Phase double-diff transformation matrix
    to_delete = [];
    st_ind = 1;
    lambda = cpt.wavelength;
    for i = 1:length(cpt.num_sv)
        ind_i = find(cpt.svprn_mark == i);
        if isempty(ind_i)
            continue;
        end
        elev_i = cpt.elev(ind_i);
        prn_i = cpt.prn_record(ind_i);
        [elev_max, ind_pivot] = max(elev_i);
        if elev_max < p.pivot_elev
            % No suitable pivot sat, remove this system
            to_delete = [to_delete,st_ind:st_ind+cpt.num_sv(i)-1];
            T_dd(:,st_ind:st_ind+cpt.num_sv(i)-1) = 0;
            num_sv_rtk(i) = 0;
            st_ind = st_ind + cpt.num_sv(i);
            continue;
        else
            to_delete = [to_delete,st_ind+ind_pivot(1)-1];
            T_dd(st_ind:st_ind+cpt.num_sv(i)-1,st_ind+ind_pivot(1)-1) = -1;
        end
        st_ind = st_ind + cpt.num_sv(i);
        sysnumToPivotPrn(i) = prn_i(ind_pivot(1));
    end
    T_dd(to_delete,:) = [];
    lambda(to_delete) = [];
    errx_minus = [zeros(length(p.error_state),1);zeros(length(lambda),1)];
    H_phi = zeros(num_phi, length(errx_minus));
    num_dd_meas = size(T_dd,1);
    if sum(num_dd_meas) < 3
        % p.i
        % warning('Unconfident number of sat to perform RTK, skip this epoch');
        estState.pos = [];
        if p.state_mode == p.ins_mode && p.save_ins_pos == true
            R_e2b_plus = convertQuatToRot(p.state0(7:10));
            lever_arm_e = R_e2b_plus' * p.imu_lever_arm;  % Lever arm in ECEF frame.
            estState.pos = p.state0(1:3) + lever_arm_e;
            estState.vel = p.state0(4:6);
            estState.clock_bias = NaN;
            estState.clock_drift = NaN;
        end
        res_rho = NaN;
        return;
    end
end

H_rho = zeros(num_rho, length(errx_minus));
H_dop = [];
if p.state_mode ~= p.pos_mode
    gyro_noise2dop = zeros(num_dop,1);
    H_dop = zeros(num_dop, length(errx_minus));
    res_dop = zeros(num_dop,1);
end

s_pos_ecef = cpt.sat_pos_Rcorr;
if p.post_mode == p.mode_ppp && p.IGS_enable == 1
    s_pos_ecef = cpt.s_pos_prc;
end
s_v_ecef = cpt.sat_v_Rcorr;

if p.state_mode == p.ins_mode
    % INS state information
    % x_minus represent IMU states
    % 1->3: pos, 4->6: vel, 7->10: quaternion, 11->13: acc bias,
    % 14->16: gyro bias, 17->end-1: clock, end: clock drift
    lever_arm_b = p.imu_lever_arm;
    % lever_arm_b = [0;0;0];
    R_e2b_hat = convertQuatToRot(x_minus(7:10));
    lever_arm_e = R_e2b_hat' * lever_arm_b;  % Lever arm in ECEF frame.
    ant_pos = x_minus(1:3)+lever_arm_e;
    imu_pos = x_minus(1:3);
    Omg_ebe_hat = vectorSkewSymMat(R_e2b_hat'*(-p.w_be_b));
    Omg_lever = Omg_ebe_hat*lever_arm_e;
    ant_vel = x_minus(4:6)+Omg_ebe_hat*lever_arm_e;
else
    ant_pos = x_minus(1:3);
    imu_pos = ant_pos;
end

geo_range_rho = zeros(num_rho,1);
geo_range_phi = zeros(num_phi,1);
lever_comp_range = zeros(length(y_rho),1); % lever arm compensation to the range residual
for j=1:num_rho
    geo_range_rho(j)=norm(s_pos_ecef(:,j)-imu_pos);
    los_r_i = (imu_pos-s_pos_ecef(:,j))'/geo_range_rho(j);
    % compute LOS for Doppler
    range_r = norm(s_pos_ecef(:,j)-ant_pos);
    los_r = (ant_pos-s_pos_ecef(:,j))'/range_r;
    H_rho(j,1:3)=los_r_i;
    if p.state_mode == p.ins_mode
        H_rho(j,7:9)=-los_r_i*vectorSkewSymMat(lever_arm_e);
        lever_comp_range(j) = los_r_i*lever_arm_e;
    end
    if p.state_mode == p.pva_mode
        H_dop(j,4:6)=los_r_i;
        res_dop(j) = y_dop(j) - los_r*(x_minus(4:6) - s_v_ecef(:,j));
    elseif p.state_mode == p.ins_mode
        H_dop(j,4:6) = los_r;
        H_dop(j,7:9) = -los_r*vectorSkewSymMat(Omg_lever);
        H_dop(j,13:15) = los_r*R_e2b_hat'*vectorSkewSymMat(lever_arm_b);
        gyro_noise2dop(j) = H_dop(j,13:15)*(p.imu_para.gyro_noise*eye(3,3))*H_dop(j,13:15)';
        % Ant vel has includes the lever comp for Dop.
        res_dop(j) = y_dop(j) - los_r*(ant_vel - s_v_ecef(:,j));
    end
    if p.post_mode == p.mode_rtkfloat
        H_phi(j,1:3)=los_r_i;
        if p.state_mode == p.ins_mode
            H_phi(j,7:9)=-los_r_i*vectorSkewSymMat(lever_arm_e);
        end
        geo_range_phi(j) = geo_range_rho(j);
    end
end
% measurement residual
% x_minus represent the state where 4 states for the Quat.
res_rho = y_rho - geo_range_rho -lever_comp_range;
res_phi = y_phi - geo_range_phi -lever_comp_range;

dd_phi = T_dd * res_phi; % double-diff phase
dd_rho = T_dd * res_rho; % double-diff pseudorange
dd_dop = T_dd * res_dop; % double-diff Doppler
res = [dd_rho;dd_phi;dd_dop];

if p.state_mode ~= p.pos_mode
    % Rdop = p.sigmaSquare_dop*eye(length(y_dop));
    R_dop = diag(p.sigmaSquare_dop+gyro_noise2dop);
    % Rdop_diag = p.sigmaSquare_dop+(0.5*sqrt(p.sigmaSquare_dop)./sin(cpt.elev)).^2+gyro_noise2dop;
    % Rdop = diag(Rdop_diag);
end

[R_rho, R_phi] = constructMeasNoise(p, cpt, 1);
Rt_rho = T_dd * R_rho * T_dd';
Rt_phi = T_dd * R_phi * T_dd';
Rt_dop = T_dd * R_dop * T_dd';
% dd_amb = (dd_phi - dd_rho)./lambda;
dd_amb = (dd_phi - dd_rho);
% cov_dd_amb = T_dd * diag((diag(R_phi)+diag(R_rho))./(cpt.wavelength.^2)) * T_dd';
cov_dd_amb = T_dd * diag((diag(R_phi)+diag(R_rho))) * T_dd'+100^2;
errx_minus(length(errx_minus)-length(lambda)+1:end) = dd_amb;
cov_minus = blkdiag(p.state_cov, cov_dd_amb);
R_dd = blkdiag(Rt_rho,Rt_phi,Rt_dop);
% R_dd = diag(diag(R_dd));
H_rho_dd = T_dd * H_rho;
H_phi_dd = T_dd * H_phi;
H_dop_dd = T_dd * H_dop;
H_phi_dd(:, length(errx_minus)-length(lambda)+1:end) = diag(ones(length(lambda),1));
H_dd = [H_rho_dd;H_phi_dd;H_dop_dd];
% measurement residual (res_rho - amb)
res_dd = res - H_dd(:,10:end)*errx_minus(10:end);

if p.i == 2519
    ii=1;
end
switch p.est_mode
    case p.ekf_est
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ekfUpdate(x_minus, ...
            p.error_state,cov_minus, res_dd, H_dd, R_dd);
        p.num_meas_used = num;
    case p.map_est
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        Rot_e2g = eye(length(errx_minus));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        H_dd = H_dd * Rot_e2g';
        cov_minus = Rot_e2g * cov_minus * Rot_e2g';
        [x_plus,dx_plus,cov_plus,infor_ned,p.augcost] = ...
            mapUpdate(p, ones(length(res_dd),1),x_minus,errx_minus,cov_minus,res_dd,...
            H_dd,R_dd,Rot_e2g);
        p.num_meas_used = length(res_dd)/3;
    case p.td_est
        [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        if (p.state_mode == p.pva_mode || p.state_mode == p.ins_mode) && flag_rapid == true
            % p.state_cov = zeros(size(p.state_cov));
            cov_minus(1:6,1:6) = cov_minus(1:6,1:6)+100^2*eye(6);
        end
        b = thresholdTest(p.td_lambda,cov_minus, res_dd, H_dd, R_dd);
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        Rot_e2g = eye(length(errx_minus));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        H_dd = H_dd * Rot_e2g';
        cov_minus = Rot_e2g * cov_minus * Rot_e2g';
        [x_plus,dx_plus,cov_plus,infor_ned,p.augcost] = ...
            mapUpdate(p, b, x_minus, errx_minus,cov_minus, res_dd, H_dd, R_dd, Rot_e2g);
        p.num_meas_used = sum(b);
    case p.raps_ned_est
        % Solve in NED frame
        % tic
        lla_deg = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla_deg);
        [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        if (p.state_mode == p.pva_mode || p.state_mode == p.ins_mode) && flag_rapid == true
            % p.state_cov = zeros(size(p.state_cov));
            p.state_cov(1:6,1:6) = p.state_cov(1:6,1:6)+50^2*eye(6);
        end
        Rot_e2g = eye(length(errx_minus));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        Ht = H_dd * Rot_e2g';
        xt_minus = zeros(length(errx_minus),1);
        cov_prior = Rot_e2g*cov_minus*Rot_e2g';
        if (p.state_mode == p.pva_mode || p.state_mode == p.ins_mode)
            num_constrain = 6;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec;...
                p.raps.velhor_cov_spec; p.raps.velhor_cov_spec;...
                p.raps.velver_cov_spec]);
            p_u = 50^2*eye(length(xt_minus));
            p_u(1:6,1:6) = cov_spec_ecef;
        elseif p.state_mode == p.pos_mode
            num_constrain = 3;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec]);
            p_u = [cov_spec_ecef, zeros(3, length(x_minus)-3);
                zeros(length(x_minus)-3, 3), p_clk];
        end
        J_l = p_u^(-1);
        % tic
        % mapRiskAverseNonBiSlackMaxJ (Non-binary DiagRAPS)
        % mapRiskAverseSlack (Binary DiagRAPS)
        % mapRiskAverseCvx (Binary DiagRAPS Globally Optimal Format)
        if cov_prior(1,1) < 150^2
            [flag,dx_plus,cov_ned,b,infor_ned,p.augcost,num_iter,constraint,p.pos_risk,p.raps_penalty] = ...
                mapRiskAverseNbSlackRtk(p, num_constrain,res_dd,Ht,cov_prior,R_dd,...
                diag(diag(J_l)),xt_minus);
            cov_plus = Rot_e2g' * cov_ned * Rot_e2g;
            x_plus = updateInsState(p, dx_plus, x_minus, Rot_e2g);
            p.num_meas_used = sum(b>0.001);
            b(b>0.01) = 1;
            b(b<=0.01) = 0;
            if p.state_mode ~= p.pos_mode
                p.raps_num_sat = sum(b(1:length(res_dd)/3));
            else
                p.raps_num_sat = sum(b);
            end
            p.raps_num_iter = num_iter;
            p.constraint = constraint;
            p.raps_flag = flag;
        else
            % if the prior is extremly unreliable, perform TD
            % Trust the measurements more than the prior
            b = thresholdTest(p.td_lambda,cov_prior, res_dd, H_dd, R_dd);
            [x_plus,dx_plus,cov_plus,infor_ned,p.augcost] = ...
                mapUpdate(p, b, x_minus, errx_minus,cov_prior, res_dd, H_dd, R_dd, Rot_e2g);
            p.num_meas_used = sum(b>0.001);
            p.raps_num_iter = NaN;
            p.constraint = NaN;
            p.raps_flag = NaN;
            p.raps_num_sat = sum(b);
        end
    otherwise
        error('Incorrect state estimation mode configuration');
end

comp_t = toc;
p.comp_t = comp_t;

num_user_states = p.modeToNumUserErrStates(p.state_mode);
p.error_state = dx_plus(1:num_user_states);
p.state0 = x_plus;
p.state_cov = cov_plus(1:num_user_states,1:num_user_states);
p.infor_ned = infor_ned(1:num_user_states,1:num_user_states);
% HH = [H_os(1:length(y_rho),1:3)];
% if p.est_mode ~= p.ekf_est && p.est_mode ~= p.map_est
%     b_rho = b(1:length(y_rho));
%     HH = diag(b_rho)*HH;
% end
% hSqrtInv = (HH'*HH)^(-1);
% p.GDOP = sqrt(trace(hSqrtInv));
p.GDOP = NaN;
% Set antenna position
if p.state_mode == p.ins_mode
    R_e2b_plus = convertQuatToRot(x_plus(7:10));
    lever_arm_e = R_e2b_plus' * lever_arm_b;  % Lever arm in ECEF frame.
    estState.pos = x_plus(1:3) + lever_arm_e;
else
    estState.pos = x_plus(1:3);
end
if p.state_mode == p.pva_mode || p.state_mode == p.ins_mode
    estState.vel = x_plus(4:6);
end
estState.clock_bias = NaN;
estState.clock_drift = NaN;

end


