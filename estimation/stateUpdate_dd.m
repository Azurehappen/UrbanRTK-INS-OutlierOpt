function [p,estState,res] = stateUpdate_dd(p, cpt, dt)
tic
%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;
% Initialize
x_minus = p.state0;
num_user_errstates = p.modeToNumUserErrStates(p.state_mode);
%------------------%
y_rho = cpt.corr_range;
p.num_sats_window = [p.num_sats_window(2:length(p.num_sats_window)), length(y_rho)];
y_dop = [];
if p.state_mode == p.pva_mode || p.state_mode == p.ins_mode
    y_dop = cpt.doppler;
end
num = length(y_rho) + length(y_dop); % The number of measurement
H = zeros(num,num_user_errstates);
if p.state_mode == p.pva_mode  || p.state_mode == p.ins_mode
    gyro_noise2dop = ones(length(y_dop),1);
    res_v = zeros(length(y_dop),1);
end
r = zeros(length(y_rho),1);
lever_comp = zeros(length(y_rho),1); % lever arm compensation to the residual
s_pos_ecef = cpt.sat_pos_Rcorr;
if p.post_mode == 1 && p.IGS_enable == 1
    s_pos_ecef = cpt.sat_posprc_Rcorr;
end
s_v_ecef = cpt.sat_v_Rcorr;

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
for j=1:length(y_rho)
    % compute LOS for code measurement
    r(j)=norm(s_pos_ecef(:,j)-imu_pos);
    los_r_i = (imu_pos-s_pos_ecef(:,j))'/r(j);
    % compute LOS for Doppler
    range_r = norm(s_pos_ecef(:,j)-ant_pos);
    los_r = (ant_pos-s_pos_ecef(:,j))'/range_r;
    H(j,1:3)=los_r_i;
    H(j,7:9)=-los_r_i*vectorSkewSymMat(lever_arm_e);
    lever_comp(j) = los_r_i*lever_arm_e;
    if p.state_mode == p.pva_mode
        H(j+length(y_rho),4:6) = los_r;
        res_v(j) = y_dop(j) - los_r*(x_minus(4:6) - s_v_ecef(:,j));
    elseif p.state_mode == p.ins_mode
        H(j+length(y_rho),4:6) = los_r;
        H(j+length(y_rho),7:9) = -los_r*vectorSkewSymMat(Omg_lever);
        H(j+length(y_rho),13:15) = los_r*R_e2b_hat'*vectorSkewSymMat(lever_arm_b);
        gyro_noise2dop(j) = H(j+length(y_rho),13:15)*(p.imu_para.gyro_noise*eye(3,3))*H(j+length(y_rho),13:15)';
        res_v(j) = y_dop(j) - los_r*(ant_vel - s_v_ecef(:,j));
    end
    %r(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),x_minus(1:3));
end
H_os = H;
[R, ~] = constructMeasNoise(p, cpt, dt); %cpt.elev, cpt.svprn_mark
% measurement residual
% x_minus represent the state where 4 states for the Quat.
res = y_rho - r -lever_comp;

if p.state_mode == p.pva_mode || p.state_mode == p.ins_mode
    % Rdop = p.sigmaSquare_dop*eye(length(y_dop));
    Rdop = diag(p.sigmaSquare_dop+gyro_noise2dop);
    % Rdop_diag = p.sigmaSquare_dop+(0.5*sqrt(p.sigmaSquare_dop)./sin(cpt.elev)).^2+gyro_noise2dop;
    % Rdop = diag(Rdop_diag);
    R = [R,zeros(length(y_rho),length(y_dop));
        zeros(length(y_dop),length(y_rho)),Rdop];
    res_all=[res;res_v];
else
    res_all=res;
end

[res_dd, H_dd, R_dd] = doubleDiff(cpt, res_all, H_os, R, 5);
R_dd = diag(diag(R_dd));
num = length(res_dd);
if isempty(res_dd)
    estState.pos = [];
    if p.save_ins_pos == true
        R_e2b_plus = convertQuatToRot(p.state0(7:10));
        lever_arm_e = R_e2b_plus' * lever_arm_b;  % Lever arm in ECEF frame.
        estState.pos = p.state0(1:3) + lever_arm_e;
        estState.vel = p.state0(4:6);
        estState.clock_bias = NaN;
        estState.clock_drift = NaN;
    end
    p.state_cov(1:6,1:6) = p.state_cov(1:6,1:6)+100^2*eye(6);
    return;
end
% y - f(x0) = H (x - x0);
switch p.est_mode
    case p.ekf_est
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ekfUpdate(x_minus, ...
            p.error_state,p.state_cov, res_dd, H_dd, R_dd);
        p.num_meas_used = num;
    case p.map_est
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        Rot_e2g = eye(length(p.error_state));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        H_dd = H_dd * Rot_e2g';
        cov_prior = Rot_e2g * p.state_cov * Rot_e2g';
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(ones(num,1),x_minus,p.error_state,cov_prior,res_dd,...
            H_dd,R_dd,Rot_e2g);
        p.num_meas_used = num;
    case p.td_est
        [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        if (p.state_mode == p.pva_mode || p.state_mode == p.ins_mode) && flag_rapid == true
            % p.state_cov = zeros(size(p.state_cov));
            p.state_cov(1:6,1:6) = p.state_cov(1:6,1:6)+100^2*eye(6);
        end
        b = thresholdTest(p.td_lambda,p.state_cov, res_dd, H_dd, R_dd);
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        Rot_e2g = eye(length(p.error_state));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        H_dd = H_dd * Rot_e2g';
        cov_prior = Rot_e2g * p.state_cov * Rot_e2g';
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(b, x_minus, p.error_state,cov_prior, res_dd, H_dd, R_dd, Rot_e2g);
        p.num_meas_used = sum(b);
    case p.raps_ned_est
        % Solve in NED frame
        % tic
        lla_deg = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla_deg);
        [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        if (p.state_mode == p.pva_mode || p.state_mode == p.ins_mode) && flag_rapid == true
            % p.state_cov = zeros(size(p.state_cov));
            p.state_cov(1:6,1:6) = p.state_cov(1:6,1:6)+100^2*eye(6);
        end
        Rot_e2g = eye(length(p.error_state));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        Ht = H_dd * Rot_e2g';
        xt_minus = zeros(length(p.error_state),1);
        cov_prior = Rot_e2g*p.state_cov*Rot_e2g';
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
        % mapRiskAverseNonBinarySdp (Non-binary Full- RAPS)
        [flag,dx_plus,cov_ned,b,J_out,p.augcost,num_iter,constraint,p.pos_risk,p.raps_penalty] = ...
            mapRiskAverseNonBiSlackMaxJ(num_constrain,res_dd,Ht,cov_prior,R_dd,...
            diag(diag(J_l)),xt_minus);
        % comp_t = toc;
        % p.comp_t = comp_t;
        % tic
        % [~,~,~,~,~,p.augcost_bcd,~,~,p.pos_risk_bcd] = ...
        %     mapRiskAverseSlack(num_constrain,res_all,Ht,Pt_minus,R,...
        %     diag(diag(J_l)),xt_minus,p.td_lambda,length(y_rho));
        % comp_t = toc;
        % p.comp_t_bcd = comp_t;
        cov_plus = Rot_e2g' * cov_ned * Rot_e2g;
        p.num_meas_used = sum(b>0.001);
        b(b>0.01) = 1;
        b(b<=0.01) = 0;
        if p.state_mode == p.pva_mode
            p.raps_num_sat = sum(b(1:length(res_dd)/2));
        else
            p.raps_num_sat = sum(b);
        end
        p.infor_ned = J_out;
        p.raps_num_iter = num_iter;
        p.constraint = constraint;
        p.raps_flag = flag;
        x_plus = updateInsState(dx_plus, x_minus, Rot_e2g);
    otherwise
        error('Incorrect state estimation mode configuration');
end

comp_t = toc;
p.comp_t = comp_t;

p.error_state = dx_plus;
p.state0 = x_plus;
p.state_cov = cov_plus;

% HH = [H_os(1:length(y_rho),1:3)];
% if p.est_mode ~= p.ekf_est && p.est_mode ~= p.map_est
%     b_rho = b(1:length(y_rho));
%     HH = diag(b_rho)*HH;
% end
% hSqrtInv = (HH'*HH)^(-1);
% p.GDOP = sqrt(trace(hSqrtInv));
p.GDOP = NaN;
% Set antenna position
R_e2b_plus = convertQuatToRot(x_plus(7:10));
lever_arm_e = R_e2b_plus' * lever_arm_b;  % Lever arm in ECEF frame.
estState.pos = x_plus(1:3) + lever_arm_e;
if p.state_mode == p.pva_mode || p.state_mode == p.ins_mode
    estState.vel = x_plus(4:6);
end
estState.clock_bias = NaN;
estState.clock_drift = NaN;

end


