function [p,estState,res] = stateUpdate_code_only(p, cpt, dt)
% Only use code measurement, no Doppler
tic
%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;
x_minus = p.state0;
num_user_errstates = p.modeToNumUserErrStates(p.state_mode);
[H_clk,x_clk] = formClkStatesAndH(cpt.num_sv);
if length(x_clk) + num_user_errstates + 1 ~= length(p.error_state)
    error('current No. of sys does not match the previous epoch');
end
%------------------%
y_rho = cpt.corr_range;
p.num_sats_window = [p.num_sats_window(2:length(p.num_sats_window)), length(y_rho)];

num = length(y_rho); % The number of measurement
H = zeros(num,num_user_errstates+length(x_clk)+1);
H(:,num_user_errstates+1:num_user_errstates+length(x_clk)) = H_clk;

r = zeros(length(y_rho),1);
lever_comp = zeros(length(y_rho),1); % lever arm compensation to the residual
s_pos_ecef = cpt.sat_pos_Rcorr;
if p.post_mode == 1 && p.IGS_enable == 1
    s_pos_ecef = cpt.sat_posprc_Rcorr;
end

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
s_v_ecef = cpt.sat_v_Rcorr;
for j=1:length(y_rho)
    % compute LOS for code measurement
    r(j)=norm(s_pos_ecef(:,j)-imu_pos);
    los_r_i = (imu_pos-s_pos_ecef(:,j))'/r(j);
    H(j,1:3)=los_r_i;
    H(j,7:9)=-los_r_i*vectorSkewSymMat(lever_arm_e);
    lever_comp(j) = los_r_i*lever_arm_e;
    %r(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),x_minus(1:3));
end
H_os = H;
[R, ~] = constructMeasNoise(p, cpt, dt); %cpt.elev, cpt.svprn_mark
% measurement residual
% x_minus represent the state where 4 states for the Quat.
res = y_rho - r - H_os(1:length(y_rho),16:end)*x_minus(17:end)-lever_comp;
[x_minus, p.state_cov, flag] = checkClockReset(p, x_minus, p.state_cov, ...
    num_user_errstates+1, res, cpt); 
if flag == true
    res = y_rho - r - H_os(1:length(y_rho),16:end)*x_minus(17:end)-lever_comp;
end
res_all=res;

% y - f(x0) = H (x - x0);
switch p.est_mode
    case p.ekf_est
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ekfUpdate(x_minus, p.error_state,p.state_cov, res_all, H_os, R);
    % case p.td_est
    %     [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res_all, H_os, R);
    case p.map_est_test
        cov_prior = p.state_cov;
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(ones(num,1),x_minus,p.error_state,cov_prior,res_all,H_os,R,eye(length(p.error_state)));
        p.num_meas_used = num;
    case p.map_est
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        Rot_e2g = eye(length(p.error_state));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        Rot_e2g(7:9,7:9) = R_e2g;
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g * p.state_cov * Rot_e2g';
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(ones(num,1),x_minus,p.error_state,cov_prior,res_all,H_os,R,Rot_e2g);
        p.num_meas_used = num;
    case p.td_est
        if p.i == 505
            oo=1;
        end
        [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        if (p.state_mode == p.pva_mode || p.state_mode == p.ins_mode) && flag_rapid == true
            p.state_cov = zeros(size(p.state_cov));
            p.state_cov(1:6,1:6) = diag(150^2*ones(1,6));
            p.state_cov(7:15,7:15) = diag(20^2*ones(1,9));
            p.state_cov(16:end-1,16:end-1) = diag(150^2*ones(1,length(p.error_state)-16));
            p.state_cov(end,end) = 20^2;
        end
        b = thresholdTest(p.td_lambda,p.state_cov, res_all, H_os, R);
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        Rot_e2g = eye(length(p.error_state));
        Rot_e2g(1:3,1:3) = R_e2g;
        Rot_e2g(4:6,4:6) = R_e2g;
        Rot_e2g(7:9,7:9) = R_e2g;
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g * p.state_cov * Rot_e2g';
        [x_plus,dx_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(b, x_minus, p.error_state,cov_prior, res_all, H_os, R, Rot_e2g);
        p.num_meas_used = sum(b);
    case p.raps_ned_est
        % Solve in NED frame
        % tic
        lla_deg = ecef2lla(x_minus(1:3)', 'WGS84');
        R_eg=computeRotForEcefToNed(lla_deg);
        if p.state_mode == p.pva_mode
            R_pva = [R_eg, zeros(3,6);
                zeros(3,3), R_eg, zeros(3,3);
                zeros(3,6), R_eg];
            [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
            % Rapid num of sat detected, may entering an open sky area,
            % reset prior covariance
            if flag_rapid == true
                p.state_cov = zeros(size(p.state_cov));
                p.state_cov(1:end-1,1:end-1) = diag(150^2*ones(1,length(x_minus)-1));
                p.state_cov(end,end) = 20^2;
            end
        else
            R_pva = R_eg;
        end
        Rot_e2g = [R_pva, zeros(num_user_errstates,length(x_minus)-num_user_errstates);
            zeros(length(x_minus)-num_user_errstates, num_user_errstates), eye(length(x_minus)-num_user_errstates)];
        Ht = H_os * Rot_e2g';
        xt_minus = Rot_e2g*(x_minus - x_minus);
        Pt_minus = Rot_e2g*p.state_cov*Rot_e2g';
        if p.state_mode == p.pva_mode
            num_constrain = 6;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec;...
                p.raps.velhor_cov_spec; p.raps.velhor_cov_spec;...
                p.raps.velver_cov_spec]);
            p_clk = diag([p.raps.va_cov_spec*ones(3,1);...
                p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(6, length(x_minus)-6);
                zeros(length(x_minus)-6, 6), p_clk];
        elseif p.state_mode == p.pos_mode
            num_constrain = 3;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec]);
            p_clk = diag([p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(3, length(x_minus)-3);
                zeros(length(x_minus)-3, 3), p_clk];
        end
        J_l = p_u^(-1);
        % tic
        % mapRiskAverseNonBiSlackMaxJ (Non-binary DiagRAPS)
        % mapRiskAverseSlack (Binary DiagRAPS)
        % mapRiskAverseCvx (Binary DiagRAPS Globally Optimal Format)
        [flag,x_ned,cov_ned,b,J_out,p.augcost,num_iter,constraint,p.pos_risk,p.raps_penalty] = ...
            mapRiskAverseNonBiSlackMaxJ(num_constrain,res_all,Ht,Pt_minus,R,...
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
            p.raps_num_sat = sum(b(1:length(y_rho)));
        else
            p.raps_num_sat = sum(b);
        end
        p.infor_ned = J_out;
        p.raps_num_iter = num_iter;
        p.constraint = constraint;
        p.raps_flag = flag;
        x_plus = Rot_e2g'*x_ned + x_minus;
    otherwise
        error('Incorrect state estimation mode configuration');
end

comp_t = toc;
p.comp_t = comp_t;

p.error_state = dx_plus;
p.state0 = x_plus;
p.state_cov = cov_plus;

HH = [H_os(1:length(y_rho),1:3),H_clk(1:length(y_rho),:)];
if p.est_mode ~= p.ekf_est && p.est_mode ~= p.map_est
    b_rho = b(1:length(y_rho));
    HH = diag(b_rho)*HH;
end
hSqrtInv = (HH'*HH)^(-1);
p.GDOP = sqrt(trace(hSqrtInv));

% Set antenna position
R_e2b_plus = convertQuatToRot(x_plus(7:10));
lever_arm_e = R_e2b_plus' * lever_arm_b;  % Lever arm in ECEF frame.
estState.pos = x_plus(1:3) + lever_arm_e;
if p.state_mode == p.pva_mode || p.state_mode == p.ins_mode
    estState.vel = x_plus(4:6);
end
estState.clock_bias = x_plus(num_user_errstates+1+1);
estState.clock_drift = x_plus(end);

clk_est = x_plus(num_user_errstates+1+1:end-1);
j = 1;
for i = 1:length(cpt.num_sv)
    if cpt.num_sv(i) == 0
        continue;
    end
    estState.clock_sys(i) = clk_est(j);
    j=j+1;
end

end


