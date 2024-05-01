function [p,ins_est] = computeGnssEpochPrior(p, last_t, next_t, imu_ind,ins_est)

I3x3 = eye(3);
Z3x3 = zeros(3,3);
Z3x1 = zeros(3,1);

acc_meas = p.imu_data.acc(:,imu_ind);
gyro_meas = p.imu_data.gyro(:,imu_ind);

time = [p.imu_data.gps_sec(imu_ind),next_t];

if (last_t ~= time(1))
    ind = [imu_ind(1)-1, imu_ind];
    acc_meas = p.imu_data.acc(:,ind);
    gyro_meas = p.imu_data.gyro(:,ind);
    time = [last_t,time];
end

num_user_errstates = p.modeToNumUserErrStates(p.state_mode);
num_clk = length(p.error_state) - num_user_errstates - 1;

Z3xCLK = zeros(3, num_clk);

state = p.state0;
error_state = p.error_state;
error_cov = p.state_cov;

% ins_est.imu_state = [ins_est.imu_state,state];
% lla = ecef2lla([-741204.316823254,-5462377.19800711,3197934.07609661], 'WGS84');
% R_e2g=computeRotForEcefToNed(lla');
% Rot_e2g = eye(length(error_state));
% Rot_e2g(1:3,1:3) = R_e2g;
% Rot_e2g(4:6,4:6) = R_e2g;
% ins_est.imu_cov_xyz = [ins_est.imu_cov_ned,diag(error_cov)];
% cov_ned = Rot_e2g * error_cov * Rot_e2g';
% ins_est.imu_cov_ned = [ins_est.imu_cov_ned,diag(cov_ned)];
% ins_est.meas_t_ind = [ins_est.meas_t_ind, size(ins_est.imu_state,2)];
for i = 1:length(time)-1
    dt = time(i+1) - time(i);
    
    % Obtain states
    p_e   = state(1:3,1); % position in ECEF
    v_e   = state(4:6,1);
    q_e2b = state(7:10,1);
    bias_a = state(11:13,1);
    bias_g = state(14:16,1);
    clk = state(17:end-1,1);
    clk_drift = state(end,1);

    lla = ecef2lla(p_e','WGS84'); % compute lat, lng, h from p_e
    R_e2n = computeRotForEcefToNed(lla); % rotation ECEF to NED
    R_n2e = R_e2n';
    R_e2b = convertQuatToRot(q_e2b);
    R_b2e = R_e2b';
    
    F_b = acc_meas(:,i)-bias_a;     % specific force in body
    w_ib_b = gyro_meas(:,i)-bias_g; % angular rate in body

    % integrate position state
    state(1:3,1) = state(1:3,1) + v_e*dt;

    % integrate velocity
    F_e = R_b2e*F_b;                            % specific force in ecef
    g_e = R_n2e*p.imu_para.g_ned;% - p.Omega_iee_mat*p.Omega_iee_mat*p_e;  % gravity in ecef, eqn 11.14
    acc_e = F_e + g_e - 2*p.Omega_iee_mat*v_e; % eqn 11.32
    state(4:6,1) = v_e + dt*acc_e;   % eqn 11.32
    % acc_used = acc_e - 2*p.Omega_ie_e*v_e;
    
    % integrate e2b quaternion
    w_be_b = R_e2b*p.omg_iee_vec - w_ib_b;       % eqn 10.24
    W   = w_be_b*dt/2;                          % Assumes constant angular rate over
    w   = norm(W);                              % the period of integration.
    Wc  = vectorSkewSymMat(W);
    W_mat = [0 -W'
        W  Wc];
    if w == 0
        sinwow = 1;
    else
        sinwow = sin(w)/w;
    end
    state(7:10,1) = (cos(w)*eye(4) + W_mat*sinwow)*q_e2b;  % eqn. D.36

    % integrate receiver clock
    state(17:end-1,1) = state(17:end-1,1) + clk_drift*dt;

    % The following section compute error state model
    % position state row
    Fpv = I3x3;
    % velocity state row
    ptp = p_e'*p_e;
    np3 = ptp^(1.5);
    np5 = ptp^(2.5);
    Fvp = -p.GM*(I3x3/np3 - 3*(p_e*p_e')/np5);  % see eqn. (11.16) in Aided Navigation
    Fvv = -2*p.Omega_iee_mat;
    Fvr = -vectorSkewSymMat(F_e);
    Fva = -R_b2e;
    % attitude error state row
    Frr = -p.Omega_iee_mat;
    Frg = -R_b2e;
    % Faa = p.imu_para.acc_lam*I3x3;
    % Fgg = p.imu_para.gyro_lam*I3x3;
    Faa = 0*I3x3;
    Fgg = 0*I3x3;
    
    if p.double_diff == false
        F = [ Z3x3 Fpv  Z3x3 Z3x3 Z3x3 Z3xCLK Z3x1;  % p - pos error in ecef
              Fvp  Fvv  Fvr  Fva  Z3x3 Z3xCLK Z3x1;  % v - velocity error in ecef
              Z3x3 Z3x3 Frr  Z3x3 Frg  Z3xCLK Z3x1;  % r - rho, attitude error in ecef
              Z3x3 Z3x3 Z3x3 Faa  Z3x3 Z3xCLK Z3x1;  % a - accelerometer bias error in body
              Z3x3 Z3x3 Z3x3 Z3x3 Fgg  Z3xCLK Z3x1; % g - gyro bias error in body
              zeros(num_clk, length(p.error_state)-1), ones(num_clk,1); % clks
              zeros(1, length(p.error_state))]; % clk drfit
    
        G = [Z3x3 Z3x3 Z3x3 Z3x3 Z3x1;
             -Fva Z3x3 Z3x3 Z3x3 Z3x1;
             Z3x3 -Frg Z3x3 Z3x3 Z3x1;
             Z3x3 Z3x3 I3x3 Z3x3 Z3x1;
             Z3x3 Z3x3 Z3x3 I3x3 Z3x1;
             Z3xCLK' Z3xCLK' Z3xCLK' Z3xCLK' zeros(num_clk,1);
             zeros(1, 3*4),1];
        Qdiag = [p.imu_para.acc_noise*ones(3,1);
            p.imu_para.gyro_noise*ones(3,1);
            p.imu_para.acc_bias*ones(3,1);
            p.imu_para.gyro_bias*ones(3,1);
            p.ekf_para.q_clkDrift];
    else
        F = [ Z3x3 Fpv  Z3x3 Z3x3 Z3x3;  % p - pos error in ecef
              Fvp  Fvv  Fvr  Fva  Z3x3;  % v - velocity error in ecef
              Z3x3 Z3x3 Frr  Z3x3 Frg ;  % r - rho, attitude error in ecef
              Z3x3 Z3x3 Z3x3 Faa  Z3x3;  % a - accelerometer bias error in body
              Z3x3 Z3x3 Z3x3 Z3x3 Fgg ]; % g - gyro bias error in body
    
        G = [Z3x3 Z3x3 Z3x3 Z3x3;
             -Fva Z3x3 Z3x3 Z3x3;
             Z3x3 -Frg Z3x3 Z3x3;
             Z3x3 Z3x3 I3x3 Z3x3;
             Z3x3 Z3x3 Z3x3 I3x3];
        Qdiag = [p.imu_para.acc_noise*ones(3,1);
            p.imu_para.gyro_noise*ones(3,1);
            p.imu_para.acc_bias*ones(3,1);
            p.imu_para.gyro_bias*ones(3,1)];
    end
    [Phi,Qd] = computeDiscreteModel(F, G, diag(Qdiag), dt);
    
    error_state = Phi*error_state;
    error_cov = Phi*error_cov*Phi'+Qd;

    % ins_est.imu_state = [ins_est.imu_state,state];
    % ins_est.imu_cov_xyz = [ins_est.imu_cov_ned,diag(error_cov)];
    % cov_ned = Rot_e2g * error_cov * Rot_e2g';
    % ins_est.imu_cov_ned = [ins_est.imu_cov_ned,diag(cov_ned)];
end
% 
% ins_est.imu_state(:,end) = [];
% ins_est.imu_cov_ned(:,end) = [];
% ins_est.imu_cov_xyz(:,end) = [];

p.state0 = state;
p.error_state = error_state;
p.state_cov = error_cov;
p.w_be_b = w_be_b;