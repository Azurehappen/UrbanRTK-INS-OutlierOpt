function [state, error_state, cov] = obtainInitInsStateAndCov(p, estState, tr)
% Initialize INS states

% IMU is stationary at the beginning
% Find acc. meas. near tr
ind = find(p.imu_data.gps_sec>=tr-1 & p.imu_data.gps_sec<=tr+1);
% Obtain acc. meas.
acc = p.imu_data.acc(:,ind);
acc_mean = -mean(acc,2);
% Euler angles from body to nav (NED)
phi = atan2(acc_mean(2),acc_mean(3)); % Initial roll, Eqn. 10.14, Farrell Book 
theta = atan2(acc_mean(2),acc_mean(3)); % Initial pitch
% No magnetic field data
psi = 0; % Initial yaw, should come with a large error covariance

estState.lla_deg = ecef2lla(estState.pos', 'WGS84')';
R_e2n = computeRotForEcefToNed(estState.lla_deg);
R_b2n = computeDirectionCos(phi,theta,psi);
% quat = directionCosToQuat(R_e2n'*R_b2n); % Intial quaternion
quat = p.imu_para.init_quat;

% State: position, velocity, quaternion, acc_bias, gyro_bias
% clock biases (m), clock drift (m/s)
% Error state: position, velocity, attitude errors, acc_bias, gyro_bias
% clock biases (m), clock drift (m/s)
state = [estState.pos;estState.vel;quat;0.0*ones(6, 1)];
% Error covariance
cov_diag = [p.ekf_para.q_pos*ones(1,3),p.ekf_para.q_vel*ones(1,3),...
            deg2rad(10)^2,deg2rad(10)^2,deg2rad(10)^2,2.5e-3*ones(1, 3),4e-6*ones(1, 3)];

if p.double_diff == false && ~isnan(estState.clock_sys(p.gps.sys_num)) 
    state = [state; estState.clock_sys(p.gps.sys_num)];
    cov_diag = [cov_diag, 40^2];
end
if p.double_diff == false && ~isnan(estState.clock_sys(p.glo.sys_num))
    state = [state; estState.clock_sys(p.glo.sys_num)];
    cov_diag = [cov_diag, 40^2];
end
if p.double_diff == false && ~isnan(estState.clock_sys(p.gal.sys_num))
    state = [state; estState.clock_sys(p.gal.sys_num)];
    cov_diag = [cov_diag, 40^2];
end
if p.double_diff == false && ~isnan(estState.clock_sys(p.bds.sys_num))
    state = [state; estState.clock_sys(p.bds.sys_num)];
    cov_diag = [cov_diag, 40^2];
end

if p.double_diff == false
    state = [state; estState.clock_drift];
    cov_diag = [cov_diag, 5^2];
end
error_state = zeros(length(state)-1,1); % In error state, 3 attitude errors.
cov = diag(cov_diag);
