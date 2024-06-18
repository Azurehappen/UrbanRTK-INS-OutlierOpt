function [imu_data, imu_para] = readBoschImuData(filename)

% Load data from file
data = readtable(filename, 'Delimiter', ',');

% Initialize arrays
imu_data.gps_week = data.x_Week';
imu_data.gps_sec = (data.wholesec+data.fracsec)';
imu_data.acc = [data.acc1';data.acc2';data.acc3'];
imu_data.gyro = [data.gyro1';data.gyro2';data.gyro3'];
imu_para.R_i2b = diag([1,-1,-1]); % The IMU was place at inverse direction of z-axis.

init_pos_ecef = [-741204.520;-5462376.740;3197933.705];
lla = ecef2lla(init_pos_ecef', 'WGS84');
lat = lla(1);
g_h = 9.7803267715*(1+0.001931851353*(sin(lat))^2) /...
    sqrt(1-0.0066943800229*(sin(lat))^2);  %navigation frame: NED, unit: m/s^2
imu_para.g_ned = [0;0;g_h];
fact_noise = 1;
imu_para.freq = 150; % Hz
acc_ss_std = 10*g_h/1e3; % lifetime zero offset, unit: m/s^2
taua_ = 10; % Accelerometer bias time constant, in seconds
taug_ = 200; % Gyro bias time constant, in seconds
imu_para.acc_lam = -1/taua_;
imu_para.acc_noise = fact_noise*(0.1*g_h^2/1e6); % acc noise, mg^2/Hz -> (m/s^2/sqrt(Hz))^2
% imu_para.acc_bias =  acc_offset^2*2*abs(imu_para.acc_lam); % Angular random walk
imu_para.acc_bias =  fact_noise*acc_ss_std^2;
gyro_ss_std = 200*pi/180/3600; % lifetime zero offset, unit: deg/h -> rad/sec
imu_para.gyro_lam = -1/taug_;
imu_para.gyro_noise = fact_noise*0.02*(pi/180)^2; % gyro noise, unit: (deg/s)^2/Hz
% imu_para.gyro_bias =  gyro_ss_std^2*2*abs(imu_para.gyro_lam); % Eqn. 4.102 in Farrell's book
imu_para.gyro_bias = fact_noise*gyro_ss_std^2;

R_e2n=computeRotForEcefToNed(lla');
init_vel = [-0.246;-0.527;-0.0204]; % NED
imu_para.init_yaw = -atan2(init_vel(2),init_vel(1));
R_b2n = eulerToRot(0,0,imu_para.init_yaw);
R_e2b = R_b2n'*R_e2n;
imu_para.init_quat = R2quat(R_e2b);

end