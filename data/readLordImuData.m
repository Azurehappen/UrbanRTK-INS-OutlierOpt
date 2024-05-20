function [imu_data, imu_para] = readLordImuData(filename)
warning off
% Load data from file
data = readtable(filename, 'Delimiter', ',');

% Initialize arrays
imu_data.gps_week = data.x_Week';
imu_data.gps_sec = (data.wholesec+data.fracsec)';
imu_data.acc = [data.acc1';data.acc2';data.acc3']; % m/s^2
imu_data.gyro = [data.gyro1';data.gyro2';data.gyro3']; % rad/s

init_pos_ecef = [-741204.520;-5462376.740;3197933.705];
lla = ecef2lla(init_pos_ecef', 'WGS84');
lat = lla(1);
g_h = 9.7803267715*(1+0.001931851353*(sin(lat))^2) /...
    sqrt(1-0.0066943800229*(sin(lat))^2);  %navigation frame: NED, unit: m/s^2
imu_para.g_ned = [0;0;g_h];
imu_para.freq = 100; % Hz
% acc_offset = 0.04*g_h/1e3; % lifetime zero offset, unit: m/s^2
% imu_para.acc_lam = -0.01;
% fact_noise = 2*10^4;
% imu_para.acc_noise = (25*g_h/1e6)^2*fact_noise; % acc noise, 25 ug/sqrt(Hz) -> m/s^2/sqrt(Hz)
% % imu_para.acc_bias =  acc_offset^2*2*abs(imu_para.acc_lam)*100; % Angular random walk
% imu_para.acc_bias =  acc_offset^2*fact_noise;
% gyro_offset = 8*pi/180/3600; % Bias instability, unit: deg/h -> deg/sec
% imu_para.gyro_lam = -0.01;
% imu_para.gyro_noise = (0.005*pi/180)^2*fact_noise; % gyro noise, unit: (rad/sec/sqrt(Hz))^2
% % imu_para.gyro_bias =  gyro_offset^2*2*abs(imu_para.gyro_lam)*100; % Eqn. 4.102 in Farrell's book
% imu_para.gyro_bias = gyro_offset^2*fact_noise;

%% Some Paremeter informaiton from TEX-CUP author.
fact_noise = 50.00*10^4;
acc_ss_std = 0.5*g_h/1e3; % Accelerometer bias steady-state standard deviation: milli-g -> m/s^2
taua_ = 100; % Accelerometer bias time constant, in seconds
taug_ = 100; % Gyro bias time constant, in seconds
imu_para.acc_lam = -1/taua_;
imu_para.acc_noise = fact_noise*(0.1*g_h/1e3)^2; % acc noise, 0.1 mg/sqrt(Hz) -> (m/s^2/sqrt(Hz))^2
% imu_para.acc_bias =  acc_ss_std^2*2*abs(imu_para.acc_lam); % Angular random walk
imu_para.acc_bias =  fact_noise*acc_ss_std^2;
gyro_ss_std = 8*pi/180/3600; % Bias instability, unit: deg/h -> rad/sec
imu_para.gyro_lam = -1/taug_;
imu_para.gyro_noise = fact_noise*(0.01*pi/180)^2; % gyro noise, unit: (rad/sec/sqrt(Hz))^2
% imu_para.gyro_bias =  gyro_ss_std^2*2*abs(imu_para.gyro_lam); % Eqn. 4.102 in Farrell's book
imu_para.gyro_bias = fact_noise*gyro_ss_std^2;

R_e2n=computeRotForEcefToNed(lla');
init_vel = [-0.246;-0.527;-0.0204]; % NED
imu_para.init_yaw = -atan2(init_vel(2),init_vel(1));
R_b2n = eulerToRot(0,0,imu_para.init_yaw);
R_e2b = R_b2n'*R_e2n;
imu_para.init_quat = R2quat(R_e2b);
end