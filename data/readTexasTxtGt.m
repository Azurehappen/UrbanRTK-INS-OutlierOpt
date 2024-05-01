function [pos_ecef, vel_ned, sow, dtime] = readTexasTxtGt(file_name, lever_arm)

% Load data from file
data = readtable(file_name, 'Delimiter', '\t');

% Initialize arrays
sow = zeros(1, size(data,1));
pos_ecef = zeros(3, size(data,1));
vel_ned = zeros(3, size(data,1));
dtime = [];
% Convert date and time to gps seconds of week
% The time in the GT is UTC
for i = 1:size(data,1)
    date = datetime(data.x_date(i), 'InputFormat', 'yyyy/MM/dd');
    date_time = date + data.time(i) + seconds(18);
    date_vector = datevec(date_time);
    [~, ~, sow(i)] = date2gnsst(date_vector);
    dtime = [dtime,datetime(date_vector)];
end

% Convert latitude, longitude, and ellipsoidHeight to ECEF coordinates
for i = 1:size(data,1)
    lla_deg = [data.latitude(i), data.longitude(i), data.ellipsoidHeight(i)];
    pos_imu_e = lla2ecef(lla_deg)';
    % rot_body2enu = euler_R_body2enu(data.heading(i), data.pitch(i), data.roll(i));
    % rot_ecef2enu = computeRotForEcefToEnu(lla_deg);
    % lever_arm_ecef = rot_ecef2enu'*rot_body2enu*lever_arm;
    pos_ecef(:,i) = pos_imu_e; %+ lever_arm_ecef;
    vel_ned(:,i) = [data.speedNorth(i);data.speedEast(i);data.speedVertical(i)];
end

% Rotation definitions from Texas dataset author
    function r = euler_R_body2enu(h_deg, p_deg, r_deg)
        heading = deg2rad(h_deg);
        pitch = deg2rad(p_deg);
        roll = deg2rad(r_deg);
        R1 = [cos(heading),sin(heading),0;
              -sin(heading),cos(heading),0;
              0,0,1];
        R2 = [1,0,0;
              0,cos(pitch),-sin(pitch);
              0,sin(pitch),cos(pitch)];
        R3 = [cos(roll),0,sin(roll);
              0,1,0;
              -sin(roll),0,cos(roll)];
        r = R1*R2*R3;
    end
end