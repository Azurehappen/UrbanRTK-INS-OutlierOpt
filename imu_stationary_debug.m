% load('out_debug.mat')
grd_data = readtable('data/univOfTexas/ground_truth.log', 'Delimiter', '\t');
grd = p.Grdpos.pos(:,1:100);
grd = mean(grd,2);

lla_gt_deg = ecef2lla(grd', 'WGS84');
wgs84 = wgs84Ellipsoid('meter');
R_e2g=computeRotForEcefToNed(lla_gt_deg');

ins_state = output.ins_est.imu_state;
ins_cov_ned = output.ins_est.imu_cov_ned;
ins_cov_xyz = output.ins_est.imu_cov_xyz;
meas_ind = output.ins_est.meas_t_ind;

[n,m] = size(ins_state);
xyz_err = ins_state(1:3,:) - grd;
xyz_err_std = sqrt(ins_cov_xyz(1:3,:));
ned_err = zeros(3,m);
ned_err_std = zeros(3,m);
vel_ned = zeros(3,m);
vel_ned_std = zeros(3,m);
vel_xyz = ins_state(4:6,:);
vel_xyz_std = sqrt(ins_cov_xyz(4:6,:));

acc_bias = ins_state(11:13,:);
acc_std = sqrt(ins_cov_xyz(10:12,:));
gyro_bias = ins_state(14:16,:);
gyro_std = sqrt(ins_cov_xyz(13:15,:));

quat_e2b = ins_state(7:10,:);
euler = zeros(3,m);

for i=1:m
    vel_ned(:,i) = R_e2g*ins_state(4:6,i);
    vel_ned_std(:,i) = sqrt(ins_cov_ned(4:6,i));
    [xNorth,yEast,zDown] = ecef2ned(ins_state(1,i),ins_state(2,i),ins_state(3,i),...
        lla_gt_deg(1),lla_gt_deg(2),lla_gt_deg(3),wgs84);
    ned_err_tp = [xNorth;yEast;zDown];
    ned_err(:,i) = ned_err_tp;
    ned_err_std(:,i) = sqrt(ins_cov_ned(1:3,i));
    R_e2b = convertQuatToRot(quat_e2b(:,i));
    R_n2b = R_e2b*R_e2g';
    [roll_phi_rad,pitch_theta_rad,yaw_psi_rad] = dcmToEuler(R_n2b');
    euler(:,i) = rad2deg([yaw_psi_rad;pitch_theta_rad;roll_phi_rad]);
end

% Time vector for plotting
% If you have specific time data, use that instead
time = (1:m)*0.01; % Example: 1 to m time steps

% Creating the 3x2 subplot
figure(10);

% Defining subplots for Position and Velocity Error in XYX directions with STD and measurement lines
for i = 1:6
    subplot(3,2,i);
    hold on;
    grid on;
    if mod(i,2) == 1 % Odd indices for position error and STD
        dir_index = (i + 1) / 2;
        plot(time, xyz_err(dir_index,:), 'b.', time, xyz_err_std(dir_index,:), 'r.', time, -xyz_err_std(dir_index,:), 'r.');
        ylim([-2,2])
    else % Even indices for velocity error and STD
        dir_index = i / 2;
        plot(time, vel_xyz(dir_index,:), 'b.', time, vel_xyz_std(dir_index,:), 'r.', time, -vel_xyz_std(dir_index,:), 'r.');
        ylim([-1,1])
    end
    
    % Adding vertical lines for measurement indices
    for ind = 1:length(meas_ind)
        xline = [meas_ind(ind), meas_ind(ind)];
        ylimits = ylim; % Getting current y-axis limits
        line(xline, [ylimits(1), ylimits(2)], 'Color', 'k', 'LineStyle', '-');
    end
    
    hold off;
    xlabel('Time, sec');
    ylabel('Error/STD');
    
    % Setting titles based on the subplot
    titles = {'Position Error in X and STD', 'Velocity Error in X and STD', ...
              'Position Error in Y and STD', 'Velocity Error in Y and STD', ...
              'Position Error in Z and STD', 'Velocity Error in Z and STD'};
    title(titles{i});
end


figure(12)
for i = 1:6
    subplot(3,2,i);
    hold on;
    grid on;
    if mod(i,2) == 1 % Odd indices for position error and STD
        dir_index = (i + 1) / 2;
        plot(time, acc_bias(dir_index,:), 'b.', time, acc_std(dir_index,:), 'r.', time, -acc_std(dir_index,:), 'r.');
        ylim([-2,2])
    else % Even indices for velocity error and STD
        dir_index = i / 2;
        plot(time, gyro_bias(dir_index,:), 'b.', time, gyro_std(dir_index,:), 'r.', time, -gyro_std(dir_index,:), 'r.');
        ylim([-1,1])
    end
    
    % Adding vertical lines for measurement indices
    for ind = 1:length(meas_ind)
        xline = [meas_ind(ind), meas_ind(ind)];
        ylimits = ylim; % Getting current y-axis limits
        line(xline, [ylimits(1), ylimits(2)], 'Color', 'k', 'LineStyle', '-');
    end
    
    hold off;
    xlabel('Time, sec');
    ylabel('Error/STD');
    
    % Setting titles based on the subplot
    titles = {'Acc bias in X and STD', 'Gyro bias in X and STD', ...
              'Acc bias in Y and STD', 'Gyro bias in Y and STD', ...
              'Acc bias in Z and STD', 'Gyro bias in Z and STD'};
    title(titles{i});
end

figure(13)
subplot(311);
plot(time,euler(1,:), 'b.')
grid on;
hold on;
for ind = 1:length(meas_ind)
    xline = [meas_ind(ind), meas_ind(ind)]*0.01;
    ylimits = ylim; % Getting current y-axis limits
    line(xline, [ylimits(1), ylimits(2)], 'Color', 'k', 'LineStyle', '-');
end
hold off;
xlabel('Time, sec');
ylabel('Euler, yaw');
subplot(312);
plot(time,euler(2,:), 'b.')
grid on;
hold on;
for ind = 1:length(meas_ind)
    xline = [meas_ind(ind), meas_ind(ind)]*0.01;
    ylimits = ylim; % Getting current y-axis limits
    line(xline, [ylimits(1), ylimits(2)], 'Color', 'k', 'LineStyle', '-');
end
hold off;
xlabel('Time, sec');
ylabel('Euler, pitch');
subplot(313);
plot(time,euler(3,:), 'b.')
grid on;
hold on;
for ind = 1:length(meas_ind)
    xline = [meas_ind(ind), meas_ind(ind)]*0.01;
    ylimits = ylim; % Getting current y-axis limits
    line(xline, [ylimits(1), ylimits(2)], 'Color', 'k', 'LineStyle', '-');
end
hold off;
xlabel('Time, sec');
ylabel('Euler, roll');
%%

% Creating the 3x2 subplot
figure(11);

% Defining subplots for Position and Velocity Error in NED directions with STD and measurement lines
for i = 1:6
    subplot(3,2,i);
    hold on;
    grid on;
    if mod(i,2) == 1 % Odd indices for position error and STD
        dir_index = (i + 1) / 2;
        plot(time, ned_err(dir_index,:), 'b.', time, ned_err_std(dir_index,:), 'r.', time, -ned_err_std(dir_index,:), 'r.');
        ylim([-2,2])
    else % Even indices for velocity error and STD
        dir_index = i / 2;
        plot(time, vel_ned(dir_index,:), 'b.', time, vel_ned_std(dir_index,:), 'r.', time, -vel_ned_std(dir_index,:), 'r.');
        ylim([-1,1])
    end
    
    % Adding vertical lines for measurement indices
    for ind = 1:length(meas_ind)
        xline = [meas_ind(ind), meas_ind(ind)];
        ylimits = ylim; % Getting current y-axis limits
        line(xline, [ylimits(1), ylimits(2)], 'Color', 'k', 'LineStyle', '-');
    end
    
    hold off;
    xlabel('Time, sec');
    ylabel('Error/STD');
    
    % Setting titles based on the subplot
    titles = {'Position Error in North and STD', 'Velocity Error in North and STD', ...
              'Position Error in East and STD', 'Velocity Error in East and STD', ...
              'Position Error in Down and STD', 'Velocity Error in Down and STD'};
    title(titles{i});
end

