% Estimate Bosch BMI055 Correlation time

% Load data from file
data = readtable('univOfTexas\bosch_imu.log', 'Delimiter', ',');

% Initialize arrays
imu_data.gps_week = data.x_Week';
imu_data.gps_sec = (data.wholesec+data.fracsec)';
imu_data.acc = [data.acc1';data.acc2';data.acc3'];
imu_data.gyro = [data.gyro1';data.gyro2';data.gyro3'];

ind = find(imu_data.gps_sec>410538 & imu_data.gps_sec<410759);

acc_0 = imu_data.acc(:,ind);
gyro_0 = imu_data.gyro(:,ind);

acc_mean_x = mean(acc_0(1,:));
acc_mean_y = mean(acc_0(2,:));
acc_mean_z = mean(acc_0(3,:));

%[avar,tau] = allanvar(acc_0(1,:),'octave',150);

acc_x_nm = acc_0(1,:)-acc_mean_x;
%%
Fs = 150;  % Sampling frequency in Hz

% Design a high-pass Butterworth filter
fc = 30;  % Cut-off frequency
n = 2;  % Order of the filter
[b, a] = butter(n, fc/(Fs/2), 'high');

% Apply the filter
filtered_data = filter(b, a, acc_x_nm);

% FFT of the filtered data
Y = fft(filtered_data);
N = length(Y);
T = 1/Fs;
f = Fs*(0:(N/2))/N;

% Compute the magnitude
P2 = abs(Y);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Plotting
figure(3)
plot(f, P1);
title('Filtered Frequency Domain (Above 30 Hz)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%%
[acf,lags] = autocorr(acc_0(1,:)-acc_mean_x);

figure(3)
plot(lags,acf,'.')

figure(1)
plot(acc_0(1,:)-acc_mean_x,'.')
hold on
plot(acc_0(2,:)-acc_mean_y,'.')
plot(acc_0(3,:)-acc_mean_z,'.')
grid on

figure(2);
loglog(tau,avar)
xlabel('\tau')
ylabel('\sigma^2(\tau)')
title('Allan Variance')
grid on