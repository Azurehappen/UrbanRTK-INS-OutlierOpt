data = readtable('lord_imu.log', 'Delimiter', ',');

acc = [data.acc1';data.acc2';data.acc3'];
acc = acc(:,1:200*200);
gyro = [data.gyro1';data.gyro2';data.gyro3'];
gyro = gyro(:,1:200*200);

acc1(1,:) = acc(1,:)-mean(acc(1,:));
acc1(2,:) = acc(2,:)-mean(acc(2,:));
acc1(3,:) = acc(3,:)-mean(acc(3,:));

acc_sum = zeros(3,size(acc1,2));
for i=2:size(acc1,2)
    acc_sum(:,i)=acc1(:,i-1)*0.01+acc_sum(:,i-1);
end

gyro1(1,:) = gyro(1,:)-mean(gyro(1,:));
gyro1(2,:) = gyro(2,:)-mean(gyro(2,:));
gyro1(3,:) = gyro(3,:)-mean(gyro(3,:));

figure(1)
plot(acc_sum(1,:),'.')
hold on


figure(2)
plot(gyro1(1,:),'.')
hold on

