% figure
% subplot(211)
% plot(p.t,output.raps_constraint(1,:))
% hold on
% plot(p.t,output.raps_constraint(2,:))
% hold on
% plot(p.t,output.raps_constraint(3,:))
% grid on
% legend('North', 'East', 'Vertical')
% subplot(212)
% plot(p.t,output.raps_constraint(4,:))
% hold on
% plot(p.t,output.raps_constraint(5,:))
% hold on
% plot(p.t,output.raps_constraint(6,:))
% grid on
% legend('North', 'East', 'Vertical')
%% POS
% Indices based on hor_err
inds1 = output.hor_err <= 1.5;
inds2 = (output.hor_err > 1.5) & (output.hor_err <= 5);
inds3 = (output.hor_err > 5) & (output.hor_err <= 50);
inds4 = output.hor_err > 50;

% Extracting data based on indices
const1 = output.raps_constraint(:, inds1);
const2 = output.raps_constraint(:, inds2);
const3 = output.raps_constraint(:, inds3);
const4 = output.raps_constraint(:, inds4);

% Plotting
figure;

% North
subplot(3, 1, 1);
hold on;
histogram(const1(1, :), 'FaceColor', 'r', 'DisplayName', '<=1.5');
histogram(const2(1, :), 'FaceColor', 'g', 'DisplayName', '1.5-5');
% histogram(const3(1, :), 'FaceColor', 'b', 'DisplayName', '5-50');
% histogram(const4(1, :), 'FaceColor', 'y', 'DisplayName', '>50');
title('North');
legend;
hold off;

% East
subplot(3, 1, 2);
hold on;
histogram(const1(2, :), 'FaceColor', 'r', 'DisplayName', '<=1.5');
histogram(const2(2, :), 'FaceColor', 'g', 'DisplayName', '1.5-5');
% histogram(const3(2, :), 'FaceColor', 'b', 'DisplayName', '5-50');
% histogram(const4(2, :), 'FaceColor', 'y', 'DisplayName', '>50');
title('East');
legend;
hold off;

% Vertical
subplot(3, 1, 3);
hold on;
histogram(const1(3, :), 'FaceColor', 'r', 'DisplayName', '<=1.5');
histogram(const2(3, :), 'FaceColor', 'g', 'DisplayName', '1.5-5');
% histogram(const3(3, :), 'FaceColor', 'b', 'DisplayName', '5-50');
% histogram(const4(3, :), 'FaceColor', 'y', 'DisplayName', '>50');
title('Vertical');
legend;
hold off;

%% Vel

% Plotting
figure;

% North
subplot(3, 1, 1);
hold on;
histogram(const1(4, :), 'FaceColor', 'r', 'DisplayName', '<=1.5');
histogram(const2(4, :), 'FaceColor', 'g', 'DisplayName', '1.5-5');
histogram(const3(4, :), 'FaceColor', 'b', 'DisplayName', '5-50');
histogram(const4(4, :), 'FaceColor', 'y', 'DisplayName', '>50');
title('North');
legend;
hold off;

% East
subplot(3, 1, 2);
hold on;
histogram(const1(5, :), 'FaceColor', 'r', 'DisplayName', '<=1.5');
histogram(const2(5, :), 'FaceColor', 'g', 'DisplayName', '1.5-5');
histogram(const3(5, :), 'FaceColor', 'b', 'DisplayName', '5-50');
histogram(const4(5, :), 'FaceColor', 'y', 'DisplayName', '>50');
title('East');
legend;
hold off;

% Vertical
subplot(3, 1, 3);
hold on;
histogram(const1(6, :), 'FaceColor', 'r', 'DisplayName', '<=1.5');
histogram(const2(6, :), 'FaceColor', 'g', 'DisplayName', '1.5-5');
histogram(const3(6, :), 'FaceColor', 'b', 'DisplayName', '5-50');
histogram(const4(6, :), 'FaceColor', 'y', 'DisplayName', '>50');
title('Vertical');
legend;
hold off;