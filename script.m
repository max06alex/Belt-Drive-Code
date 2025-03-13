clc, clear all, close all;
Data_90 = csvread("90_beta_belt.csv", 1, 0);
Data_180 = csvread('180_beta_belt.csv', 1, 0);
Data_270 = csvread('270_beta_belt.csv', 1, 0);
Data_360 = csvread('360_beta_belt.csv', 1, 0);

T2_90 = Data_90(:, 1);
T1_90 = Data_90(:, 2);
omega_90 = Data_90(:, 4);
current_90 = Data_90(:, 5);

T2_180 = Data_180(:, 1);
T1_180 = Data_180(:, 2);
omega_180 = Data_180(:, 4);
current_180 = Data_180(:, 5);

T2_270 = Data_270(:, 1);
T1_270 = Data_270(:, 2);
omega_270 = Data_270(:, 4);
current_270 = Data_270(:, 5);

T2_360 = Data_360(:, 1);
T1_360 = Data_360(:, 2);
omega_360 = Data_360(:, 4);
current_360 = Data_360(:, 5);
Voltage = 12;

p_90 = polyfit(T1_90, T2_90,1);
bestfit_90 = p_90(1)*T1_90 + p_90(2);

p_180 = polyfit(T1_180, T2_180, 1);
bestfit_180 = p_180(1)*T1_180 + p_180(2);

p_270 = polyfit(T1_270, T2_270, 1);
bestfit_270 = p_270(1)*T1_270 + p_270(2);

p_360 = polyfit(T1_360, T2_360, 1);
bestfit_360 = p_360(1)*T1_360 +p_360(2);
actual_points = [p_90(1) p_180(1) p_270(1) p_360(1)];
figure
plot(T1_90, T2_90, 'bx');
hold on
plot(T1_90, bestfit_90, 'b--', 'LineWidth', 1.5);
plot(T1_180, T2_180, 'gx');
plot(T1_180, bestfit_180, 'g--', 'LineWidth', 1.5);
plot(T1_270, T2_270, 'mx');
plot(T1_270, bestfit_270, 'm--', 'LineWidth',1.5);
plot(T1_360, T2_360, 'kx');
plot(T1_360, bestfit_360, 'k--', 'LineWidth',1.5);
xlabel('T1 (N)', 'FontSize', 16);
ylabel('T2 (N)', 'FontSize', 16);
legend('90 Deg.', '90 Deg. Fit','180 Deg.', '180 Deg. Fit', '270 Deg.', '270 Deg. Fit', '360 Deg.', '360 Deg. Fit',  'Location', 'southeast');
grid on
hold off

beta = [1.570796 3.141592 4.712389 6.283185];
ln_data = [log(p_90(1)) log(p_180(1)) log(p_270(1)) log(p_360(1))];
p_ln = polyfit(beta, ln_data, 1);
fit_ln = p_ln(1)*beta + p_ln(2);
tension = [p_90(1) p_180(1) p_270(1) p_360(1)];
figure(2)
plot(beta, ln_data, 'bx');
hold on
plot(beta, fit_ln, 'b--', 'LineWidth', 1.5);
xlabel('Contact Angle (Radians)', 'FontSize', 16);
ylabel('ln(T2/T1)', 'FontSize', 16)
xlim([0 7]);
ylim([0 2]);
legend('Data Points', 'Best Fit Line');
xticks(0:pi/2:2*pi);
grid on
hold off

num = [0:0.1:7];
ideal_ratio = exp(num*p_ln(1));
ideal_points = exp(beta*p_ln(1));
figure(3)
plot(beta,ideal_points, 'mx', 'MarkerSize', 6);
hold on
plot(num, ideal_ratio, 'k-', 'LineWidth', 1);
plot(beta, tension ,'bx', 'MarkerSize',6);
grid on
xlim([0 7])
ylim([0 8])
xticks(0:pi/2:2*pi)
xlabel('Contact Angle (Radians)', 'FontSize', 16);
ylabel('Belt Tension Ratio', 'FontSize', 16);
legend('Theoretical Points', 'Theortical Line', 'Experimental Data');
hold off

difference = ideal_points - actual_points;
percentage_error = (abs(difference)./ abs(ideal_points))*100;
percentage_error_ave = mean(percentage_error)
percentage_error_max = max(percentage_error)

input_power = current_180(2:end)*Voltage;
torque = (T2_180(2:end) - T1_180(2:end))*0.05;
output_power = torque.*omega_180(2:end);
efficiency = output_power./input_power;
p_eff = polyfit(torque, efficiency, 2);
num2 = [0:0.001:1];
fit_eff = polyval(p_eff, num2);
figure(4)
plot(torque, efficiency, 'bx');
hold on
plot(num2, fit_eff, 'b-', 'LineWidth', 1.5)
xlabel('Torque (Nm)', 'FontSize', 16);
ylabel('Efficiency', 'FontSize',16);
xlim([0 0.27])
ylim([0 1])
legend('Data Points', 'Best Fit Line - Second Order Model')
grid on
hold off

index = find(efficiency == max(efficiency))
torque_max_eff = torque(index)
