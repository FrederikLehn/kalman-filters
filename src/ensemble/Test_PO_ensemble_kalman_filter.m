%% Loading data
clear
clc
z = csvread('Satelliteorbit.csv');
z = z';

n_mes = length(z(1,:));

u = [0 0 0]';
Q = [500^2 0 0; 0 0.005^2 0; 0 0 0.005^2];
R = [2000^2 0; 0 0.03^2];

%start guess
x_0 = [z(:,1); 0];
P_0 = Q;
N = 100;

[x_hat, P] = PO_ensemble_kalman_filter(z, u, @f, @h, Q, R, x_0, P_0, N);

%calculating confidence intervals
r_conf_upper = x_hat(1,:)+1.96.*sqrt(P(1,1,:));
r_conf_upper = r_conf_upper(1:n_mes);
r_conf_lower = x_hat(1,:)-1.96.*sqrt(P(1,1,:));
r_conf_lower = r_conf_lower(1:n_mes);

theta_conf_upper = x_hat(2,:)+1.96.*sqrt(P(2,2,:));
theta_conf_upper = theta_conf_upper(1:n_mes);
theta_conf_lower = x_hat(2,:)-1.96.*sqrt(P(2,2,:));
theta_conf_lower = theta_conf_lower(1:n_mes);

v_conf_upper = x_hat(3,:)+1.96.*sqrt(P(3,3,:));
v_conf_upper = v_conf_upper(1:n_mes);
v_conf_lower = x_hat(3,:)-1.96.*sqrt(P(3,3,:));
v_conf_lower = v_conf_lower(1:n_mes);

%% Plotting
pl = plot(1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10);
c = get(pl,'Color');

% Plotting radius
t = 1:(length(z(1,:)));
plot(t, z(1,:))
hold on
plot(t, x_hat(1,:), 'LineWidth',1.2)
plot(t, r_conf_upper, 'k--', 'LineWidth',1.2)
plot(t, r_conf_lower, 'k--', 'LineWidth',1.2)
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('Radius, $r_t$ [km]', 'interpreter', 'latex')
leg = legend('Measurements, $r_t^m$', 'Estimate, $\hat{r}_t^p$','95\% Conf. int.', 'Location', 'NW');
set(leg,'interpreter','latex')
grid on

% Plotting angle
figure
plot(t, z(2,:))
hold on
plot(t, x_hat(2,:), 'LineWidth',1)
plot(t, theta_conf_upper, 'k--', 'LineWidth',1)
plot(t, theta_conf_lower, 'k--', 'LineWidth',1)
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('Angle, $\theta_t$ [radian]', 'interpreter', 'latex')
leg = legend('Measurements, $\theta_t^m$', 'Estimate, $\hat{\theta}_t^p$','95\% Conf. int.', 'Location', 'NW');
set(leg,'interpreter','latex')
grid on

% Plotting angular velocity
figure
hold on
plot(t, x_hat(3,:), 'LineWidth',1, 'Color', c{2})
plot(t, v_conf_upper, 'k--', 'LineWidth',1)
plot(t, v_conf_lower, 'k--', 'LineWidth',1)
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('Angular Velocity, $v_{\theta,t}$ [$^o/hour$]', 'interpreter', 'latex')
leg = legend('Estimate, $\hat{v}_{\theta,t}^p$','95\% Conf. int.', 'Location', 'NW');
set(leg,'interpreter','latex')
grid on

% plotting in cartesian coordinates
mes_x_car = z(1,:).*cos(z(2,:));
mes_y_car = z(1,:).*sin(z(2,:));

est_x_car = x_hat(1,:).*cos(x_hat(2,:));
est_y_car = x_hat(1,:).*sin(x_hat(2,:));
x_car_conf_upper = r_conf_upper.*cos(theta_conf_upper);
y_car_conf_upper = r_conf_upper.*sin(theta_conf_upper);
x_car_conf_lower = r_conf_lower.*cos(theta_conf_lower);
y_car_conf_lower = r_conf_lower.*sin(theta_conf_lower);

figure
plot(mes_x_car, mes_y_car)
hold on
plot(est_x_car, est_y_car, 'LineWidth', 1.1)
plot(x_car_conf_upper,y_car_conf_upper, 'k--', 'LineWidth', 1.1)
plot(x_car_conf_lower,y_car_conf_lower, 'k--', 'LineWidth', 1.1)
xlabel('$1^{st}$ lateral distance from earth, [km]', 'interpreter', 'latex')
ylabel('$2^{nd}$ lateral distance from earth, [km]', 'interpreter', 'latex')
leg = legend('Measurement, $(x,y)$', 'Estimate, $(\hat{x}, \hat{y})$','95\% Conf. int.', 'Location', 'SW');
set(leg,'interpreter','latex')
ylim([-10000, 45000])
grid on

%Root mean sqauerd error as function of t (RMS)
RMS_r = zeros(1, n_mes); %RMS of radius
RMS_theta = zeros(1, n_mes); %RMS of theta (angle)
for i = 1:n_mes
    RMS_r(i) = 1/i*sum((z(1,1:i)-x_hat(1,1:i)).^2);
    RMS_theta(i) = 1/i*sum((z(2,1:i)-x_hat(2,1:i)).^2);
end
RMS_r = sqrt(RMS_r);
RMS_theta = sqrt(RMS_theta);

figure
plot(t(2:end), RMS_r(2:end))
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('RMS', 'interpreter', 'latex')
leg = legend('RMS$_r$', 'Location', 'NE');
set(leg,'interpreter','latex')
grid on

figure
plot(t(2:end), RMS_theta(2:end))
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('RMS', 'interpreter', 'latex')
leg = legend('RMS$_{\theta}$', 'Location', 'NE');
set(leg,'interpreter','latex')
grid on


%% Plots of state variables at multiple ensemble sizes
[x_hat_10, ~] = PO_ensemble_kalman_filter(z, u, @f, @h, Q, R, x_0, P_0, 10);
[x_hat_100, ~] = PO_ensemble_kalman_filter(z, u, @f, @h, Q, R, x_0, P_0, 100);
[x_hat_1000, ~] = PO_ensemble_kalman_filter(z, u, @f, @h, Q, R, x_0, P_0, 1000);

% Plotting radius
t = 1:(length(z(1,:)));
% plot(t, z(1,:), 'LineWidth',1.1)
hold on
plot(t, x_hat(1,:), 'LineWidth',1.2, 'Color', c{2})
plot(t, x_hat_10(1,:), 'LineWidth',1.1, 'Color', c{3})
plot(t, x_hat_100(1,:), 'LineWidth',1.1, 'Color', c{4})
plot(t, x_hat_1000(1,:), 'LineWidth',1.1, 'Color', c{5})
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('Radius, $r_t$ [km]', 'interpreter', 'latex')
leg = legend('KF, $\hat{r}_t^p$', 'EnKF, $\hat{r}_t^p (N=10)$', 'EnKF, $\hat{r}_t^p (N=100)$', 'EnKF, $\hat{r}_t^p (N=1000)$', 'Location', 'NW');
set(leg,'interpreter','latex')
grid on

% Plotting angle
figure
% plot(t, z(2,:), 'LineWidth',1.1)
hold on
plot(t, x_hat(2,:), 'LineWidth',1.2, 'Color', c{2})
plot(t, x_hat_10(2,:), 'LineWidth',1.1, 'Color', c{3})
plot(t, x_hat_100(2,:), 'LineWidth',1.1, 'Color', c{4})
plot(t, x_hat_1000(2,:), 'LineWidth',1.1, 'Color', c{5})
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('Angle, $\theta_t$ [radian]', 'interpreter', 'latex')
leg = legend('Kf, $\hat{\theta}_t^p$', 'EnKF, $\hat{\theta}_t^p (N=10)$', 'EnKF, $\hat{\theta}_t^p (N=100)$', 'EnKF, $\hat{\theta}_t^p (N=1000)$', 'Location', 'NW');
set(leg,'interpreter','latex')
grid on

% Plotting angular velocity
figure
hold on
plot(t, x_hat(3,:), 'LineWidth',1.2, 'Color', c{2})
plot(t, x_hat_10(3,:), 'Color', c{3}, 'LineWidth',1.1)
plot(t, x_hat_100(3,:), 'Color', c{4}, 'LineWidth',1.1)
plot(t, x_hat_1000(3,:), 'Color', c{5}, 'LineWidth',1.1)
xlabel('Time-step, $t$ [-]', 'interpreter', 'latex')
ylabel('Angular Velocity, $v_{\theta,t}$ [$^o/hour$]', 'interpreter', 'latex')
leg = legend('KF, $\hat{v}_{\theta,t}^p$','EnKF, $\hat{v}_{\theta,t}^p (N=10)$','EnKF, $\hat{v}_{\theta,t}^p (N=100)$','EnKF, $\hat{v}_{\theta,t}^p (N=1000)$', 'Location', 'SE');
set(leg,'interpreter','latex')
grid on

% plotting in cartesian coordinates
mes_x_car = z(1,:).*cos(z(2,:));
mes_y_car = z(1,:).*sin(z(2,:));

est_x_car_10 = x_hat_10(1,:).*cos(x_hat_10(2,:));
est_y_car_10 = x_hat_10(1,:).*sin(x_hat_10(2,:));
est_x_car_100 = x_hat_100(1,:).*cos(x_hat_100(2,:));
est_y_car_100 = x_hat_100(1,:).*sin(x_hat_100(2,:));
est_x_car_1000 = x_hat_1000(1,:).*cos(x_hat_1000(2,:));
est_y_car_1000 = x_hat_1000(1,:).*sin(x_hat_1000(2,:));
figure
% plot(mes_x_car, mes_y_car, 'LineWidth',1.1)
hold on
plot(est_x_car, est_y_car, 'LineWidth', 1.2, 'Color', c{2})
plot(est_x_car_10, est_y_car_10, 'LineWidth',1.1, 'Color', c{3})
plot(est_x_car_100, est_y_car_100, 'LineWidth',1.1, 'Color', c{4})
plot(est_x_car_1000, est_y_car_1000, 'LineWidth',1.1, 'Color', c{5})

xlabel('$1^{st}$ lateral distance from earth, [km]', 'interpreter', 'latex')
ylabel('$2^{nd}$ lateral distance from earth, [km]', 'interpreter', 'latex')
leg = legend('KF, $(\hat{x}, \hat{y})$', 'EnKF, $(\hat{x}, \hat{y}) (N=10)$', 'EnKF, $(\hat{x}, \hat{y}) (N=100)$', 'EnKF, $(\hat{x}, \hat{y}) (N=1000)$', 'Location', 'SW');
set(leg,'interpreter','latex')
ylim([-10000, 45000])
grid on

