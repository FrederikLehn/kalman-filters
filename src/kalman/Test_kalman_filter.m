%% Loading data
clear
clc
z = csvread('Satelliteorbit.csv');
z = z';

n_mes = length(z(1,:));

F = [1 0 0; 0 1 1; 0 0 1];
B = [0 0 0; 0 0 0; 0 0 0];
u = [0 0 0]';
Q = [500^2 0 0; 0 0.005^2 0; 0 0 0.005^2];

H = [1 0 0; 0 1 0];
R = [2000^2 0; 0 0.03^2];

%start guess
x_0 = [z(:,1); 0];
P_0 = Q;

[x_hat, P, S] = kalman_filter(z, u, F, B, H, Q, R, x_0, P_0);

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
h = plot(1:10,1:10,1:10,1:10,1:10,1:10);
c = get(h,'Color');

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
