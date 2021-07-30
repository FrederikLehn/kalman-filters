%% Grid
clear
clc
rng(12345)
[nx, ny, nz] = deal( 11,  11, 3);
[Dx, Dy, Dz] = deal(200, 200, 60);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

N = 50;
rockEnsemble = cell(N, 1);
referencePerm = 0;
for i = 1:N
    K = convertFrom(gaussianField(G.cartDims, [0.1, 5]), milli*darcy);
    rock.perm = reshape(K(G.cells.indexMap), [], 1);
    rock.poro  = .3*ones(G.cells.num, 1);
    rock.ntg   = ones(G.cells.num, 1);
    rockEnsemble{i} = rock;
    
    referencePerm = referencePerm + rock.perm;
end
referencePerm = referencePerm./N;
% referencePerm = rockEnsemble{1}.perm;
figure
title('Permeability Field (Reference Model)')
cla, plotCellData(G, referencePerm)%, view(3);
h = colorbar;
set(get(h,'label'),'string','Permeability, $\textbf{K}$ [mD]','Rotation',90.0, 'interpreter', 'latex');
caxis([2*10^(-15), 6*10^(-15)])
xlabel('x-direction, $D_x$ [m]', 'interpreter', 'latex')
ylabel('y-direction, $D_y$ [m]', 'interpreter', 'latex')
zlabel('z-direction, $D_z$ [m]', 'interpreter', 'latex')
%% Assume a single horizontal well
layer = fix(nz/2);

W = addWell([], G, rock, (layer*nx*ny+round(nx/2)) : ny : (layer+1)*nx*ny, ...
            'InnerProduct', 'ip_tpf', ...
            'Type', 'bhp' , 'Val', 1.0e5, ...
            'Radius', 0.1, 'Dir', 'y');

% To check if the wells are placed as we wanted them, we plot them
clf
plotGrid(G, 'FaceColor', 'none'); view(3);
[ht, htxt, hs] = plotWell(G, W, 'radius', 0.1, 'height', 2);
set(htxt, 'FontSize', 16);
xlabel('x-direction, $D_x$ [m]', 'interpreter', 'latex')
ylabel('y-direction, $D_y$ [m]', 'interpreter', 'latex')
zlabel('z-direction, $D_z$ [m]', 'interpreter', 'latex')

%% create aritifical data
%perform a generic reservoir simulation
n_mes = 30;
T      = 365*day();
dt     = T/n_mes;
x = [rockEnsemble{1}.perm; ones(2*G.cells.num, 1); 1];

[~, sol, ~] = reservoir_simulator(x(:,1), G, W, rockEnsemble{1}, T, dt);

x = repmat(x, 1, N);
for i = 2:N
    x(1:G.cells.num,i) = rockEnsemble{i}.perm;
end

times = convertTo([sol(2:end).time], day)*day();
referenceProduction = convertTo(-[sol(2:end).qS], meter^3/day)';

R = 10;
offset = 10;
measuredProduction = (referenceProduction + offset + mvnrnd(0, R, n_mes))';

%% Perform ensemble kalman filtering
x_0 = x;
Q = zeros(length(x)); %assuming no model error

[x_hat, x_hat_i, RMSE, ens_prod] = HM_ensemble_kalman_filter(measuredProduction, times, dt, @reservoir_simulator, G, W, rockEnsemble, @limit_fun, Q, R, x_0, N);

%% Plot of mean production curves
clf
plot(times./day(), measuredProduction, 'o', 'MarkerSize', 3)
hold on
plot(times./day(), referenceProduction, 'LineWidth', 1.5)
plot(times./day(), x_hat(end,:), 'LineWidth', 1.5)
xlabel('Time [days]', 'interpreter', 'latex')
ylabel('Oil Flowrate $\left[\frac{m^3}{day}\right]$', 'interpreter', 'latex')
leg = legend('Measured Production', 'Initial Forecast', 'History Matched Forecast', 'Location', 'NE');
set(leg,'interpreter','latex')
grid on

%% Plot of RMSE
RMSE_ref = sqrt(1/n_mes*sum((measuredProduction-referenceProduction').^2));
plot(0:n_mes,[RMSE_ref, RMSE]./RMSE_ref, 'LineWidth', 1.5)
hold on
minRMSE = min(RMSE./RMSE_ref);
plot([0, n_mes], [minRMSE, minRMSE], 'k--')
xlabel('Assimilation number [-]', 'interpreter', 'latex')
ylabel('RMSE', 'interpreter', 'latex')
leg = legend('RMSE', 'Asymptotic Minimum', 'Location', 'NE');
set(leg,'interpreter','latex')
grid on
ylim([0, 1])

%% plot of permeability fields
figure
subplot(2,1,1)
title('Permeability Field (Reference Model)')
cla, plotCellData(G, referencePerm), view(3);
h = colorbar;
set(get(h,'label'),'string','Permeability, $\textbf{K}$ [mD]','Rotation',90.0, 'interpreter', 'latex');
caxis([2*10^(-15), 6*10^(-15)])
subplot(2,1,2)
title('Permeability Field (Mean History Matched Model)')
cla, plotCellData(G, x_hat(1:G.cells.num, end)), view(3);   
h = colorbar;
set(get(h,'label'),'string','Permeability, $\textbf{K}$ [mD]','Rotation',90.0, 'interpreter', 'latex');
caxis([2*10^(-15), 6*10^(-15)])
%% Plot of change in parameters
n_theta = length(referencePerm);
param_norm = zeros(1, n_mes-1);
for i = 2:n_mes
    param_norm(i-1) = norm(x_hat(1:n_theta,i)-x_hat(1:n_theta,i-1))./norm(x_hat(1:n_theta,i-1));
end
figure
plot(2:n_mes, param_norm,'LineWidth', 1.5)
xlabel('Assimilation number [-]', 'interpreter', 'latex')
ylabel('Relative change $\frac{||\theta_i-\theta_{i-1}||}{||\theta_{i-1}||}$', 'interpreter', 'latex')
grid on

%% Plot of ensemble production curves
figure
plot(times./day(), measuredProduction, 'o', 'MarkerSize', 3)
hold on
for i = 1:N
    [~, ~, q_o] = reservoir_simulator(x(:,i), G, W, rockEnsemble{i}, T, dt);
    plot(times./day(), q_o, 'LineWidth', 0.5)
end
xlabel('Time [days]', 'interpreter', 'latex')
ylabel('Oil Flowrate $\left[\frac{m^3}{day}\right]$', 'interpreter', 'latex')
leg = legend('Measured Production', ['N=', num2str(N), ' Initial Forecasts'], 'Location', 'NE');
set(leg,'interpreter','latex')
grid on


figure
plot(times./day(), measuredProduction, 'o', 'MarkerSize', 3)
hold on
for i = 1:N
    plot(times./day(), ens_prod(i,:), 'LineWidth', 0.5)
end
xlabel('Time [days]', 'interpreter', 'latex')
ylabel('Oil Flowrate $\left[\frac{m^3}{day}\right]$', 'interpreter', 'latex')
leg = legend('Measured Production', ['N=', num2str(N), ' History Matched Forecasts'], 'Location', 'NE');
set(leg,'interpreter','latex')
grid on
