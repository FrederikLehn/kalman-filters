function [x_hat, x_hat_i, RMSE, ens_prod] = HM_ensemble_kalman_filter(z, times, dt, res_sim, G, W, rockEnsemble, h, Q, R, x_0, N)
    %% Pertubed Observation Ensemble Kalman Filter Implementation
    % The Pertubed Observation Ensemble Kalman filter as applied to the stationary stochastic state-space model:
    %   x_t = f(x_{t-1},u_t) + w_t
    %   z_t = h(x_t) + v_t
    %
    % Implemented observation matrix-free to allow non-linear h(x_t)
    %
    % Inputs:
    %   - z: Observations (n_sta x n_mes vector)
    %   - res_sim: reservoir simulation function (n_con x 1 vector)
    %   - G: MRST geometry structures
    %   - W: MRST well structure
    %   - rockEnsemble: cell-array of MRST rock structures
    %   - h: Observation model (function handle of u)
    %   - Q: Covariance of the system (n_sta x n_sta matrix)
    %   - R: Covariance of the observations (n_obs x n_obs matrix)
    %   - x_0: Start guess for the state variable (n_obs x 1 vector)
    %   - P_0: Start guess for variance of the state (n_sta x n_sta matrix)
    %   - N: The number of realisations in the ensemble
    %
    %   n_sta is the number of state variables, n_mes is the number of
    %   measurements, n_con is the number of control variables, n_obs is
    %   the number of observed variables
    %
    % Outputs:
    %   - x_hat: Estimated state variable at each time-step (n_sta x n_mes vector)
    %   - P: Variance of the state at each time-step (n_sta x n_sta x n_mes 3D matrix)
    %   - S: Variance of the observations at each time-step (n_obs x n_obs x n_mes 3D matrix)
    %
%% Preallocating memory
    %Preallocating
    n_mes       = length(z(1,:)); %number of measurements
    n_sta_var   = length(x_0); %number of state variables
    n_obs_var   = length(z(:,1)); %number observed variables 
    I           = eye(n_obs_var);
    mu0_obs     = zeros(n_obs_var, 1);
    mu0_sta     = zeros(n_sta_var, 1);
    x_hat       = zeros(n_sta_var,n_mes);
    y_hat       = zeros(n_obs_var, n_mes);
    y_hat_i     = zeros(n_obs_var, N);
    ens_prod    = zeros(N, n_mes);
    RMSE        = zeros(1, n_mes);
    
%% Initiazation
    %Assigning ensemble
    x_hat_i = x_0;
   
%% Kalman Iterations
    for t = 1:n_mes
        w = mvnrnd(mu0_sta, Q, N)';
        v = mvnrnd(mu0_obs, R, N)';
        R_ek = (v*v')./(N-1);
        
        %predicting ensemble
        for i = 1:N
            [xsol, ~, ~] = res_sim(x_hat_i(:,i), G, W, rockEnsemble{i}, times(t), dt);
            x_hat_i(:,i) = xsol + w(:,i);
        end
        x_hat(:,t) = mean(x_hat_i, 2);
        
        %Updating ensemble
        for i = 1:N
            y_hat_i(:,i) = h(x_hat_i(:,i));
        end
        y_hat(:,t) = mean(y_hat_i, 2);
        P_yt = 1/(N-1).*(y_hat_i-y_hat(:,t))*(y_hat_i-y_hat(:,t))'+R_ek;
        P_xyt = 1/(N-1).*(x_hat_i-x_hat(:,t))*(y_hat_i-y_hat(:,t))';
        K_t = P_xyt*(P_yt\I);
        
        for i = 1:N
            x_hat_i(:,i) = x_hat_i(:,i) + K_t*(z(:,t) + v(:,i) - y_hat_i(:,i));
        end
        x_hat(:,t) = mean(x_hat_i, 2);
        
        %Finding the RMSE
        [~, ~, pred_t] = res_sim(x_hat(:,t), G, W, rockEnsemble{i}, times(end), dt);
        RMSE(t) = sqrt(1/n_mes*sum((z-pred_t).^2));
        
        disp(['Assimilation number: ', num2str(t), ' completed'])
    end
    
%% Additional output
    for i = 1:N
        [~, ~, q_o] = res_sim(x_hat_i(:,i), G, W, rockEnsemble{i}, times(end), dt);
        ens_prod(i,:) = q_o;
    end
end
