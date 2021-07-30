function [x_hat, P, S] = extended_kalman_filter(z, u, f, J_f, h, J_h, Q, R, x_0, P_0)
    %% Extended Kalman Filter Implementation
    % The Extended Kalman filter as applied to the stationary stochastic state-space model:
    %   x_t = f(x_{t-1},u_t) + w_t
    %   z_t = h(x_t) + v_t
    %
    % Inputs:
    %   - z: Observations (n_sta x n_mes vector)
    %   - u: Control-vector (n_con x 1 vector)
    %   - f: State-transition model (function handle of x and u)
    %   - J_f: Jacobian of f (function handle of x and u)
    %   - h: Observation model (function handle of u)
    %   - J_h: Jacobian of h (function handle of x)
    %   - Q: Covariance of the system (n_sta x n_sta matrix)
    %   - R: Covariance of the observations (n_obs x n_obs matrix)
    %   - x_0: Start guess for the state variable (n_obs x 1 vector)
    %   - P_0: Start guess for variance of the state (n_sta x n_sta matrix)
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
    x_hat       = zeros(n_sta_var,n_mes);
    y_tilde     = zeros(n_obs_var, n_mes);
    P           = zeros(n_sta_var,n_sta_var,n_mes);
    S           = zeros(n_obs_var,n_obs_var,n_mes);

%% Initiazation
    %Assigning start guesses
    x_hat(:,1)      = x_0;
    P(:,:,1)        = P_0;

    %Updating Jacobians
    H = J_h(x_hat(:,1));
    Ht = H';
    
    %Update
    y_tilde(:,1)    = z(:,1) - h(x_hat(:,1));
    S(:,:,1)        = H*P(:,:,1)*Ht + R;
    K_t             = P(:,:,1)*Ht*(S(:,:,1)\I);
    x_hat(:,1)      = x_hat(:,1) + K_t*y_tilde(:,1);
    P(:,:,1)        = P(:,:,1) - K_t*H*P(:,:,1);
    y_tilde(:,1)    = z(:,1) - H*x_hat(:,1);
    
%% Kalman Iterations
    for t = 2:n_mes
        %Prediction
        F = J_f(x_hat(:,t-1), u);
        x_hat(:,t)  = f(x_hat(:,t-1), u);
        P(:,:,t)    = F*P(:,:,t-1)*F' + Q;
       
        %Update
        H = J_h(x_hat(:,t));
        Ht = H';
        
        y_tilde(:,t)    = z(:,t) - H*x_hat(:,t);
        S(:,:,t)        = H*P(:,:,t)*Ht + R;
        K_t             = P(:,:,t)*Ht*(S(:,:,t)\I);
        x_hat(:,t)      = x_hat(:,t) + K_t*y_tilde(:,t);
        P(:,:,t)        = P(:,:,t) - K_t*H*P(:,:,t); 
        y_tilde(:,t)    = z(:,t) - H*x_hat(:,t);
    end
end
