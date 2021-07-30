function out = f(x_hat, u)
    F = [1 0 0; 0 1 1; 0 0 1];
    B = [0 0 0; 0 0 0; 0 0 0];
    out = F*x_hat + B*u;
end