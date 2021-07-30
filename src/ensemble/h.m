function out = h(x_hat)
    H = [1 0 0; 0 1 0];
    out = H*x_hat;
end