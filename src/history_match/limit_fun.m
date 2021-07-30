function out = limitFun(x_hat)
H = zeros(1, length(x_hat));
H(end) = 1;
out = H*x_hat;
end

