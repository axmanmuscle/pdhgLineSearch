function xstar = testProx(x0, f, proxf)
%testProx tests proximal operators
%   prox operators are defined as
% prox_{tf}(u) = argmin_{x} f(x) + 1/(2t) || u - x ||_2^2

rng(2024);
t = 2;
prox_tfx = proxf(x0, t);

options = optimoptions('fminunc', 'MaxFunctionEvaluations', 100, ...
                       'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-12, ...
                       'Display', 'iter-detailed');

prox_f_def = @(x, u) f(x) + 1/(2*t) * norm(x - u, 'fro')^2;

iters = 1;
for i = 1:iters
    test_x = randn(size(x0));
    test_fun = @(x) prox_f_def(x, test_x);

    xStar = fminunc(test_fun, randn(size(x0)), options);
end

end