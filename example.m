% event example
P = [0 0 1 3 6 3 1 0 0];     % mm/h
Q = [0 0 0.1 0.4 0.9 0.7 0.3 0.1 0]; % mm/h
uh = struct('type','refh','Up',0.65,'Uk',0.8);   % standrad ReFH shape
res = estimate_Tp_timevariance(P, Q, uh);
fprintf('Tp_est = %.2f h (c=%.3f)\n', res.Tp_est, res.c);
