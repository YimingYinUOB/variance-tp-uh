% ReFH (standard)
uh = struct('type','refh','Up',0.65,'Uk',0.8);
out = estimate_Tp_timevariance(P, Q, uh);

% ReFH triangular
uh = struct('type','triangular');
out = estimate_Tp_timevariance(P, Q, uh);

% Gamma (Nash) with alpha=5
uh = struct('type','gamma','alpha',5);
out = estimate_Tp_timevariance(P, Q, uh);
