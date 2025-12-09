function out = estimate_Tp_timevariance(rainfall, runoff, uh, t)
% estimate_Tp_timevariance
% Estimate Tp from rainfall & runoff using the variance-based method.
%
% USAGE:
%   out = estimate_Tp_timevariance(P, Q, uh)
%   out = estimate_Tp_timevariance(P, Q, uh, t)
%
% INPUTS
%   rainfall : 1xN or Nx1 numeric, rainfall intensity (e.g., mm/h), non-negative
%   runoff   : 1xN or Nx1 numeric, direct runoff (mm/h), non-negative
%   uh       : struct specifying UH shape (only ONE free timing parameter Tp):
%              % -- allowed values only:
%              %   uh.type = 'refh'       (ReFH standard; Up, Uk fixed)
%              %   uh.type = 'triangular' (ReFH triangular; Up=0.5, Uk=1)
%              %   uh.type = 'gamma'      (Nash–Gamma; alpha>1)
%              %
%              % fields per type:
%              %   'refh'      -> uh.Up (default 0.65), uh.Uk (default 0.8)
%              %   'triangular'-> no extra fields (uses Up=0.5, Uk=1)
%              %   'gamma'     -> uh.alpha (required, >1)
%   t        : (optional) 1xN or Nx1 time in hours or datetime; if omitted, uses 0:(N-1)
%
% OUTPUT (struct)
%   out.Tp_est : estimated Tp (hours)
%   out.Var_tP : time variance of rainfall (h^2)
%   out.Var_tQ : time variance of runoff   (h^2)
%   out.Var_TU : inferred UH variance = Var(tQ) - Var(tP) (h^2)
%   out.c      : shape constant c = sqrt(Var(T_U^(1)))
%
% NOTES
% - Scaling rainfall by a constant (e.g., an event runoff coefficient) does not change Var(t_P).
% - Ensure 'runoff' is direct runoff (baseflow removed at event scale if needed).

    % ---------- inputs & time vector ----------
    P = rainfall(:)'; Q = runoff(:)';
    if numel(P) ~= numel(Q)
        error('rainfall and runoff must have the same length.');
    end
    N = numel(P);

    if nargin < 4 || isempty(t)
        t = 0:(N-1);                 % hours, starting at 0
    else
        t = t(:)';                   % row
        if numel(t) ~= N, error('t must have the same length as rainfall/runoff.'); end
        if isdatetime(t)
            t = hours(t - t(1));     % convert and start from 0
        else
            t = t - t(1);            % start from 0
        end
    end

    % ---------- time variance (volume-weighted) ----------
    Var_tP = local_time_variance(t, P);
    Var_tQ = local_time_variance(t, Q);
    Var_TU = max(Var_tQ - Var_tP, 0);  % numerical safeguard

    % ---------- shape constant c (only 'refh' | 'triangular' | 'gamma') ----------
    c = local_shape_constant_strict(uh);

    % ---------- Tp ----------
    Tp_est = sqrt(Var_TU) / c;

    % ---------- pack output ----------
    out = struct('Tp_est',Tp_est,'Var_tP',Var_tP,'Var_tQ',Var_tQ,'Var_TU',Var_TU,'c',c);
end

% ===== helper: time variance (volume-weighted) =====
function v = local_time_variance(t, x)
    t = t(:)'; x = max(0, x(:))';
    if numel(t) < 2 || all(x==0)
        v = NaN; return;
    end
    % robust dt (handles non-even spacing)
    dt_last = max(1e-12, t(end) - t(max(1,end-1)));
    dt = [diff(t), dt_last];
    w  = x .* dt;                     % volume weights
    W  = sum(w);
    if W <= 0, v = NaN; return; end
    Et  = sum(t .* w) / W;
    Et2 = sum((t.^2) .* w) / W;
    v   = Et2 - Et^2;                 % h^2
end

% ===== helper: UH shape constant c = sqrt(Var(T_U^(1))) with strict types =====
function c = local_shape_constant_strict(uh)
    if ~isfield(uh,'type')
        error('uh.type is required and must be one of: ''refh'', ''triangular'', ''gamma''.');
    end

    tname = lower(string(uh.type));
    valid = ["refh","triangular","gamma"];
    if ~any(tname == valid)
        error('Invalid uh.type = ''%s''. Allowed: ''refh'', ''triangular'', ''gamma''.', uh.type);
    end

    switch char(tname)
        case 'refh'
            if ~isfield(uh,'Up'), uh.Up = 0.65; end
            if ~isfield(uh,'Uk'), uh.Uk = 0.8;  end
            % Hourly discrete grid for ReFH (Up=0.65, Uk=0.8): Var ≈ 0.533 -> c ≈ 0.73
            c = sqrt(0.533);

        case 'triangular'
            % ReFH triangular (Up=0.5, Uk=1): Var = 13/18
            c = sqrt(13/18);   % ≈ 0.84984

        case 'gamma'
            if ~isfield(uh,'alpha')
                error('For uh.type=''gamma'', uh.alpha (>1) is required.');
            end
            alpha = uh.alpha;
            if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 1
                error('For uh.type=''gamma'', uh.alpha must be a finite scalar > 1.');
            end
            % Var(T_U^(1)) = alpha/(alpha-1)^2  => c = sqrt(alpha)/(alpha-1)
            c = sqrt(alpha) / (alpha - 1);
    end
end

