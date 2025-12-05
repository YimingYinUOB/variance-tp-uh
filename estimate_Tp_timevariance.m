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
%   runoff   : 1xN or Nx1 numeric, direct runoff (e.g., mm/h), non-negative
%   uh       : struct specifying UH shape (only ONE free timing parameter Tp):
%              - ReFH (standard): uh.type = 'refh',   uh.Up = 0.65, uh.Uk = 0.8
%              - Triangular     : uh.type = 'triangular'            (Up=0.5, Uk=1)
%              - Gamma (Nash)   : uh.type = 'gamma', uh.alpha > 1
%   t        : (optional) 1xN or Nx1 time in hours; if omitted, uses 0:(N-1)
%
% OUTPUT (struct)
%   out.Tp_est   : estimated Tp (hours)
%   out.Var_tP   : time variance of rainfall (h^2)
%   out.Var_tQ   : time variance of runoff   (h^2)
%   out.Var_TU   : inferred UH variance = Var(tQ) - Var(tP) (h^2)
%   out.c        : shape constant c = sqrt(Var(T_U^(1)))
%
% NOTES
% - Multiplying rainfall by a constant (e.g., using an event runoff coefficient)
%   does NOT change Var(t_P) because the weight normalization cancels the factor.
% - Ensure runoff is direct runoff。

    % ---------- inputs & time vector ----------
    P = rainfall(:)'; Q = runoff(:)';
    if numel(P) ~= numel(Q)
        error('rainfall and runoff must have the same length.');
    end
    N = numel(P);

    if nargin < 4 || isempty(t)
        t = 0:(N-1);             % hours, starting at 0
    else
        t = t(:)';               % to row
        if numel(t) ~= N, error('t must have the same length as rainfall/runoff.'); end
        % make it start from 0 hours for numerical stability
        if isdatetime(t)
            t = hours(t - t(1));
        else
            t = t - t(1);
        end
    end

    % ---------- time variance (volume-weighted) ----------
    Var_tP = local_time_variance(t, P);
    Var_tQ = local_time_variance(t, Q);
    Var_TU = max(Var_tQ - Var_tP, 0);  % numerical safeguard

    % ---------- shape constant c ----------
    c = local_shape_constant(uh);

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
    % dt for possibly non-even spacing
    dt = [diff(t), max(1e-12, t(end)-t(end-1))];
    w  = x .* dt;                     % volume weights
    W  = sum(w);
    if W <= 0, v = NaN; return; end
    Et  = sum(t .* w) / W;
    Et2 = sum((t.^2) .* w) / W;
    v   = Et2 - Et^2;                 % h^2
end

% ===== helper: UH shape constant c = sqrt(Var(T_U^(1))) =====
function c = local_shape_constant(uh)
    if ~isfield(uh,'type'); error('uh.type is required: refh | triangular | gamma'); end
    switch lower(uh.type)
        case 'refh'
            % Standard UK recommended values if not provided
            if ~isfield(uh,'Up'), uh.Up = 0.65; end
            if ~isfield(uh,'Uk'), uh.Uk = 0.8;  end
            % For hourly discrete grid, Var(T_U^(1)) ≈ 0.533, so c ≈ 0.73
            c = sqrt(0.533);

        case 'triangular'
            % ReFH triangular (Up=0.5, Uk=1): Var = 13/18
            c = sqrt(13/18);   % ≈ 0.84984

        case 'gamma'
            if ~isfield(uh,'alpha'), error('gamma UH: uh.alpha (>1) is required.'); end
            alpha = uh.alpha;
            if alpha <= 1, error('gamma UH: alpha must be > 1.'); end
            % Var(T_U^(1)) = alpha/(alpha-1)^2
            c = sqrt(alpha) / (alpha - 1);

        otherwise
            error('Unknown uh.type: %s', uh.type);
    end
end
