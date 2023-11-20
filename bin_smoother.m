function [u_hat, res] = bin_smoother(ys, ts, tv, nbin)

% receives as input the time vector (ts) and signal measurements vector (ys),
% the virtual grid (tv), and the number of bins (nbin). Returns the smoother value
% on the virtual grid (u_hat) and the absolute residuals on the sampling grid (res)

nv = length(tv);               % calculation of the number of samples on the virtual grid
Tv = tv(2) - tv(1);             % calculation of the virtual grid step
u_hat = zeros(length(tv), 1);   % smoother initialization
t_end = tv(end);
dt = round(t_end/nbin);         % calculation of the temporal bin width
i = 1;

% loop for constructing the bin
while dt*i <= t_end
    u_hat(tv <= dt*i & tv > dt*i - dt) = mean(ys(ts <= dt*i & ts > dt*i - dt));
    i = i + 1;
end

% calculation of instants on the virtual grid corresponding to those on the sampling grid
idx = zeros(nv, 1);
for i = 1:length(ts)
    k = ts(i)/Tv;
    idx(k) = 1;
end

idx = logical(idx);
res = u_hat(idx) - ys;           % calculation of absolute residuals

end
