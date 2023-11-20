%%

close all; clear; clc
load('test_a.mat')
Tv = 1;                % virtual grid step
tv = (1:Tv:ts(end));   % construction of virtual grid
ns = length(ts);       % number of samples on the sampling grid


%% examples of smoothing

[u_hat_10, res_10] = bin_smoother(ys, ts, tv, 10); % over
[u_hat_20, res_20] = bin_smoother(ys, ts, tv, 20); % ok
[u_hat_75, res_75] = bin_smoother(ys, ts, tv, 75); % under

figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, title('10 bins, oversmoothing')
        plot(ts, ys, 'bo'), plot(tv, u_hat_10, 'r','linewidth',1.25')
        legend('Data', 'Bin smoother', 'location', 'best'), xlabel('t [min]')
    subplot(212), hold on, grid minor, title('10 bins, residuals'), xlabel('t [min]')
        plot(ts, res_10, 'bo--'), plot(ts, zeros(1,ns), 'r--','linewidth',1.5)
figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, title('20 bins, reasonable smoothing')
        plot(ts, ys, 'bo'), plot(tv, u_hat_20, 'r','linewidth',1.25')
        legend('Data', 'Bin smoother', 'location', 'best'), xlabel('t [min]')
    subplot(212), hold on, grid minor, title('20 bins, residuals'), xlabel('t [min]')
        plot(ts, res_20, 'bo--'), plot(ts, zeros(1,ns), 'r--','linewidth',1.5)
figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, title('75 bins, undersmoothing')
        plot(ts, ys, 'bo'), plot(tv, u_hat_75, 'r','linewidth',1.25')
        legend('Data', 'Bin smoother', 'location', 'best'), xlabel('t [min]')
    subplot(212), hold on, grid minor, title('75 bins, residuals'), xlabel('t [min]')
        plot(ts, res_75, 'bo--'), plot(ts, zeros(1,ns), 'r--','linewidth',1.5)
        
        
%% search for optimal nbin with discrepancy

% idea: find the bin number that gives the sum of absolute residuals closest to ns*sd^2

p_max = 200;
p_min = 1;
n_iter = 0;
arss = inf;
tol = 2.5e-2;

% calculation of bin number with discrepancy criterion
while abs((arss-ns*sd^2)/(ns*sd^2)) > tol
    
    p = round(10^((log10(p_min)+log10(p_max))/2));
    [~, res] = bin_smoother(ys, ts, tv, p);
    arss = res'*res;
    n_iter = n_iter+1;
    disp('_____________________________________________')
    disp(['# bins = ', num2str(p)])
    disp(['Quadratic norm of absolute residuals ', num2str(arss)])
    disp(['Optimal value ', num2str(sd^2*ns)])
    disp(['# iterations = ', num2str(n_iter)])
    disp('_____________________________________________')
    if arss > ns*sd^2
        p_min = p;
    else
        p_max = p;
    end
    
end

nbin_vec = 1:ns;                           % vector with bin numbers to test
res_abs_vec = zeros(1, length(nbin_vec));   % vector with corresponding residuals for bin numbers
% calculation of smoother and respective residuals for different bin number values,
% inspecting all possible values
for i = 1:length(nbin_vec)
    [~, res] = bin_smoother(ys, ts, tv, nbin_vec(i));
    res_abs_vec(i) = res'*res;
end
% search for bin number with sum of squared absolute residuals closest to the ideal value ns*sd^2
[min_res_norm, idx_min] = min(abs(res_abs_vec-ns*sd^2));
best_nbin = nbin_vec(idx_min);

figure, hold on, grid minor, title('Sum of squared absolute residuals')
    stem(nbin_vec, res_abs_vec, 'b.'), xlabel('# bins')
    plot(nbin_vec, ns*sd^2*ones(length(ts)), 'r--', 'linewidth', 1.25)
    legend('ARSS', 'target')

% calculation of smoother and residuals with optimal bin number
[u_hat, res] = bin_smoother(ys, ts, tv, best_nbin);

figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, xlabel('t [min]')
        plot(ts, ys, 'bo', tv, u_hat, 'r'), legend('Data', 'Smoother')
        title(['Optimal bin smoother discrepancy, # bins = ', num2str(best_nbin)])
    subplot(212), hold on, grid minor, xlabel('t [min]')
        plot(ts, res, 'bo', 'LineWidth', 1.25)
        plot(ts, res, 'b--'), plot(ts, zeros(ns), 'r--', 'linewidth', 1.25)
        title(['Absolute residuals, ARSS = ', num2str(res'*res)])
        
disp(['Minimum value of norm quad absolute error:   ', num2str(res_abs_vec(idx_min)), ' (# bins = ', num2str(best_nbin), ')'])
disp(['Optimal value of norm quad absolute error: ', num2str(ns*sd^2)])


%%
