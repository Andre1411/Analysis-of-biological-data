function [u_hat, ys_hat, res] = consistency(sigma2, B, F, G, g_max, g_min, ys, ns, nv, crit)

% Construction of the matrix B^(-1/2)
B_root = diag(sqrt(diag(B)).^(-1));
% Singular Value Decomposition of the matrix 'H = B^(-1/2)*G*F^(-1)'
disp(' ')
disp('Computing SVD')
[U, D, V] = svd(B_root * G * inv(F));
disp(' ')
disp('SVD complete')

convergence = 0;
n_iter = 0;
% Change of coordinates for the measurement vector
psi = U' * B_root * ys;
d = diag(D(1:ns, 1:ns));

switch crit

    case 3 % Consistency criterion 3
        disp(' ')
        disp('Consistency criterion 3')
        while convergence == 0
            gamma = 10^((log10(g_min) + log10(g_max)) / 2);
            ni_hat = d .* psi ./ (d.^2 + gamma);
            pred = (d.^2) .* psi ./ (d.^2 + gamma);
            rho = psi - pred;
            wess = sum((d .* psi ./ (d.^2 + gamma)).^2);
            wrss = sum((gamma * psi ./ (d.^2 + gamma)).^2);
            q = sum((d.^2) ./ (d.^2 + gamma));
            n_iter = n_iter + 1;
            disp('***********************')
            disp(['# iterations = ', num2str(n_iter), ' ;'])
            disp(['WRSS = ', num2str(wrss), ' ;'])
            disp(['WESS = ', num2str(wess), ' ;'])
            disp(['left-hand side =  ', num2str(wrss / (ns - q)), ' ;'])
            disp(['right-hand side = ', num2str(gamma * wess / q), ' ;'])
            disp(' ')
            if wrss / (ns - q) > gamma * wess / q
                g_max = gamma;
            else
                g_min = gamma;
            end
            if abs(wrss / (ns - q) - gamma * wess / q) < 0.001
                convergence = 1;
            end
        end

    case 1 % Consistency criterion 1
        disp(' ')
        disp('Consistency criterion 1')
        while convergence == 0
            gamma = 10^((log10(g_min) + log10(g_max)) / 2);
            ni_hat = d .* psi ./ (d.^2 + gamma);
            pred = (d.^2) .* psi ./ (d.^2 + gamma);
            rho = psi - pred;
            wess = sum((d .* psi ./ (d.^2 + gamma)).^2);
            q = sum((d.^2) ./ (d.^2 + gamma));
            n_iter = n_iter + 1;
            disp('***********************')
            disp(['# iterations = ', num2str(n_iter), ' ;'])
            disp(['WESS =   ', num2str(wess), ' ;'])
            disp(['target = ', num2str(sigma2 * q / gamma), ' ;'])
            disp(' ')
            if wess > sigma2 * q / gamma
                g_max = gamma;
            else
                g_min = gamma;
            end
            if abs((wess - sigma2 * q / gamma) / wess) < 0.001
                convergence = 1;
            end
        end

    case 2 % Consistency criterion 2
        disp(' ')
        disp('Consistency criterion 2')
        while convergence == 0
            gamma = 10^((log10(g_min) + log10(g_max)) / 2);
            ni_hat = d .* psi ./ (d.^2 + gamma);
            pred = (d.^2) .* psi ./ (d.^2 + gamma);
            rho = psi - pred;
            wrss = sum((gamma * psi ./ (d.^2 + gamma)).^2);
            q = sum((d.^2) ./ (d.^2 + gamma));
            n_iter = n_iter + 1;
            disp('***********************')
            disp(['# iterations = ', num2str(n_iter), ' ;'])
            disp(['WRSS   = ', num2str(wrss), ' ;'])
            disp(['target = ', num2str(sigma2 * (ns - q)), ' ;'])
            disp(' ')
            if wrss > sigma2 * (ns - q)
                g_max = gamma;
            else
                g_min = gamma;
            end
            if abs((wrss - sigma2 * (ns - q)) / wrss) < 0.001
                convergence = 1;
            end
        end

end

% Return to the original coordinates
ni_hat(end:nv) = 0;
u_hat = filter(1, [1 -2 1], V * ni_hat);
ys_hat = G * u_hat;
res = ys_hat - ys;

end
