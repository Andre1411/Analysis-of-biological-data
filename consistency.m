function [u_hat,ys_hat,res] = consistency(sigma2,B,F,G,g_max,g_min,ys,ns,nv,crit)

% costrzione della matrice B^(-1/2)
B_root = diag(sqrt(diag(B)).^-1) ;
% single value decompositione della matrice 'H = B^(-1/2)*G*F^(-1)'
disp(' ')
disp('Computing SVD')
[U,D,V] = svd(B_root*G*inv(F)) ;
disp(' ')
disp('SVD complete')

convergenza = 0 ;
n_iter = 0 ;
% cambio di coordinate del vettore delle misure
psi = U'*B_root*ys ;
d = diag(D(1:ns,1:ns)) ;

switch crit
    
    case 3 % criterio di consistenza 3
        disp(' ')
        disp('criterio di consistenza 3')
        while convergenza == 0
            gamma = 10^((log10(g_min)+log10(g_max))/2) ;
            ni_hat = d.*psi./(d.^2+gamma) ;
            pred = (d.^2).*psi./(d.^2+gamma) ;
            rho = psi-pred ;
            wess = sum((d.*psi./(d.^2+gamma)).^2) ;
            wrss = sum((gamma*psi./(d.^2+gamma)).^2) ;
            q = sum((d.^2)./(d.^2+gamma)) ;
            n_iter = n_iter+1 ;
            disp('***********************')
            disp(['# iterazioni = ',num2str(n_iter),' ;'])
            disp(['WRSS = ',num2str(wrss),' ;'])
            disp(['WESS = ',num2str(wess),' ;'])
            disp(['left hand side =  ',num2str(wrss/(ns-q)),' ;'])
            disp(['right hand side = ',num2str(gamma*wess/q),' ;'])
            disp(' ')
            if wrss/(ns-q) > gamma*wess/q
                g_max = gamma ;
            else
                g_min = gamma ;
            end
            if abs(wrss/(ns-q)-gamma*wess/q) < 0.001
                convergenza = 1 ;
            end
        end
        
    case 1 % criterio di consistenza 1
        disp(' ')
        disp('criterio di consistenza 1')
        while convergenza == 0
            gamma = 10^((log10(g_min)+log10(g_max))/2) ;
            ni_hat = d.*psi./(d.^2+gamma) ;
            pred = (d.^2).*psi./(d.^2+gamma) ;
            rho = psi-pred ;
            wess = sum((d.*psi./(d.^2+gamma)).^2) ;
            q = sum((d.^2)./(d.^2+gamma)) ;
            n_iter = n_iter+1 ;
            disp('***********************')
            disp(['# iterazioni = ',num2str(n_iter),' ;'])
            disp(['WESS =   ',num2str(wess),' ;'])
            disp(['target = ',num2str(sigma2*q/gamma),' ;'])
            disp(' ')
            if wess > sigma2*q/gamma
                g_max = gamma ;
            else
                g_min = gamma ;
            end
            if abs((wess-sigma2*q/gamma)/wess) < 0.001
                convergenza = 1 ;
            end
        end

    case 2 % criterio di consistenza 2
        disp(' ')
        disp('criterio di consistenza 2')
        while convergenza == 0
            gamma = 10^((log10(g_min)+log10(g_max))/2) ;
            ni_hat = d.*psi./(d.^2+gamma) ;
            pred = (d.^2).*psi./(d.^2+gamma) ;
            rho = psi-pred ;
            wrss = sum((gamma*psi./(d.^2+gamma)).^2) ;
            q = sum((d.^2)./(d.^2+gamma)) ;
            n_iter = n_iter+1 ;
            disp('***********************')
            disp(['# iterazioni = ',num2str(n_iter),' ;'])
            disp(['WRSS   = ',num2str(wrss),' ;'])
            disp(['target = ',num2str(sigma2*(ns-q)),' ;'])
            disp(' ')
            if wrss > sigma2*(ns-q)
                g_max = gamma ;
            else
                g_min = gamma ;
            end
            if abs((wrss-sigma2*(ns-q))/wrss) < 0.001
                convergenza = 1 ;
            end
        end
        
end

% ritorno alle coordinate originali
ni_hat(end:nv) = 0 ;
u_hat = filter(1,[1 -2 1],V*ni_hat) ;
ys_hat = G*u_hat ;
res = ys_hat-ys ;

end