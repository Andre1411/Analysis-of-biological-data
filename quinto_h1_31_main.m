

%%

close all ; clear ; clc
load('test_a.mat')
Tv = 1 ;                % passo virtual grid
tv = (1:Tv:ts(end)) ;   % costruzione virtual grid
ns = length(ts) ;       % numero di campioni sulla sampling grid


%% esempi di smoothing

[u_hat_10,res_10] = bin_smoother(ys,ts,tv,10) ; % over
[u_hat_20,res_20] = bin_smoother(ys,ts,tv,20) ; % ok
[u_hat_75,res_75] = bin_smoother(ys,ts,tv,75) ; % under

figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, title('10 bin, oversmoothing')
        plot(ts,ys,'bo'),plot(tv,u_hat_10,'r','linewidth',1.25')
        legend('Dati','Bin smoother','location','best'), xlabel('t [min]')
    subplot(212), hold on, grid minor, title('10 bin, residui'), xlabel('t [min]')
        plot(ts,res_10,'bo--'), plot(ts,zeros(1,ns),'r--','linewidth',1.5)
figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, title('20 bin, smoothing ragionevole')
        plot(ts,ys,'bo'),plot(tv,u_hat_20,'r','linewidth',1.25')
        legend('Dati','Bin smoother','location','best'), xlabel('t [min]')
    subplot(212), hold on, grid minor, title('20 bin, residui'), xlabel('t [min]')
        plot(ts,res_20,'bo--'), plot(ts,zeros(1,ns),'r--','linewidth',1.5)
figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, title('75 bin, undersmoothing')
        plot(ts,ys,'bo'),plot(tv,u_hat_75,'r','linewidth',1.25')
        legend('Dati','Bin smoother','location','best'), xlabel('t [min]')
    subplot(212), hold on, grid minor, title('75 bin, residui'), xlabel('t [min]')
        plot(ts,res_75,'bo--'), plot(ts,zeros(1,ns),'r--','linewidth',1.5)
        
        
%% ricerca nbin ottimo con discrepanza

% idea: cercare il numero di bin che dà la somma dei residui assoluti più
% vicina a ns*sd^2

p_max = 200 ;
p_min = 1 ;
n_iter = 0 ;
arss = inf ;
tol = 2.5e-2 ;

% calcolo del numero di bin con criterio di discrepanza
while abs((arss-ns*sd^2)/(ns*sd^2)) > tol
    
    p = round(10^((log10(p_min)+log10(p_max))/2)) ;
    [~,res] = bin_smoother(ys,ts,tv,p) ;
    arss = res'*res ;
    n_iter = n_iter+1 ;
    disp('_____________________________________________')
    disp(['# bin = ',num2str(p)])
    disp(['Norma quadratica residui assoluti ',num2str(arss)])
    disp(['Valore ottimale ',num2str(sd^2*ns)])
    disp(['# iterazoini = ',num2str(n_iter)])
    disp('_____________________________________________')
    if arss > ns*sd^2
        p_min = p ;
    else
        p_max = p ;
    end
    
end

nbin_vec = 1:ns ;                           % vettore con numero di bin da testare
res_abs_vec = zeros(1,length(nbin_vec)) ;   % veottore con residui corrispondenti al numero di bin
% calcolo dello smoother e rispettivi residui per diversi valori del numero
% di bin, ispezionando tutti i possibili valori
for i = 1:length(nbin_vec)
    [~,res] = bin_smoother(ys,ts,tv,nbin_vec(i)) ;
    res_abs_vec(i) = res'*res ;
end
% ricerca del numero di bin con somma quadrtica dei residui assoluti più
% vicina al valore ideale ns*sd^2
[min_res_norm,idx_min] = min(abs(res_abs_vec-ns*sd^2)) ;
best_nbin = nbin_vec(idx_min) ;

figure, hold on, grid minor, title('Somma quadratica residui assoluti')
	stem(nbin_vec,res_abs_vec,'b.'), xlabel('#bin')
	plot(nbin_vec,ns*sd^2*ones(length(ts)),'r--','linewidth',1.25)
    legend('ARSS','target')

% calcolo smoother e residui con numero ottimale di bin
[u_hat,res] = bin_smoother(ys,ts,tv,best_nbin) ;

figure('units','normalized','outerposition',[0 0 0.5 1])
    subplot(211), hold on, grid minor, xlabel('t [min]')
        plot(ts,ys,'bo',tv,u_hat,'r'), legend('Dati','Smoother')
        title(['Bin smoother ottimale discrepanza, #bin = ',num2str(best_nbin)])
    subplot(212), hold on, grid minor, xlabel('t [min]')
        plot(ts,res,'bo','LineWidth',1.25)
        plot(ts,res,'b--'), plot(ts,zeros(ns),'r--','linewidth',1.25)
        title(['Residui assoluti, ARSS = ',num2str(res'*res)])
        
disp(['Valore minimo norm quad err abs:   ',num2str(res_abs_vec(idx_min)),' (#bin = ',num2str(best_nbin),')'])
disp(['Velore ottimale norm quad err abs: ',num2str(ns*sd^2)])


%%

