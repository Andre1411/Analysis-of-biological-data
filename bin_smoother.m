function [u_hat,res] = bin_smoother(ys,ts,tv,nbin)

% riceve in ingersso il vettore dei tempi (ts) e delle misure (ys) del segnale,
% la virtual grid (tv) ed il numero di bin (nbin), restituisce il valore dello
% smoother sulla vg (u_hat) e i residui assoluti sulla sampling grid (res)

nv = length(tv) ;               % calcolo del numero di camioni della vg
Tv = tv(2)-tv(1) ;              % calcollo del passo della vg
u_hat = zeros(length(tv),1) ;   % inizializzazione smoother
t_end = tv(end) ;
dt = round(t_end/nbin) ;        % calcolo dell'ampiezza temporale del bin
i = 1 ;

% ciclo di costruzione del bin
while dt*i <= t_end
    u_hat(tv<=dt*i & tv>dt*i-dt) = mean(ys(ts<=dt*i & ts>dt*i-dt)) ;
    i = i+1 ;
end

% calcolo dgli istanti sulla vg corrispondenti a quelli sulla sampling grid
idx = zeros(nv,1) ;
for i = 1:length(ts)
    k = ts(i)/Tv ;
    idx(k) = 1 ;
end

idx = logical(idx) ;
res = u_hat(idx)-ys ;           % calcolo dei residui assoluti

end