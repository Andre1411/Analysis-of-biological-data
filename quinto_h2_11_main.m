

%% caricamento dei dati

close all ; clear ; clc
load('data_LFP.mat') ;
ts = time ; clear time ; % [ms]
ys = data ; clear data ; % [mV]
sd = 0.005 ; % sd costante del rumore di misura

    
%% crazione matrici e griglie temporali

% creazione virtual grid
Tv = 1/Fs*1e3 ;
tv = (Tv:Tv:ts(end))' ;
nv = length(tv) ;
ns = length(ts) ;

% costruzione di F
m = 2 ;
sigma2 = sd^2 ;
B = diag(ones(1,ns)) ;
if m == 0
    F = eye(nv) ;
else
    c = zeros(nv,1) ;
    c(1) = 1 ;
    c(2) = -1 ;
    r = zeros(nv,1) ;
    r(1) = 1 ;
    delta = toeplitz(c,r) ;
    F = delta^m ;
end

% trovo i virtually missing samples
idx = zeros(nv,1) ;
for i = 1:ns
    k = round(ts(i)/Tv) ;
    idx(k) = 1 ;
end
idx = logical(idx) ;

% crazione delle matrici di trasferimento
gv = zeros(1,nv) ;
for i = 0:(nv-1)
    gv(i+1) = Tv^2*(2*i+1)/2 ;
end
% creazione di Gv e di G doppia integrazione
r = zeros(1,nv) ;
r(1) = gv(1) ;
Gv = toeplitz(gv,r) ;
G = Gv(idx,:) ;
% creazione di Gv1 e di G1 singola integrazione
gv1 = Tv*ones(1,nv) ;
r1 = zeros(1,nv) ;
r1(1) = gv1(1) ;
Gv1 = toeplitz(gv1,r1) ;
G1 = Gv1(idx,:) ;


%% calcolo delle stime

disp(' ')
crit = input('Inserire numero del criterio di consistenza da utilizzare: ') ;
[u_hat,ys_hat,res] = consistency(sigma2,B,F,G,1e10,1e-10,ys,ns,nv,crit) ;

% calcolo della stima su vg
y_reg = Gv*u_hat ;
% calcolo della derivata prima su vg
dy_dt = Gv1*u_hat ;
% calcolo della derivata seconda su vg
d2y_dt2 = u_hat ;


%% ricerca minimo, massimo, flesso e rappresentazione

L = 5/Tv ;
R = 25/Tv ;
m = round((L+R)/2) ;
    
[val_max,idx_max] = min(abs(dy_dt(L:m))) ;
[val_min,idx_min] = min(abs(dy_dt(m:R))) ;
idx_max = idx_max+L-1 ;
idx_min = idx_min+m-1 ;

figure, hold on, grid minor, title('1st derivative')
    plot(tv,dy_dt,'b',tv,zeros(1,nv),'g--',...       
        [tv(L) tv(L)],[-0.07 0.03],'c--',...
        [tv(R) tv(R)],[-0.07 0.03],'c--',...
        tv(idx_max),dy_dt(idx_max),'ro',...
        tv(idx_min),dy_dt(idx_min),'ro','linewidth',1.25)
    ylim([-0.07 0.03]), xlabel('t [ms]'), ylabel('[mV/ms]')
    legend('dy/dt','Zero','Estremi di ricerca')

L = idx_max ;
R = idx_min ;
m = round((L+R)/2) ;

[val_f,idx_f] = min(abs(d2y_dt2(L:R))) ;
idx_f = idx_f+L-1 ;

figure, hold on, grid minor, title('2nd derivative')
    plot(tv,d2y_dt2,'b',tv,zeros(1,nv),'g--',...
        [tv(L) tv(L)],[-0.02 0.02],'c--',...
        [tv(R) tv(R)],[-0.02 0.02],'c--',...
        tv(idx_f),d2y_dt2(idx_f),'ro','linewidth',1.25)
    ylim([-0.02 0.02]), xlabel('t [ms]'), ylabel('[mV/ms^2]')
    legend('d^2y/dt^2','Zero','Estremi di ricerca')
    
figure, hold on, grid minor, title('LFP reg')
    plot(tv,y_reg,'b',...
        [tv(idx_max) tv(idx_max)],[0.05 -0.4],'r--',...
        [tv(idx_min) tv(idx_min)],[0.05 -0.4],'g--',...
        [tv(idx_f) tv(idx_f)],[0.05 -0.4],'c--','linewidth',1.25)
    ylim([-0.4 0.05]), xlabel('t [ms]'), ylabel('V [mV]')
    legend('LFP regolarizzato','Massimo','Minimo','Flesso')

figure
    subplot(211), hold on, grid minor, title(['Criterio di consistenza ',num2str(crit)])
        plot(tv,y_reg,'r','linewidth',1.25), plot(ts,ys,'b')
        legend('LFP regolarizzato','LFP originale')
        xlabel('t [ms]'), ylabel('V [mV]')
    subplot(212), hold on, grid minor, title('Residui assoluti')
        area(ts,res,'facecolor','r','edgecolor','none')
        xlabel('t [ms]'), ylabel('V [mV]')


%%

figure
    subplot(211), hold on, grid minor, title(['Criterio di consistenza ',num2str(crit)])
        plot(tv,y_reg,'r','linewidth',1.25), plot(ts,ys,'b')
        legend('LFP regolarizzato','LFP originale')
        xlabel('t [ms]'), ylabel('V [mV]')
    subplot(212), hold on, grid minor, title('Residui normalizzati')
        plot(ts,res./sqrt(sigma2),'b-')
        plot(ts,zeros(1,ns),'k-',ts,ones(1,ns),'k--',ts,-ones(1,ns),'k--')
        xlabel('t [ms]'), ylabel('V [mV]')


%%

