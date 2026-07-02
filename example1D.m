%%%%%%% Script for the one-dimensional case using velocity autocorrelation data

clear, close all

%%%%%%% 1. Read in data
dim = 1;
datavv = importdata('Data/Data1D/vv.txt');
dataxv = importdata('Data/Data1D/xv.txt');
dataxx = importdata('Data/Data1D/xx.txt');

datavv = datavv/datavv(1);
dataxx = dataxx/datavv(1);
dataxv = dataxv/datavv(1);

phiv = reshape(datavv,1,1,[]);
phix = reshape(dataxx,1,1,[]);
phixv = reshape(dataxv,1,1,[]);

tt = 0:0.005:5*4.999;

%%%%%% 2. Sample data
n = 41;
twidth = 5;

y = phiv(:,:,1:twidth:twidth*n);
t = tt(1:twidth:twidth*n);
g = 100;
r = 1.15;
opt.aaatol = 1e-4; % specify aaa tolerance (standard: 1e-4)
opt.acf = 'v'; % specify acf data used ('p' or 'v', standard: v)

%%%%%% 3. Compute Markovian approximation
% Omega known
Omega = inv(phix(:,:,1));
out = markovianAppPot(y,t,g,r,Omega,opt);
A = out.A; Sigma = out.Sigma;

% Omega unknown
Omega = inf;
outnoOmega = markovianAppPot(y,t,g,r,Omega,opt);
A2 = outnoOmega.A; Sigma2 = outnoOmega.Sigma;

% PACF data used
n = 31;
twidth = 5;

Omega = inv(phix(:,:,1));
opt.acf = 'p'; opt.aaatol = 1e-6;
r = 1.2;
yp = phix(:,:,1:twidth:twidth*n);
t = tt(1:twidth:twidth*n);
outp = markovianAppPot(yp,t,g,r,Omega,opt);
A3 = outp.A; Sigma3 = outp.Sigma;

m = size(A,1); m2 = size(A3,1);
E1 = eye(m,1); E2 = flip(E1); E1p = eye(m2,1); E2p = flip(E1p);
for i=1:length(tt)
    appv(i,:) = E1'*expm(tt(i)*A)*Sigma*E1;
    appvnoOmega(i,:) = E1'*expm(tt(i)*A2)*Sigma2*E1;
    appvp(i,:) = E1p'*expm(tt(i)*A3)*Sigma3*E1p;
    
    appx(i,:) = E2'*expm(tt(i)*A)*Sigma*E2;
    appxnoOmega(i,:) = E2'*expm(tt(i)*A2)*Sigma2*E2;
    appxp(i,:) = E2p'*expm(tt(i)*A3)*Sigma3*E2p;
    
    appxv(i,:) = E2'*expm(tt(i)*A)*Sigma*E1;
    appxvnoOmega(i,:) = E2'*expm(tt(i)*A2)*Sigma2*E1;
    appxvp(i,:) = E2p'*expm(tt(i)*A3)*Sigma3*E1p;
end

phiv = phiv(:); phix = phix(:); phixv = phixv(:); y = y(:);

cut = 500; inset = 150:500; prange = 1:cut;
yellow = [0.929,0.694,0.125];

%%% Normal plots
figure
set(gcf, 'Position',  [100, 100, 1000, 800])
subplot(2,2,1)
p1 = plot(tt(prange),appvp(prange),'b-','LineWidth',2);
hold on
p2 = plot(tt(prange),appv(prange),'-','color',yellow,'LineWidth',2);
p3 = plot(tt(prange),appvnoOmega(prange),'r-','LineWidth',2);
p4 = plot(tt(prange),phiv(prange),'k--','LineWidth',2);
ylabel('$C_{V}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)
legend([p2,p3,p1,p4],'$\Omega$ known', '$\Omega$ unknown','PACF data','ACF','Location', 'southeast',...
    'FontSize', 12,'interpreter','latex','NumColumns',2);

subplot(2,2,3)
plot(tt(prange),appxvp(prange),'b-','LineWidth',2);
hold on
plot(tt(prange),appxv(prange),'-','color',yellow,'LineWidth',2);
plot(tt(prange),appxvnoOmega(prange),'r-','LineWidth',2);
plot(tt(prange),phixv(prange),'k--','LineWidth',2);
ylabel('$C_{RV}$','interpreter', 'latex','FontSize',16);
xlabel('t','FontSize',16)
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)

subplot(2,2,2)
plot(tt(prange),-appxvp(prange),'b-','LineWidth',2);
hold on
plot(tt(prange),-appxv(prange),'-','color',yellow,'LineWidth',2);
plot(tt(prange),-appxvnoOmega(prange),'r-','LineWidth',2);
plot(tt(prange),-phixv(prange),'k--','LineWidth',2);
ylabel('$C_{VR}$','interpreter', 'latex','FontSize',16);
xlabel('t','FontSize',16)
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)

subplot(2,2,4)
plot(tt(prange),appxp(prange),'b-','LineWidth',2);
hold on
plot(tt(prange),appx(prange),'-','color',yellow,'LineWidth',2);
plot(tt(prange),appxnoOmega(prange),'r-','LineWidth',2);
plot(tt(prange),phix(prange),'k--','LineWidth',2);
ylabel('$C_{R}$','interpreter', 'latex','FontSize',16);
xlabel('t','FontSize',16)
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)

%%% Error plots
yticks_val = 10.^(-6:2:0);
ytick_labels = cellstr(num2str(round(log10(yticks_val(:))), '10^{%d}'));

figure
set(gcf, 'Position',  [100, 100, 1000, 800])
subplot(2,2,1);
p1 = semilogy(tt(prange),abs(appvp(prange)-phiv(prange)),'b-','LineWidth',2);
hold on
p2 = semilogy(tt(prange),abs(appv(prange)-phiv(prange)),'-','color',yellow,'LineWidth',2);
p3 = semilogy(tt(prange),abs(appvnoOmega(prange)-phiv(prange)),'r-','LineWidth',2);
p4 = semilogy(tt(prange),abs(phiv(prange)),'k--','LineWidth',2);
grid on;
xlim([0,2.5]);
ylim([1e-7,1e0]); yticks(yticks_val); yticklabels(ytick_labels);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{V}$','interpreter','latex','FontSize',16);
legend([p2,p3,p1,p4],'$\Omega$ known', '$\Omega$ unknown','PACF data','ACF', 'Location', 'northeast',...
    'FontSize', 12, 'interpreter', 'latex')
set(gca,'FontSize',14);

subplot(2,2,2);
semilogy(tt(prange),abs(-appxvp(prange)+phixv(prange)),'b-','LineWidth',2);
hold on
semilogy(tt(prange),abs(-appxv(prange)+phixv(prange)),'-','color',yellow,'LineWidth',2);
semilogy(tt(prange),abs(-appxvnoOmega(prange)+phixv(prange)),'r','LineWidth',2);
semilogy(tt(prange),abs(phixv(prange)),'k--','LineWidth',2);
grid on;
xlim([0,2.5]);
ylim([1e-7,1e0]); yticks(yticks_val); yticklabels(ytick_labels);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{VR}$','interpreter','latex','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
semilogy(tt(prange),abs(appxvp(prange)-phixv(prange)),'b-','LineWidth',2);
hold on
semilogy(tt(prange),abs(appxv(prange)-phixv(prange)),'-','color',yellow,'LineWidth',2);
semilogy(tt(prange),abs(appxvnoOmega(prange)-phixv(prange)),'r','LineWidth',2);
semilogy(tt(prange),abs(-phixv(prange)),'k--','LineWidth',2);
grid on;
xlim([0,2.5]);
ylim([1e-7,1e0]); yticks(yticks_val); yticklabels(ytick_labels);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{RV}$','interpreter','latex','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
semilogy(tt(prange),abs(appxp(prange)-phix(prange)),'b-','LineWidth',2);
hold on
semilogy(tt(prange),abs(appx(prange)-phix(prange)),'-','color',yellow,'LineWidth',2);
semilogy(tt(prange),abs(appxnoOmega(prange)-phix(prange)),'r','LineWidth',2);
semilogy(tt(prange),abs(phix(prange)),'k--','LineWidth',2);
grid on;
xlim([0,2.5]);
ylim([1e-7,1e0]); yticks(yticks_val); yticklabels(ytick_labels);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{R}$','interpreter','latex','FontSize',16);
set(gca,'FontSize',14);