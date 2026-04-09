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
n = 20;
n2 = 2*n;
twidth = 5;

y = phiv(:,:,1:twidth:twidth*n2);
t = tt(1:twidth:twidth*n2);
g = 100;
r = 1.15;
opt.aaatol = 1e-4;

%%%%%% 3. Compute Markovian approximation
Omega = inv(phix(:,:,1));
out = markovianAppPot(y,t,g,r,Omega,opt);
A = out.A; Sigma = out.Sigma;

Omega = inf;
outnoOmega = markovianAppPot(y,t,g,r,Omega,opt);
A2 = outnoOmega.A; Sigma2 = outnoOmega.Sigma;

m = size(A,1);
E1 = eye(m,1); E2 = flip(E1);
for i=1:length(tt)
    appv(i,:) = E1'*expm(tt(i)*A)*Sigma*E1;
    appvnoOmega(i,:) = E1'*expm(tt(i)*A2)*Sigma2*E1;
    appx(i,:) = E2'*expm(tt(i)*A)*Sigma*E2;
    appxnoOmega(i,:) = E2'*expm(tt(i)*A2)*Sigma2*E2;
    appxv(i,:) = E2'*expm(tt(i)*A)*Sigma*E1;
    appxvnoOmega(i,:) = E2'*expm(tt(i)*A2)*Sigma2*E1;
end

appvx = -appxv; appvxnoOmega = -appxvnoOmega;
phiv = phiv(:); phix = phix(:); phixv = phixv(:); y = y(:);

cut = 500; inset = 200:cut; prange = 1:cut;

%%% Normal plots
figure
set(gcf, 'Position',  [100, 100, 1000, 800])
subplot(2,2,1)
plot(tt(prange),phiv(prange),'k--','LineWidth',2);
hold on
plot(t,y,'ko','MarkerFaceColor','k','MarkerSize',3);
plot(tt(prange),appv(prange),'b-','LineWidth',2);
plot(tt(prange),appvnoOmega(prange),'r-','LineWidth',2);
ylabel('$C_{VV}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)
legend('ACF','Data points','$\Omega$ known', '$\Omega$ unknown', 'Location', 'southeast',...
    'FontSize', 12,'interpreter','latex','NumColumns',2);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.15, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), phiv(inset), 'k--','LineWidth',2);
plot(h, tt(inset), appv(inset), 'b-','LineWidth',2);
plot(h, tt(inset), appvnoOmega(inset), 'r-','LineWidth',2);
grid on
set(gca,'FontSize',11);

subplot(2,2,3)
plot(tt(prange),phixv(prange),'k--','LineWidth',2);
hold on
plot(tt(prange),appxv(prange),'b-','LineWidth',2);
plot(tt(prange),appxvnoOmega(prange),'r-','LineWidth',2);
ylabel('$C_{VR}$','interpreter', 'latex','FontSize',16);
xlabel('t','FontSize',16)
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.12, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), phixv(inset), 'k--','LineWidth',2);
plot(h, tt(inset), appxv(inset), 'b-','LineWidth',2);
plot(h, tt(inset), appxvnoOmega(inset), 'r-','LineWidth',2);
grid on
set(gca,'FontSize',11);

subplot(2,2,2)
plot(tt(prange),-phixv(prange),'k--','LineWidth',2);
hold on
plot(tt(prange),appvx(prange),'b-','LineWidth',2);
plot(tt(prange),appvxnoOmega(prange),'r-','LineWidth',2);
ylabel('$C_{RV}$','interpreter', 'latex','FontSize',16);
xlabel('t','FontSize',16)
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.1, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), -phixv(inset), 'k--','LineWidth',2);
plot(h, tt(inset), -appxv(inset), 'b-','LineWidth',2);
plot(h, tt(inset), -appxvnoOmega(inset), 'r-','LineWidth',2);
grid on
set(gca,'FontSize',11);

subplot(2,2,4)
plot(tt(prange),phix(prange),'k--','LineWidth',2);
hold on
plot(tt(prange),appx(prange),'b-','LineWidth',2);
plot(tt(prange),appxnoOmega(prange),'r-','LineWidth',2);
ylabel('$C_{RR}$','interpreter', 'latex','FontSize',16);
xlabel('t','FontSize',16)
xlim([0,2.5]);
grid on;
set(gca,'FontSize',14)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.12, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), phix(inset), 'k--','LineWidth',2);
plot(h, tt(inset), appx(inset), 'b-','LineWidth',2);
plot(h, tt(inset), appvnoOmega(inset), 'r-','LineWidth',2);
grid on
set(gca,'FontSize',11);

%%% Error plots
figure
set(gcf, 'Position',  [100, 100, 1000, 800])
subplot(2,2,1);
semilogy(tt(prange),abs(phiv(prange)),'k-','LineWidth',2);
hold on
semilogy(tt(prange),abs(appv(prange)-phiv(prange)),'b-','LineWidth',2);
semilogy(tt(prange),abs(appvnoOmega(prange)-phiv(prange)),'r-','LineWidth',2)
grid on;
xlim([0,2.5]);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{VV}$','interpreter','latex','FontSize',16);
legend('ACF','$\Omega$ known', '$\Omega$ unknown', 'Location', 'northeast',...
    'FontSize', 12, 'interpreter', 'latex')
set(gca,'FontSize',14);

subplot(2,2,2);
semilogy(tt(prange),abs(phixv(prange)),'k-','LineWidth',2);
hold on
semilogy(tt(prange),abs(appvx(prange)+phixv(prange)),'b-','LineWidth',2);
semilogy(tt(prange),abs(appvxnoOmega(prange)+phixv(prange)),'r','LineWidth',2);
grid on;
xlim([0,2.5]);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{VR}$','interpreter','latex','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
semilogy(tt(prange),abs(-phixv(prange)),'k-','LineWidth',2);
hold on
semilogy(tt(prange),abs(appxv(prange)-phixv(prange)),'b-','LineWidth',2);
semilogy(tt(prange),abs(appxvnoOmega(prange)-phixv(prange)),'r','LineWidth',2);
grid on;
xlim([0,2.5]);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{RV}$','interpreter','latex','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
semilogy(tt(prange),abs(phix(prange)),'k-','LineWidth',2);
hold on
semilogy(tt(prange),abs(appx(prange)-phix(prange)),'b-','LineWidth',2);
semilogy(tt(prange),abs(appxnoOmega(prange)-phix(prange)),'r','LineWidth',2);
grid on;
xlim([0,2.5]);
xlabel('t','FontSize',16);
ylabel('$\Delta C_{RR}$','interpreter','latex','FontSize',16);
set(gca,'FontSize',14);
