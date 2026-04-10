%%%%%%% Script for the two-dimensional example using velocity autocorrelation data

clear, close all

%%%%%% 1. Read in data
dim = 2;
datav1v1 = importdata('Data/Data2D/vv_1.txt');
datav1v2 = importdata('Data/Data2D/vv_12.txt');
datav2v2 = importdata('Data/Data2D/vv_2.txt');

datax1v1 = importdata('Data/Data2D/xv_1.txt');
datax1v2 = importdata('Data/Data2D/xv_12.txt');
datax2v2 = importdata('Data/Data2D/xv_2.txt');

datax1x1 = importdata('Data/Data2D/xx_1.txt');
datax1x2 = importdata('Data/Data2D/xx_12.txt');
datax2x2 = importdata('Data/Data2D/xx_2.txt');

range = 1:5000;
phiv(1,1,:) = datav1v1(range); phiv(1,2,:) = datav1v2(range);
phiv(2,1,:) = phiv(1,2,:); phiv(2,2,:) = datav2v2(range);

phixv(1,1,:) = datax1v1(range); phixv(1,2,:) = datax1v2(range);
phixv(2,1,:) = phixv(1,2,:); phixv(2,2,:) = datax2v2(range);

phix(1,1,:) = datax1x1(range); phix(1,2,:) = datax1x2(range);
phix(2,1,:) = phix(1,2,:); phix(2,2,:) = datax2x2(range);
 
tt = 0:0.005:5*4.999;

% normalization of vacf
L = chol(phiv(:,:,1))';
for i=1:length(tt)
    phiv(:,:,i) = L\phiv(:,:,i)/L';
    phix(:,:,i) = L\phix(:,:,i)/L';
    phixv(:,:,i) = L\phixv(:,:,i)/L';
end

%%%%%% 2. Sample data
n = 18;
n2 = 2*n;
twidth = 40;

y = phiv(:,:,1:twidth:twidth*n2);
t = tt(1:twidth:twidth*n2);
g = 100;
r = 1.3;
opt.aaatol = 5e-4;

%%%%%% 3. Compute Markovian approximation
Omega = inv(phix(:,:,1)); % Omega known
out = markovianAppPot(y,t,g,r,Omega,opt);
A = out.A; Sigma = out.Sigma;

Omega = inf; % Omega unknown
outnoOmega = markovianAppPot(y,t,g,r,Omega,opt);
A2 = outnoOmega.A; Sigma2 = outnoOmega.Sigma;

m = size(A,1);
E1 = eye(m,dim);
E2 = flip(flip(E1,1),2);

for i=1:length(tt)
    appv(:,:,i) = E1'*expm(tt(i)*A)*Sigma*E1;
    appxv(:,:,i) = E2'*expm(tt(i)*A)*Sigma*E1;
    appx(:,:,i) = E2'*expm(tt(i)*A)*Sigma*E2;
    
    appvnoOmega(:,:,i) = E1'*expm(tt(i)*A2)*Sigma2*E1;
    appxvnoOmega(:,:,i) = E2'*expm(tt(i)*A2)*Sigma2*E1;
    appxnoOmega(:,:,i) = E2'*expm(tt(i)*A2)*Sigma2*E2;
end

%%%%%%% VACF normal
pcut = 5;
cut = 15;
inset = (pcut/0.005)+1:cut/0.005;

figure
set(gcf,'Position',[100,100,1000,800]);
% V1V1
subplot(2,2,1)
plot(tt,reshape(phiv(1,1,:),[],1),'k--','LineWidth',2);
hold on
plot(t,reshape(y(1,1,:),[],1),'ko','MarkerFaceColor','k','MarkerSize',3)
plot(tt,reshape(appv(1,1,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appvnoOmega(1,1,:),[],1),'r-','LineWidth',2);
ylabel('$C_{V_1V_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)
legend('ACF','Data points','$\Omega$ known', '$\Omega$ unknown', 'Location', 'southeast',...
    'FontSize', 14,'interpreter','latex','NumColumns',2);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phiv(1,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appv(1,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appvnoOmega(1,1,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);
% V1V2
subplot(2,2,2)
plot(tt,reshape(phiv(1,2,:),[],1),'k--','LineWidth',2);
hold on
plot(t,reshape(y(1,2,:),[],1),'ko','MarkerFaceColor','k','MarkerSize',3)
plot(tt,reshape(appv(1,2,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appvnoOmega(1,2,:),[],1),'r-','LineWidth',2);
ylabel('$C_{V_1V_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phiv(1,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appv(1,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appvnoOmega(1,2,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);
% V2V1
subplot(2,2,3)
plot(tt,reshape(phiv(2,1,:),[],1),'k--','LineWidth',2);
hold on
plot(t,reshape(y(2,1,:),[],1),'ko','MarkerFaceColor','k','MarkerSize',3)
plot(tt,reshape(appv(2,1,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appvnoOmega(2,1,:),[],1),'r-','LineWidth',2);
ylabel('$C_{V_2V_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phiv(2,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appv(2,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appvnoOmega(2,1,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);
% V2V2
subplot(2,2,4)
plot(tt,reshape(phiv(2,2,:),[],1),'k--','LineWidth',2);
hold on
plot(t,reshape(y(2,2,:),[],1),'ko','MarkerFaceColor','k','MarkerSize',3)
plot(tt,reshape(appv(2,2,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appvnoOmega(2,2,:),[],1),'r-','LineWidth',2);
ylabel('$C_{V_2V_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phiv(2,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appv(2,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appvnoOmega(2,2,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

%%%%%%% PACF normal
figure
set(gcf,'Position',[100,100,1000,800]);
% X1X1
subplot(2,2,1)
plot(tt,reshape(phix(1,1,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appx(1,1,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxnoOmega(1,1,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_1R_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)
legend('ACF','$\Omega$ known', '$\Omega$ unknown', 'Location', 'south',...
    'FontSize', 14,'interpreter','latex','NumColumns',2);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phix(1,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appx(1,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxnoOmega(1,1,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

% X1X2
subplot(2,2,2)
plot(tt,reshape(phix(1,2,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appx(1,2,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxnoOmega(1,2,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_1R_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phix(1,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appx(1,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxnoOmega(1,2,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

% X2X1
subplot(2,2,3)
plot(tt,reshape(phix(2,1,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appx(2,1,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxnoOmega(2,1,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_2R_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phix(2,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appx(2,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxnoOmega(2,1,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

% X2X2
subplot(2,2,4)
plot(tt,reshape(phix(2,2,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appx(2,2,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxnoOmega(2,2,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_2R_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phix(2,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appx(2,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxnoOmega(2,2,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

%%%%%% Cross correlation normal
pcut = 5;
cut = 15;
inset = (pcut/0.005)+1:cut/0.005;

figure
set(gcf,'Position',[100,100,1000,800]);
% X1V1
subplot(2,2,1)
plot(tt,reshape(phixv(1,1,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appxv(1,1,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxvnoOmega(1,1,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_1V_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)
legend('ACF','$\Omega$ known', '$\Omega$ unknown', 'Location', 'south',...
    'FontSize', 14,'interpreter','latex','NumColumns',2);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phixv(1,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appxv(1,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxvnoOmega(1,1,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);
% X1V2
subplot(2,2,2)
plot(tt,reshape(phixv(1,2,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appxv(1,2,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxvnoOmega(1,2,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_1V_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phixv(1,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appxv(1,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxvnoOmega(1,2,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);
% X2V1
subplot(2,2,3)
plot(tt,reshape(phixv(2,1,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appxv(2,1,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxvnoOmega(2,1,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_2V_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phixv(2,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appxv(2,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxvnoOmega(2,1,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);
% X2V2
subplot(2,2,4)
plot(tt,reshape(phixv(2,2,:),[],1),'k--','LineWidth',2);
hold on
plot(tt,reshape(appxv(2,2,:),[],1),'b-','LineWidth',2);
plot(tt,reshape(appxvnoOmega(2,2,:),[],1),'r-','LineWidth',2);
ylabel('$C_{R_2V_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(phixv(2,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(appxv(2,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxvnoOmega(2,2,inset),[],1), 'r-','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

%%%%%% Error Fro-norm
errv = zeros(n,1); errxv = errv; errx = errv;
errvnok = errv; errxvnok = errv; errxnok = errv;
for i=1:length(tt)
    errv(i) = norm(appv(:,:,i)-phiv(:,:,i),'fro');
    errxv(i) = norm(appxv(:,:,i)-phixv(:,:,i),'fro');
    errx(i) = norm(appx(:,:,i)-phix(:,:,i),'fro');
    errvnok(i) = norm(appvnoOmega(:,:,i)-phiv(:,:,i),'fro');
    errxvnok(i) = norm(appxvnoOmega(:,:,i)-phixv(:,:,i),'fro');
    errxnok(i) = norm(appxnoOmega(:,:,i)-phix(:,:,i),'fro');
end

figure
set(gcf,'Position',[100,100,1200,600]);
subplot(1,2,1)
semilogy(tt,errv,'LineWidth',2)
hold on
semilogy(tt,errxv,'LineWidth',2)
semilogy(tt,errx,'LineWidth',2)
xlim([0,5]);
xlabel('t')
ylabel('$\Vert\Delta C(t) \Vert_F$','Interpreter','latex')
grid on
legend('V','RV','R','Location','southeast')
set(gca,'FontSize',16)
title('\Omega known')

subplot(1,2,2)
semilogy(tt,errvnok,'LineWidth',2)
hold on
semilogy(tt,errxvnok,'LineWidth',2)
semilogy(tt,errxnok,'LineWidth',2)
xlim([0,5]);
ylim([1e-4,1e-1]);
xlabel('t')
ylabel('$\Vert\Delta C(t) \Vert_F$','Interpreter','latex')
grid on
legend('V','RV','R','Location','southeast')
set(gca,'FontSize',16)
title('\Omega unknown')