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
 
tt = 0:0.0025:4999*0.0025;

% normalization of vacf
L = chol(phiv(:,:,1))';
for i=1:length(tt)
    phivnorm(:,:,i) = L\phiv(:,:,i)/L';
    phixnorm(:,:,i) = L\phix(:,:,i)/L';
    phixvnorm(:,:,i) = L\phixv(:,:,i)/L';
end

%%%%%% 2. Sample data
n = 31;
twidth = 40;

y = phivnorm(:,:,1:twidth:twidth*n);
t = tt(1:twidth:twidth*n);
g = 100;
r = 1.5;
opt.aaatol = 5e-4;

%%%%%% 3. Compute Markovian approximation
Omega = inv(phixnorm(:,:,1)); % Omega known
out = markovianAppPot(y,t,g,r,Omega,opt);
A = out.A; Sigma = out.Sigma;

Omega = inf; % Omega unknown
outnoOmega = markovianAppPot(y,t,g,r,Omega,opt);
A2 = outnoOmega.A; Sigma2 = outnoOmega.Sigma;

m = size(A,1);
E1 = eye(m,dim)*L';
E2 = flip(flip(eye(m,dim),1),2)*L';

for i=1:length(tt)
    appv(:,:,i) = E1'*expm(tt(i)*A)*Sigma*E1;
    appxv(:,:,i) = E2'*expm(tt(i)*A)*Sigma*E1;
    appx(:,:,i) = E2'*expm(tt(i)*A)*Sigma*E2;
    
    appvnoOmega(:,:,i) = E1'*expm(tt(i)*A2)*Sigma2*E1;
    appxvnoOmega(:,:,i) = E2'*expm(tt(i)*A2)*Sigma2*E1;
    appxnoOmega(:,:,i) = E2'*expm(tt(i)*A2)*Sigma2*E2;
end

%%%%%%% VACF, PACF diagonal entries
pcut = 4;
cut = 8;
inset = (pcut/0.0025)+1:cut/0.0025;

figure
set(gcf,'Position',[100,100,1000,800]);
% V1V1
subplot(2,2,1)
plot(tt,reshape(appv(1,1,:),[],1),'b-','LineWidth',2);
hold on
plot(tt,reshape(appvnoOmega(1,1,:),[],1),'r-','LineWidth',2);
plot(tt,reshape(phiv(1,1,:),[],1),'k--','LineWidth',2);
ylabel('$C_{V_1V_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)
legend('$\Omega$ known', '$\Omega$ unknown', 'ACF', 'Location', 'southeast',...
    'FontSize', 14,'interpreter','latex','NumColumns',2);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(appv(1,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appvnoOmega(1,1,inset),[],1), 'r-','LineWidth',2);
plot(h, tt(inset), reshape(phiv(1,1,inset),[],1), 'k--','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

% V2V2
subplot(2,2,2)
plot(tt,reshape(appv(2,2,:),[],1),'b-','LineWidth',2);
hold on
plot(tt,reshape(appvnoOmega(2,2,:),[],1),'r-','LineWidth',2);
plot(tt,reshape(phiv(2,2,:),[],1),'k--','LineWidth',2);
ylabel('$C_{V_2V_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(appv(2,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appvnoOmega(2,2,inset),[],1), 'r-','LineWidth',2);
plot(h, tt(inset), reshape(phiv(2,2,inset),[],1), 'k--','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

% R1R1
subplot(2,2,3)
plot(tt,reshape(appx(1,1,:),[],1),'b-','LineWidth',2);
hold on
plot(tt,reshape(appxnoOmega(1,1,:),[],1),'r-','LineWidth',2);
plot(tt,reshape(phix(1,1,:),[],1),'k--','LineWidth',2);
ylabel('$C_{R_1R_1}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(appx(1,1,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxnoOmega(1,1,inset),[],1), 'r-','LineWidth',2);
plot(h, tt(inset), reshape(phix(1,1,inset),[],1), 'k--','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

% R2R2
subplot(2,2,4)
plot(tt,reshape(appx(2,2,:),[],1),'b-','LineWidth',2);
hold on
plot(tt,reshape(appxnoOmega(2,2,:),[],1),'r-','LineWidth',2);
plot(tt,reshape(phix(2,2,:),[],1),'k--','LineWidth',2);
ylabel('$C_{R_2R_2}$','interpreter','latex','FontSize',16);
xlabel('t','FontSize',16);
xlim([0,cut]);
grid on;
set(gca,'FontSize',18)

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.2]);
hold(h);
plot(h, tt(inset), reshape(appx(2,2,inset),[],1), 'b-','LineWidth',2);
plot(h, tt(inset), reshape(appxnoOmega(2,2,inset),[],1), 'r-','LineWidth',2);
plot(h, tt(inset), reshape(phix(2,2,inset),[],1), 'k--','LineWidth',2);
xlim([pcut,cut])
grid on
set(gca,'FontSize',11);

%%%%%% Error Fro-norm
errv = zeros(n,1); errxv = errv; errx = errv;
errvnok = errv; errxvnok = errv; errxnok = errv;
phivfro = zeros(n,1);
for i=1:length(tt)
    errv(i) = norm(appv(:,:,i)-phiv(:,:,i),'fro');
    errxv(i) = norm(appxv(:,:,i)-phixv(:,:,i),'fro');
    errx(i) = norm(appx(:,:,i)-phix(:,:,i),'fro');
    phivfro(i) = norm(phiv(:,:,i),'fro');
    errvnok(i) = norm(appvnoOmega(:,:,i)-phiv(:,:,i),'fro');
    errxvnok(i) = norm(appxvnoOmega(:,:,i)-phixv(:,:,i),'fro');
    errxnok(i) = norm(appxnoOmega(:,:,i)-phix(:,:,i),'fro');
end

figure
set(gcf,'Position',[100,100,1000,400]);
subplot(1,2,1)
semilogy(tt,errv,'LineWidth',2)
hold on
semilogy(tt,errxv,'LineWidth',2)
semilogy(tt,errx,'LineWidth',2)
semilogy(tt,phivfro,'k--','LineWidth',2)
xlim([0,10]);
ylim([1e-5,1e0])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel('t')
ylabel('$\Vert\Delta C(t) \Vert_F$','Interpreter','latex')
grid on
legend('V','RV','R', 'VACF','Location','northeast')
set(gca,'FontSize',16)
title('\Omega known')

subplot(1,2,2)
semilogy(tt,errvnok,'LineWidth',2)
hold on
semilogy(tt,errxvnok,'LineWidth',2)
semilogy(tt,errxnok,'LineWidth',2)
semilogy(tt,phivfro,'k--','LineWidth',2)
xlim([0,10]);
ylim([1e-5,1e0]);
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel('t')
ylabel('$\Vert\Delta C(t) \Vert_F$','Interpreter','latex')
grid on
legend('V','RV','R','VACF', 'Location','northeast')
set(gca,'FontSize',16)
title('\Omega unknown')
