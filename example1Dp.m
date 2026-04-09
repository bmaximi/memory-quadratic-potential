%%%%%%% Script for the one-dimensional case using position autocorrelation data

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
K = inf;

%%%%%% 2. Sample data
n = 25;
n2 = 2*n;
twidth = 2;

y = phix(:,:,1:twidth:twidth*n2);
t = tt(1:twidth:twidth*n2);
g = 100;
r = 1.1;
opt.aaatol = 1e-4; % specify aaa tolerance (standard: 1e-4)
opt.acf = 'p'; % specify acf data used ('p' or 'v', standard: v)

%%%%%% 3. Compute Markovian approximation
out = markovianAppPot(y,t,g,r,K,opt);
A = out.A; Sigma = out.Sigma;

app = zeros(2,2,length(tt));
phi = [phiv, -phixv; phixv, phix];

m = size(A,1);
E = [eye(m,1), flip(eye(m,1))];
for i=1:length(tt)
    app(:,:,i) = E'*expm(tt(i)*A)*Sigma*E;
end

cut = 1000; inset = 400:cut;

figure
set(gcf,'Position',[100 100 1000 800]);
subplot(2,2,1)
plot(tt(1:cut),reshape(phi(1,1,1:cut),[],1),'k--','LineWidth',2);
hold on
plot(tt(1:cut),reshape(app(1,1,1:cut),[],1),'r-','LineWidth',2);
xlabel('t');
ylabel('C_{VV}');
xlim([0,5]);
grid on
set(gca,'FontSize',18);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.21]);
hold(h);
plot(h, tt(inset), reshape(phi(1,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(app(1,1,inset),[],1), 'r-','LineWidth',2);
xlim([2,5])
grid on
set(gca,'FontSize',11);

subplot(2,2,2)
plot(tt(1:cut),reshape(phi(1,2,1:cut),[],1),'k--','LineWidth',2);
hold on
plot(tt(1:cut),reshape(app(1,2,1:cut),[],1),'r-','LineWidth',2);
xlabel('t');
ylabel('C_{VR}');
xlim([0,5]);
grid on
set(gca,'FontSize',18);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.21]);
hold(h);
plot(h, tt(inset), reshape(phi(1,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(app(1,2,inset),[],1), 'r-','LineWidth',2);
xlim([2,5])
grid on
set(gca,'FontSize',11);

subplot(2,2,3)
plot(tt(1:cut),reshape(phi(2,1,1:cut),[],1),'k--','LineWidth',2);
hold on
plot(tt(1:cut),reshape(app(2,1,1:cut),[],1),'r-','LineWidth',2);
xlabel('t');
ylabel('C_{RV}');
xlim([0,5]);
grid on
set(gca,'FontSize',18);

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.21]);
hold(h);
plot(h, tt(inset), reshape(phi(2,1,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(app(2,1,inset),[],1), 'r-','LineWidth',2);
xlim([2,5])
grid on
set(gca,'FontSize',11);

subplot(2,2,4)
plot(tt(1:cut),reshape(phi(2,2,1:cut),[],1),'k--','LineWidth',2);
hold on
plot(t,y(:),'ko','MarkerFaceColor','k','MarkerSize',3);
plot(tt(1:cut),reshape(phi(2,2,1:cut),[],1),'r-','LineWidth',2);
xlabel('t');
ylabel('C_{RR}');
xlim([0,5]);
grid on
set(gca,'FontSize',18);
legend('ACF', 'Data points', 'Approximation','FontSize',14,'Location','southeast')

p = get(gca,'Position');
h = axes('Parent',gcf,'Position',[p(1)+0.13, p(2)+0.17, p(3)-0.15, p(4)-0.21]);
hold(h);
plot(h, tt(inset), reshape(phi(2,2,inset),[],1), 'k--','LineWidth',2);
plot(h, tt(inset), reshape(app(2,2,inset),[],1), 'r-','LineWidth',2);
xlim([2,5])
grid on
set(gca,'FontSize',11);