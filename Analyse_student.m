eval(sprintf('!make student'))
eval(sprintf(['!patrick.exe']))

%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier = 'output';
data  = load([fichier,'_obs.out']);
t     = data(:,1);
P1    = data(:,2);
P2    = data(:,3);
E     = data(:,4);
xmoy  = data(:,5);
x2moy = data(:,6);
pmoy  = data(:,7);
p2moy = data(:,8);
data  = load([fichier,'_pot.out']);
x     = data(:,1);
V     = data(:,2);

wave  = reshape(load([fichier,'_psi2.out']), length(t), 3, length(x));
psi2  = squeeze(wave(:, 1, :));
real_psi=squeeze(wave(:, 2, :));



% Uncertainty in x, p; eqs. (9-10-11)
dx = sqrt(x2moy - xmoy.^2);
dp = sqrt(p2moy - pmoy.^2);
uncertainty = dx .* dp;

%% Figures %%
%%%%%%%%%%%%%
% figure('Name',['Analyse de ' fichier],'Position',[300,500,1200,1800])
% subplot(3,2,1)
% hold on;
% plot(x,V)
% plot(x, E(1) * ones(size(x)));
% grid on
% xlabel('x')
% ylabel('V(x)')
% legend('V(x)', 'E(t=0)');
% 
% subplot(3,2,2)
% plot(t,P1,t,P2,t,P1+P2)
% grid on
% xlabel('t')
% legend('P_{x<x_c}','P_{x>=x_c}','P_{tot}')
% 
% subplot(3,2,3)
% [X,T] = meshgrid(x,t);
% pcolor(X,T,psi2)
% xlabel('x')
% ylabel('t')
% shading interp
% 
% subplot(3,2,4)
% plot(t,xmoy,t,x2moy,t,pmoy,t,p2moy);
% xlabel('t')
% legend('\langlex\rangle','\langlex^2\rangle','\langlep\rangle','\langlep^2\rangle')
% grid on
% 
% subplot(3,2,5)
% plot(t, uncertainty);
% hold on
% plot(t,0.5*t./t);
% xlabel('t')
% legend('\langle \Delta x \rangle \cdot \langle \Delta p \rangle','h_{bar}/2')
% grid on
% 
% subplot(3,2,6)
% plot(t, E);
% xlabel('t[s]')
% ylabel('E [J]')
% grid on

% figure
% plot(t,xmoy)
% set(gca,"FontSize",14)
% xlabel('t [s]')
% ylabel(['\langle x \rangle [m]'])
% grid on

% figure
% plot(t,pmoy)
% set(gca,"FontSize",14)
% xlabel('t [s]')
% ylabel('\langle p \rangle [kg m s^{-1}]')
% grid on
% 
figure
plot(t,dx)
set(gca,"FontSize",14)
xlabel('t [s]')
ylabel('\Delta x [m]')
grid on

figure
plot(t,dp)
set(gca,"FontSize",14)
xlabel('t [s]')
ylabel('\Delta p [m]')
grid on


% figure
% [X,T] = meshgrid(x,t);
% pcolor(X,T,real_psi)
% set(gca,"FontSize",14)
% xlabel('x')
% ylabel('t')
% shading interp

% [X,T] = meshgrid(x,t);
% pcolor(X,T,psi2)
% set(gca,"FontSize",14)
% xlabel('x [x]')
% ylabel('t [s]')
% shading interp

% [X,T] = meshgrid(x,t);
% pcolor(X,T,real_psi)
% set(gca,"FontSize",14)
% xlabel('x [x]')
% ylabel('t [s]')
% shading interp

%%
%part 2
%simplify because m=1
w1= 0.01;
m=1;

x_classic=sqrt(2.*E./m)*1./w1.*sin(w1.*t);
p_classic=m*sqrt(2.*E./m).*cos(w1.*t);
% 
% figure
% plot(t,x_classic)
% set(gca,"FontSize",14)
% xlabel('t [s]')
% ylabel('x_{classic} [m]')
% grid on

% figure
% plot(t,p_classic)
% set(gca,"FontSize",14)
% xlabel('t [s]')
% ylabel('x_{classic} [m]')
% grid on
% % 
% figure
% plot(t,xmoy)
% set(gca,"FontSize",14)
% xlabel('t [s]')
% ylabel('\Delta x [m]')
% grid on
% 
% figure
% plot(t,pmoy)
% set(gca,"FontSize",14)
% xlabel('t [s]')
% ylabel('\Delta p [m]')
% grid on