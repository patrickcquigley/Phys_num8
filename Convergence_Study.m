%% Convergence Study 
% script to perform the convergence study for Exercise 8 - by varying nstep

% params 

tfin=1000;
xl=-150;
xr=150;
w=1e-2;
w2=2e-3;
delta=0;
x0=0;
n=24;
sigma_norm=0.04;

Nintervals=200;

%% Variables fo simulation + config file 
repertoire = './'; % './' on Linux, '' on Windows
executable = 'Exercice_2023_V5'; % Nom de l'executable


input = 'configuration.in.example';
Nsteps       = [2000,4000,8000,16000,32000,64000,128000];

nsimul   = numel(Nsteps);
paramstr = 'Nsteps'; % Nom du parametre a scanner, par exemple dt, w, x0, etc
param    = Nsteps; % Valeurs du parametre a scanner

dt= tfin./Nsteps ;
%% Run simulations varyiong the number of steps 
output = cell(1, nsimul);

for ii = 1:nsimul
    % Variant to scan N1 and N2 together:
    filename  = [paramstr, '=', num2str(param(1,ii))];
    output{ii} = [filename, '.out'];
    eval(sprintf('!%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(ii), output{ii}));
    disp('Done.')
end

%% Load data for plotting - Extract momentum and position at the final time 

x_tfin=zeros(nsimul,1);
p_tfin=zeros(nsimul,1);
deltax_tfin=zeros(nsimul,1);
deltap_tfin=zeros(nsimul,1);


for ii= 1:nsimul
    
    fichier = output{ii};
    data  = load([fichier,'_obs.out']);
   
    xmoy  = data(:,5);
    x2moy = data(:,6);
    pmoy  = data(:,7);
    p2moy = data(:,8);
    data  = load([fichier,'_pot.out']);


    dx = sqrt(x2moy - xmoy.^2);
    dp = sqrt(p2moy - pmoy.^2);

    x_tfin(ii)=xmoy(end);
    p_tfin(ii)=pmoy(end);
    deltax_tfin(ii)=dx(end);
    deltap_tfin(ii)=dp(end);
    
end 

%% Plotting 
% fit for x moy 
dt_cont=linspace(0,max(dt),100);

p_x=polyfit(dt.^2,x_tfin,1);
plot_x=p_x(1)*dt_cont.^2+p_x(2);

label=sprintf(' y= %.3f x + %.3f', p_x(1), p_x(2));

figure 
hold on 
box on 
grid on
axis padded 
plot(dt.*dt,x_tfin,'x','markersize',8,'color','#0072BD','LineWidth', 2)
plot(dt_cont.*dt_cont,plot_x,'linestyle','--','color','#0072BD','linewidth',1)
xlabel('$\Delta t^2$','interpreter','latex');
ylabel('$\langle x \rangle (t_{fin}) $','interpreter','latex');
set(gca,'fontsize',15)
legend('$\langle x \rangle (t_{fin}) $',label,'location','SE','interpreter','latex');
%%

p_p=polyfit(dt.^2,p_tfin,1);
plot_p=p_p(1)*dt_cont.^2+p_p(2);

label=sprintf(' y= %.3f x + %.3f', p_p(1), p_p(2));

figure 
hold on
box on 
grid on
axis padded 
plot(dt.*dt,p_tfin,'x','markersize',8,'color','#D95319','LineWidth', 2)
plot(dt_cont.*dt_cont,plot_p,'linestyle','--','color','#D95319','linewidth',1)
xlabel('$\Delta t^2$','interpreter','latex');
ylabel('$\langle p \rangle (t_{fin}) $ ','interpreter','latex');
set(gca,'fontsize',15)
legend('$\langle p \rangle (t_{fin}) $',label,'location','NE','interpreter','latex');

%%
p_dx=polyfit(dt.^2,deltax_tfin,1);
plot_dx=p_dx(1)*dt_cont.^2+p_dx(2);

label=sprintf(' y= %.3f x + %.3f', p_dx(1), p_dx(2));

figure
hold on
box on 
grid on
axis padded 
xlabel('$\Delta t^2$','interpreter','latex');
plot(dt.*dt,deltax_tfin,'x','markersize',8,'color','#77AC30','LineWidth', 2)
plot(dt_cont.*dt_cont,plot_dx,'linestyle','--','color','#77AC30','linewidth',1)
ylabel('$\langle \Delta x \rangle (t_{fin}) $ ','interpreter','latex');
set(gca,'fontsize',15)
legend('$\langle \Delta x \rangle (t_{fin})$',label,'location','SE','interpreter','latex')

%%
p_dp=polyfit(dt.^2,deltap_tfin,1);
plot_dp=p_dp(1)*dt_cont.^2+p_dp(2);

label=sprintf(' y= %.3f x + %.3f', p_dp(1), p_dp(2));


figure 
hold on
box on 
grid on
axis padded 
plot(dt.*dt,deltap_tfin,'x','markersize',8,'color','#7E2F8E','LineWidth', 2)
plot(dt_cont.*dt_cont,plot_dp,'linestyle','--','color','#7E2F8E','linewidth',1)
xlabel('$\Delta t^2$','interpreter','latex');
ylabel('$\langle \Delta p \rangle (t_{fin}) $','interpreter','latex');
set(gca,'fontsize',15)
legend('$\langle \Delta p \rangle (t_{fin}) $',label,'location','SE','interpreter','latex')