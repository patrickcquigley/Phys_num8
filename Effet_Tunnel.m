%% Script for processing the Tunnel effect 

%% 
plot_case='b_left';
 
%cases 
% a_right - paticle on the right plots wave 
% b_right - particle on the right plots prob of tunneling 
% a_left  - particle on the left plots wave 
% b_left  - particle on the left plots probe of tunneling 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change x0 in the config file if you switch to the right/left cases 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% params
tfin=1000;
xl=-150;
xr=150;
w=1e-2;
w2=2e-3;
delta=50;
x0=-50;% change as needed
%n=24; changed between simulations 
sigma_norm=0.04;
Nintervals=200;
Nsteps=1000;

dt= tfin./Nsteps ;

%% Paramters for 
repertoire = './'; % './' on Linux, '' on Windows
executable = 'Exercice_2023_V5'; % Nom de l'executable


input = 'configuration.in.example';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change n according to the simulation
switch plot_case 
     
    case {'a_left','a_right'}

        n = [20,30,40] ;% n values for part a 
    case {'b_right','b_left'}  
        n       = [20,22,24,26,28,30,32,34,36,38,40]; % n values for part b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
nsimul   = numel(n);
paramstr = 'n'; % Nom du parametre a scanner, par exemple dt, w, x0, etc
param    = n; % Valeurs du parametre a scanner

%%
output = cell(1, nsimul);

for ii = 1:nsimul
    % Variant to scan N1 and N2 together:
    filename  = [paramstr, '=', num2str(param(1,ii))];
    output{ii} = [filename, '.out'];
    eval(sprintf('!%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(ii), output{ii}));
    disp('Done.')
end

%% Load the data 

t_trans=zeros(nsimul,1);
prob_trans=zeros(nsimul,1);
E_mean=zeros(nsimul,1);
V_max=zeros(nsimul,1);



for ii= 1:nsimul
    
    fichier = output{ii};
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
    psiRe  = squeeze(wave(:, 2, :));
    psiIm = squeeze(wave(:, 3, :));
    
    psiMag=sqrt(psiRe.^2+psiIm.^2);

    % Uncertainty in x, p; eqs. (9-10-11)
    dx = sqrt(x2moy - xmoy.^2);
    dp = sqrt(p2moy - pmoy.^2);
    uncertainty = dx .* dp;

    t_trans(ii)=t(301);
    
    switch plot_case
        
        case 'b_left'
    
            prob_trans(ii)=P2(301);
        
        case 'b_right'
        
            prob_trans(ii)=P1(301);
    end
    
    E_mean(ii)=mean(E);
    V_max=max(V(50:150));
    
    
    switch plot_case
        
        case {'a_right','a_left'}
    
        figure 
        hold on
        box on
        [X,T] = meshgrid(x,t);
        pcolor(X,T,psiMag)
        colormap('default')
        colorbar
        plot([16.5,16.5],[0,1000],'color','red','Linewidth',1,'linestyle','--')
        xlabel('x')
        ylabel('t')
        shading interp
        set(gca,'fontsize',15)
        legend('|\Psi|','x_c')

        figure
        plot(t,P1,t,P2,t,P1+P2,'linewidth',1)
        grid on
        box on
        xlabel('t')
        ylabel('Probability')
        legend('P_{x<x_c}','P_{x>x_c}','P_{tot}')
        set(gca,'fontsize',15)

        figure
        hold on;
        box on
        plot(x,V,'linewidth',1)
        plot(x, E(1) * ones(size(x)),'linewidth',1)
        grid on
        xlabel('x')
        ylabel('V(x)')
        legend('V(x)', 'E(t=0)');
        set(gca,'fontsize',15)
    end

end 

switch plot_case
    
    case {'b_left','b_right'}

    x_axis=E_mean./V_max;


    figure 
    grid on 
    hold on 
    box on
    plot([1,1],[0,1],'color','k','linestyle',':')
    plot(x_axis,prob_trans,'linewidth',2)
    xlabel('$\langle E \rangle$/$V_{max}$','interpreter','latex')
    ylabel('$P_{trans}$','interpreter','latex')
    set(gca,'fontsize',15)

end

