%=========================================================================
% This code simulates the Burgers-type SPDE driven by linear multiplicative noise 
% See Eq.(10.12) in Section 10.2 of the reference
% given below for more details.
%=========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference: 
% M. D. Chekroun, H. Liu, and S. Wang, On stochastic parameterizing manifolds: Pullback
% characterization and Non-Markovian reduced equations, Preprint, 143 pp., 2014. 
% arXiv link: http://arxiv.org/abs/1310.3896.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dW,u_snap,u_fem]=run_stoch_Burgers(sig,save_index)

%clc; 
%clear; 
%close all;
%addpath('./code_mollifier');


initial_flag = 1; % =0, step func; = 1, mollified version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Parameters 
L = 1;
NGrids = 2^10-1; % grid points for SPDE
nu=10^-3;
%sig = 0;
%sig = 0.1;
gam = 1;
lam = 0;

%delta = 0.2;  % percent of lam away from its critical value.
%lam = nu*pi^2/L^2*(1+delta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1e-4;  % dt = time step size for SPDE            
N =  ceil(1/dt);   % N = total number of iterations for SPDE

dx = L/(NGrids+1); % spatial grid size
dis_x = (linspace(dx, L-dx, NGrids))';  

dx2 = dx^2;
eignLap = 2/dx^2*(cos(pi*(1:1:NGrids)'/(NGrids+1)) - 1); % eigenvalues for discrete Laplacian
SeignLap = nu*dt*eignLap; % scaled eigenvalues for discrete Laplacian


% initial data for SPDE
%u0 = 0.1*eigens(:,1) + 0.1*eigens(:,2);

if initial_flag == 0
    u0 = zeros(NGrids,1);
    u0(1:floor(NGrids/2)) = 1;
elseif initial_flag == 1
    romberg_par = 10;
    myfunc = @step_func;
    x_grids = dis_x;
    epsilon = 0.01;
    y_m = mollifier(myfunc, x_grids, epsilon,romberg_par);
    y = step_func(x_grids);
    u0 = y_m;
end


% figure;
% plot(x_grids,y,'.-k');
% hold on; grid on;
% plot(x_grids,y_m,'.-r');
% str = legend('step func','mollified version');
% set(str, 'fontsize',17);

% Generate the noise
%seed = 100; randn('state',seed); % using a fixed noise path 

dW = sqrt(dt)*randn(1,N);             

% Simulate the SPDE solution
u = Burgers_eqn_func(lam, gam, nu, sig, L, u0, dt, N, NGrids, dW);


% plot the solution field
%figure; 
%tt = (0:1:(N-1))'*dt; % time vecto
uu= [zeros(1,N);u;zeros(1,N)]; % add the boundary value back
%[XX,TT] = meshgrid([0;dis_x;L],tt); 
%mesh(XX,TT,uu);
%view([-19 38])
%colormap jet
%set(gca,'fontsize',15,'fontweight','b')
%xlabel('$x$', 'interpreter', 'latex', 'fontsize', 20);
%ylabel('$t$', 'interpreter', 'latex', 'fontsize', 20);
%zlabel('$u$', 'interpreter', 'latex', 'fontsize', 20);
%if sig == 0
%    zlim([0 1.5])
%end

 
 %u_fem=uu;
% u_snap=u_fem(:,1:save_index:end);
 %save_index=100;
 u_fem=uu;
 %  for k=1:N
 %      if rem(k,save_index)==0
 %            u_snap(:,k/save_index+1)=uu(:,k);
 %      end          
 %  end
 %u_snap=u_fem(1:100:end);
 u_snap=zeros(NGrids+2, 101);
 u_snap(:,1:100)=u_fem(:,1:100:end);
 u_snap(:,end)=u_fem(:,end);
 end

