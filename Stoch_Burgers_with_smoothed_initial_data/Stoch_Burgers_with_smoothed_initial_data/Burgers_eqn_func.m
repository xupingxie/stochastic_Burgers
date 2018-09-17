%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code integrates the Burgers-type equation 
% du =(nu*u_xx + lam*u - gam*u*u_x)*dt + sig*u*dW.
% See Sect. 10.1 of the Reference given below for more details about 
% the description of the numerical method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference: 
% M. D. Chekroun, H. Liu, and S. Wang, On stochastic parameterizing manifolds: Pullback
% characterization and Non-Markovian reduced equations, Preprint, 143 pp., 2014. 
% arXiv link: http://arxiv.org/abs/1310.3896.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function uu = Burgers_eqn_func(lam, gam, nu, sig, L, u0, dt, N, NGrids, dW)



dx = L/(NGrids+1); % spatial grid size

eignLap = 2/dx^2*(cos(pi*(1:1:NGrids)'/(NGrids+1)) - 1); % eigenvalues for discrete Laplacian
SeignLap = nu*dt*eignLap; % scaled eigenvalues for discrete Laplacian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop
uu = zeros(NGrids, N); % to store the SPDE solution
uu(:,1) = u0;
u_cur = u0;

fprintf('Solving the SPDE ...\n');  
for t = 2:N
    % print to screen
    if mod(t, 1e4) == 0
        fprintf(['     processing i = %d ... ','\n'], t);        
    end
    
        
    % dealing with the nonlinear term in physical space
    u2 = u_cur.^2;
    u2_extend = [0;u2;0];
    u2_x = u2_extend;
    for j = 2:NGrids+1
        u2_x(j) = (u2_extend(j) - u2_extend(j-1))/(dx);
    end
    u2_x = - 0.5*gam*u2_x(2:(NGrids+1))*dt;
        
    % transform the RHS to frequence space
    %%%RHS = dst(u_cur + u2_x + sig*dW(t-1)*u_cur + 0.5*sig^2*u_cur*dt); 
    tempy=u_cur + u2_x + sig*dW(t-1)*u_cur + 0.5*sig^2*u_cur*dt;
    RHS = dstn(tempy, 1); 
    % note that the term 0.5*sig^2*u_cur in the above line is due to the
    % conversion between the Ito and Stratonovich stochastic integral.
    
    % compute the coef. for LHS    
    LHS = (1-SeignLap - lam*dt);    
    
    % u_next in freq space   
    u_hat = RHS./LHS;
    
    % u_next in physical space
    %%%u_next = idst(u_hat);
    u_next = idstn(u_hat,1);
    
    u_cur = u_next;
  
    % save data
    uu(:,t) = u_next;
    
    % check if solution blows up
    if isnan(uu(2,t)) == 1
        error('SPDE solution blows up')
    end
end
fprintf('End of solving the SPDE.\n\n'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




