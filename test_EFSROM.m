

clear all
clc

x_l=0;
x_r=1;

t=[0,1];

N=10000;
%n=1024;
NGrids=2^10-1;
n=NGrids+1;

opt.nu=0.001;
opt.theta=0;
opt.del=1e-4;
podopt.p=6;            %-----number of POD basis
podopt.L2_space=1;      % 1 L2, 0 H1 basis
podopt.plot_basis=0;
%sigma=0.3;
sigma=0.05:0.05:0.6;
de_efrom=0.0011;                   %---=0 GROM, otherwise, LROM
save_index=100;
n_gp=5;                 %---quadrature points
degree=1;
load_mesh=1;

if load_mesh==1;
    %%------load mesh
    %load fem_mesh_data.mat
    %Mass=M1;
    %Stiff=S1;
    load fem_mesh_data_1023.mat
else
    %%---generate mesh
    a_name=1;
    mesh=mesh_generator_1D(x_l,x_r,n);
    GDOF=global_dof_1D_fe_Lagrange(mesh,degree);
    FEM=struct('mesh',mesh,'GDOF',GDOF,'degree',degree);
    [Mass,~,Stiff]=FE_matrix_1D_Lagrange(a_name,FEM,n_gp);       %-----FE Mass, Stiff matrix
    fprintf('FE Mesh info complete !\n')
end


run_time=1000;
n_sig=length(sigma);
KE_efROM=zeros(n_sig,run_time);
KE_GROM=zeros(n_sig,run_time);
KE_dns=zeros(n_sig,run_time);
Err_GROM=zeros(n_sig,run_time);
Err_efROM=zeros(n_sig,run_time);

for j=1:n_sig
    sig=sigma(j);
    
    err_grom=[];
    err_efrom=[];
    err_lrom=[];
    
    ke_grom=[];
    ke_efrom=[];
    ke_dns=[];
    fprintf('compute the new sigma\n')
    for i=1:run_time
        
       
        [dW,u_snap,u_fem]=run_stoch_Burgers(sig,save_index);
        
        %----sample data
        %u_snap=u_dns(:,1:save_index:end);
        Eng_avg_dns = sum(diag(u_snap'*Mass*u_snap))/size(u_snap,2);
        %-----compute POD basis, matrix, tensor
        if i==1
            fprintf('generate POD\n')
            fprintf('     \n')
            [pod_u,CumEng,CumEng_ratio,POD_all,Diag_S,lp]=POD_basis(FEM,u_snap,Mass,Stiff,podopt);             %-------generate POD basis;
            r=podopt.p;
            [~,Ten_0x,~]=POD_tensor_assemble_1D(r,pod_u,FEM,n_gp);        %---get the tensor;
            Mr=pod_u'*Mass*pod_u;
            Sr=pod_u'*Stiff*pod_u;
            save Mr Sr Ten_0x pod_u
        else
            fprintf('use same POD info, load...\n')
            load Mr Sr Ten_0x pod_u
        end
        
        u0=u_fem(:,1);
        C0=Mr\pod_u'*Mass*u0;
        
        [Err_L2_avg_grom,Err_L2_avg_efrom,C_grom,C_efrom]=Stochastic_EFROMs(N,podopt,...
    opt,de_efrom,Mass,Stiff,Eng_avg_dns,pod_u,u_snap,dW,sig,Mr,Sr,Ten_0x,C0,save_index);

        %[Err_L2_avg_grom,Err_L2_avg_efrom,Eng_avg_dns,Eng_avg_grom,Eng_avg_efrom,C_grom,C_efrom]=Stochastic_EFROMs(N,podopt,opt,...
        %    de_efrom,Mass,Stiff,pod_u,u_fem,dW,sig,Mr,Sr,Ten_0x,C0,save_index);
        %[Err_L2_avg_grom,Err_L2_avg_efrom,Err_L2_avg_lrom,Eng_avg_dns,Eng_avg_grom,Eng_avg_lrom,C_grom,C_efrom,C_lrom]=Stochastic_EFROMs(N,podopt,...
    %opt,de_efrom,de_lrom,Mass,Stiff,pod_u,u_fem,dW,sig,Mr,Sr,Ten_0x,C0,save_index);
        
        %ke_grom=[ke_grom,Eng_avg_grom];
        %ke_efrom=[ke_efrom,Eng_avg_efrom];
        %ke_dns=[ke_dns,Eng_avg_dns];
        err_grom=[err_grom,Err_L2_avg_grom];
        err_efrom=[err_efrom,Err_L2_avg_efrom];
        
        fprintf('sigma=%2f, runtime=%d\n',sig,i)
    end
    %fprintf('Complete one sigma, delete the POD..\n')
    %clear Mr Sr Ten_0x pod_u
    
    %KE_efROM(j,:)=ke_efrom;
    %KE_GROM(j,:)=ke_grom;
    %KE_dns(j,:)=ke_dns;
    Err_GROM(j,:)=err_grom;
    Err_efROM(j,:)=err_efrom;  
    j,
end

fprintf('Complete !\n')
u_grom=pod_u*C_grom;
Zmax= max(max(u_grom))*1;
Zmin= min(min(u_grom))*1;
t_initial=0;
t_final=1;

% GDOF=FEM.GDOF;
% dt=0.01;
% t_checkpts=0:dt:1;
% x=GDOF.P_g;
% figure(1)
% mesh(x,t_checkpts', u_grom');
% colormap jet
% %xlabel('x'); ylabel('t');
% %view([-1 -1 1]);
% axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
% %title('DNS')
% %title(['GROM r=', num2str(Tr)],'Fontsize',10);
% %title(['EFROM r=', num2str(Tr)],'Fontsize',10);
% %set(gca,'fontsize',10,'fontweight','b')
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', 16);
% ylabel('$t$', 'interpreter', 'latex', 'fontsize', 16);
% zlabel('$u$', 'interpreter', 'latex', 'fontsize', 16);
% if sig == 0
%     zlim([0 1.5])
% end

% figure(2)
% u_efrom=pod_u*C_efrom;
% Zmax= max(max(u_efrom))*1;
% Zmin= min(min(u_efrom))*1;
% mesh(x,t_checkpts', u_efrom');
% colormap jet
% %xlabel('x'); ylabel('t');
% %view([-1 -1 1]);
% axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', 16);
% ylabel('$t$', 'interpreter', 'latex', 'fontsize', 16);
% zlabel('$u$', 'interpreter', 'latex', 'fontsize', 16);

% figure(3)
% mesh(x,t_checkpts', u_snap');
% colormap jet
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', 16);
% ylabel('$t$', 'interpreter', 'latex', 'fontsize', 16);
% zlabel('$u$', 'interpreter', 'latex', 'fontsize', 16);

% figure(4)
% hold on
% plot(Err_efROM,'r')
% plot(Err_GROM,'b')
% hold off
% legend('EF-ROM','G-ROM')
% title('Error')
% figure(5)
% hold on
% plot(KE_dns,'b')
% plot(KE_efROM,'r')
% plot(KE_GROM,'g')
% hold off
% title('KE evolution')
% tt = (0:1:(N-1))'*dt; % time vecto
% %uu= [zeros(1,N);u;zeros(1,N)]; % add the boundary value back
% [XX,TT] = meshgrid([0;dis_x;L],tt); 
% mesh(XX,TT,uu);
% view([-19 38])
% colormap jet

% figure(1)
% plot(KE_dns,'ro')
% hold on
% plot(KE_dns_1time,'b*')
% hold off
% title('Average KE for DNS')
% legend('update','no-update')
% xlabel('runtime')
% ylabel('KE')
% figure(2)
% plot(KE_rom,'ro')
% hold on
% plot(KE_rom_1time,'b*')
% hold off
% title('Average KE for ROM')
% legend('update','no-update')
% xlabel('runtime')
% ylabel('KE')
% figure(3)
% plot(Err,'ro')
% hold on
% plot(Err_1time,'b*')
% hold off
% title('Average L2 error for ROM')
% legend('update','no-update')
% xlabel('runtime')
% ylabel('Error')
% 
%   plot(ke_lrom,'bo')
%   hold on
%   plot(KE_GROM,'r*')
%   plot(KE_dns,'ks')
%   hold off
%   plot(Err_LROM,'bo')
%   hold on
%   plot(Err_GROM,'r*')
%   hold off