%----------Differential fitler process
%----
% format long
%% -------May 17, 2016
%---latest version
clc;
clear all;

time_star = tic;

x_l=0;
x_r=1;

t=[0,1];

N=10000;
n=1025;

opt.nu=0.001;
opt.theta=0;
opt.del=1e-4;
%sig=0.5;
%Tr=2:4:50;      % # of pod basis
Tr=10;
podopt.L2_space=1;  % 1 L2, 0 H1 basis
podopt.plot_basis=1;
EL2=[];
EH1=[];
EHs=[];
save_index=100;

degree=1;
n_gp=5;


podopt.p=Tr;

load fem_mesh_data.mat
Mass=M1;
Stiff=S1;
%load u_fem_data_sig05_0.mat

run_time=1;
sig=0.3;
ab_sig=3;
for ll=0:run_time-1 
    if sig==0
        fileID=['u_fem_data_sig0.mat'];
        fprintf('This is determinstic Burgers\n')
        load(fileID);
    else
        fileID=['u_fem_data_sig0',num2str(ab_sig),'_',num2str(ll),'.mat'];
        fprintf('load data form file\n')
        fprintf('%s\n',fileID);
        load(fileID,'dW','u_fem');
        fprintf('simulate sig=%f, run_time=%d\n',sig,ll);
    end
u0=u_fem(:,1);


[pod_u,CumEng,CumEng_ratio,POD_all,Diag_S,lp]=POD_basis(FEM,u_fem,Mass,Stiff,podopt);             %-------generate POD basis;   
r=podopt.p;

[Ten,~,~]=POD_tensor_assemble_1D(r,pod_u,FEM,n_gp);        %---get the tensor;

Nonlin0 = zeros(r,r^2);
Nonlin1 = zeros(r,r^2);
Ones = zeros(r^2,1);
C=sparse(r,r);

for i=1:r
    Ones((i-1)*r+1:i*r,1:r)=eye(r);
end

Mr=pod_u'*Mass*pod_u;
Sr=pod_u'*Stiff*pod_u;

C0=Mr\pod_u'*Mass*u0;

C_temp=C0;

 
C(:,1)=C_temp;
%C0=C(:,1);


%de = 0.15;
de=0;

p=0;
E=[];
%Wt=dW./sqrt(opt.del);
%Wt=dW;
%for j=1:length(de)
 
    C0=C(:,1);
    for k=1:N
        
        C_bar=(Mr+(de^2)*Sr)\(Mr*C0);
        %C_bar=(C0);
         Nh=zeros(r,1);
%         for i=1:r
%             Nonlin0(i,(i-1)*r+1:i*r)=(Ten(:,:,i)*C0)';
%         end
%         
%         Nh=Nh+Nonlin0*Ones*C_bar;
        
            for i=1:r
                Nh(i) = (C_bar)'*Ten(:,:,i)*C0; %
            end

        %Nh=C_bar'*Nonlin
        F=Mr*C0-opt.nu*opt.del*(Sr*C0)-opt.del*Nh+sig*dW(k)*Mr*C0+(sig^2/2)*opt.del*Mr*C0;       %---F vector
        
        C1=Mr\F;
        if rem(k,save_index)==0
            C(:,k/save_index+1)=C1;
        end
        C0=C1;
        k=k+1;
    end
    
    u_rom=pod_u*C;
    
    %figure(2)
    %plot(GDOF.P_g,u_rom(:,end),'ob');
    
   %Av=POD_AverDNS_new(u_rom,Mass);
   %E=[E,Av];
    
    %-----Caculate the error btw  FEM and POD ROM
%     
     Err_L2=zeros(101,1);
     Err_H1=Err_L2;
     Err_H1_semi=Err_L2;
% % %     Err_max=Err_L2;
     
    for i=1:101
        
        [Err_max(i),Err_L2(i),Err_H1_semi(i),Err_H1(i)]=...
            POD_error_1D(Mass,Stiff,u_fem(:,i),u_rom(:,i));
    end
     
     Err_L2_avg=sqrt(sum(Err_L2.^2)/101)
     Err_H1_avg=sqrt(sum(Err_H1.^2)/101);
     Err_H1_semi_avg=sqrt(sum(Err_H1_semi.^2)/101);     
     
     EL2=[EL2,Err_L2_avg];
     EH1=[EH1,Err_H1_avg];
     EHs=[EHs,Err_H1_semi_avg];
%    p=p+1,
    ll,
end

ddd=POD_AverDNS_new(u_fem,Mass);

cputime = toc(time_star);

Zmax= max(max(u_rom))*1;
Zmin= min(min(u_rom))*1;
t_initial=0;
t_final=1;

GDOF=FEM.GDOF
dt=100/N;
t_checkpts=0:dt:1;
x=GDOF.P_g;
figure(3)
mesh(x,t_checkpts', u_fem');
colormap jet
xlabel('x'); ylabel('t');
%view([-1 -1 1]);
axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
title('DNS')
%title(['GROM r=', num2str(Tr)],'Fontsize',10);
%title(['LROM r=', num2str(Tr)],'Fontsize',10);
set(gca,'fontsize',10,'fontweight','b')
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 15);
ylabel('$t$', 'interpreter', 'latex', 'fontsize', 15);
zlabel('$u$', 'interpreter', 'latex', 'fontsize', 15);
if sig == 0
    zlim([0 1.5])
end

figure(4)
 h=pcolor(x,t_checkpts',u_rom');
 set(h,'edgecolor','none')
 colormap jet
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 15);
ylabel('$t$', 'interpreter', 'latex', 'fontsize', 15);

figure(5)
x=zeros(run_time,1);
x(:,1)=0.3;
plot(x,EL2,'ro')
hold on
title('Error')

%colorbar
% r=zeros(length(E));
% r=Av;
% figure(2)
% subplot(2,1,1)
% plot(de,E,'bo--','Linewidth',2.5)
% hold on
% plot([0 max(de)],[ddd ddd],'r-','Linewidth',2.5)
% hold off
% legend('Leray-DF','DNS');
% title(['Average energy vs delta at r= ', num2str(Tr)])
% xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel('Average energy','Fontsize',15);
% subplot(2,1,2)
% plot(de,EL2,'bo--','Linewidth',2.5)
% title(['Average errors vs delta at r= ', num2str(Tr)])
% 
% %%-----
% figure(1)
% plot(dr,E,'b','Linewidth',2.5)
% hold on
% plot([0 max(dr)],[ddd ddd],'r','Linewidth',2.5)
% hold off
% legend('EF-Proj','DNS');
% title('Average L2-norm vs delta with r=6','Fontsize',15);
% xlabel('r1','Fontsize',15);
% %xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel(['$\frac{1}{101}\sum_{j=1}^{101}\|u_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% 
% % figure(2)
% % plot(dr,EL2,'b','Linewidth',1.5)
% % %legend('L-ROM-DF','Fontsize',5);
% % hgd = legend('EF-ROM-Proj');
% % set(hgd,'Position',([0.18,0.85,0.27,0.06]))
% % %title(['Average errors vs delta at r= 6'],'Fontsize',15);
% % xlabel('r1','Fontsize',12);
% % %xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% % ylabel(['$\frac{1}{101}\sum_{j=1}^{101}\|e_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% % %ylabel('Average Error','Fontsize',15);
% % 
% % 
% % 
% % 
% % 
% % figure(1)
% % %r=1:1:101;
% % r=0:0.01:1;
% % rgb_value=[0,127,38];
% % rgb_value=rgb_value/255;
% % plot(r,Err_L2,'b','Linewidth',2)
% % hold on
% % plot(r,Err_L2grom,'g','Linewidth',2,'color',rgb_value)
% % hold off
% % hgd = legend('EF-ROM-Proj','G-ROM');
% % %hgd = legend('G-ROM');
% % set(hgd,'Position',[0.15,0.79,0.27,0.12],'fontsize',11);
% % xlabel(['$t','$'],'Interpreter','latex','Fontsize',15);
% % %xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% % ylabel(['$\|e_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',18);
% % %ylabel('Error','Fontsize',15);
% % 
% % 
% % 
% % 
% % 
% % figure(1)
% % r=0:0.01:1;
% % rgb_value=[0,127,38];
% % rgb_value=rgb_value/255;
% % %plot(r,Err_L2,'b','Linewidth',1.5)
% % hold on
% % plot(r,Err_L2,'g','Linewidth',1.8,'color',rgb_value)
% % %plot(r,Err_L2grom,'g','Linewidth',1.5)
% % hold off
% % %hgd = legend('EF-ROM-Proj','G-ROM');
% % hgd = legend('G-ROM');
% % set(hgd,'Position',[0.15,0.85,0.27,0.06],'fontsize',11);
% % set(gca,'box','on')
% % %title('Time evolution of error of L-ROM with r=6','Fontsize',15);
% % xlabel(['$t','$'],'Interpreter','latex','Fontsize',15);
% % %xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% % ylabel(['$\|e_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',18);
% % 
% % 

% figure(3)
% mesh(x,t_checkpts', u_fem');
% colormap jet
% xlabel('x'); ylabel('t');
% view([-1 -1 1]);
% axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
% title(['GROM r=', num2str(Tr)],'Fontsize',10);
% title(['DNS'],'Fontsize',10);
% set(gca,'fontsize',10,'fontweight','b')
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', 15);
% ylabel('$t$', 'interpreter', 'latex', 'fontsize', 15);
% zlabel('$u$', 'interpreter', 'latex', 'fontsize', 15);
% if sig == 0
%     zlim([0 1.5])
% end
% 
% figure(4)
%  h=pcolor(x,t_checkpts',u_fem');
%  set(h,'edgecolor','none')
%  colormap jet
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', 15);
% ylabel('$t$', 'interpreter', 'latex', 'fontsize', 15);



