%----------Differential fitler process
%----
% format long
%%--------May 26, 2016, compute GROM, LROM at the same time

function [Err_L2_avg_grom,Err_L2_avg_lrom,Eng_avg_dns,Eng_avg_grom,Eng_avg_lrom,C_grom,C_lrom]=Stochastic_2ROMs(N,podopt,...
    opt,de_lrom,Mass,Stiff,pod_u,u_fem,dW,sig,Mr,Sr,Ten,rom_u0,save_index)
%time_star = tic;
%u0=u_fem(:,1);
%C0=Mr\pod_u'*Mass*u0;
%C_temp=C0;

%---C0 is the initial condition for the ROM time integration

r=podopt.p;

C_grom(:,1)=rom_u0;
C_lrom(:,1)=rom_u0;

C0_grom=C_grom(:,1);
C0_lrom=C_lrom(:,1);

for k=1:N-1
    
    
    C_bar=(Mr+(de_lrom^2)*Sr)\(Mr*C0_lrom);
    Nh_lrom=zeros(r,1); 
    Nh=zeros(r,1);
    for i=1:r
        Nh_lrom(i) = (C_bar)'*Ten(:,:,i)'*C0_lrom;        %----This is for LROM nolinear
        Nh(i)         = (C0_grom)'*Ten(:,:,i)'*C0_grom;            %----This is for gROM nolinear
    end    
    F_lrom=Mr*C0_lrom-opt.nu*opt.del*(Sr*C0_lrom)-opt.del*Nh_lrom+sig*dW(k)*Mr*C0_lrom+(sig^2/2)*opt.del*Mr*C0_lrom;       %---F vector    
    F_grom=Mr*C0_grom-opt.nu*opt.del*(Sr*C0_grom)-opt.del*Nh+sig*dW(k)*Mr*C0_grom+(sig^2/2)*opt.del*Mr*C0_grom;       %---F vector
    C1_lrom=Mr\F_lrom;
    C1_grom=Mr\F_grom;
    
    %if rem(k,save_index)==0
    %   C_grom(:,k/save_index+1)=C1_grom;
    %   C_lrom(:,k/save_index+1)=C1_lrom;
    %end
    C_grom(:,k+1)=C1_grom;
    C_lrom(:,k+1)=C1_lrom;
    
    C0_lrom=C1_lrom;
    C0_grom=C1_grom;
    k=k+1;
end

u_grom=pod_u*C_grom;
u_lrom=pod_u*C_lrom;


%---postprocessing data

%----compute the average KE for DNS and ROM

Eng_avg_grom=POD_AverDNS_new(u_grom,Mass);
Eng_avg_lrom=POD_AverDNS_new(u_lrom,Mass);
Eng_avg_dns=POD_AverDNS_new(u_fem,Mass);


%-------Caculate the error btw  FEM and POD ROM
time_step=size(u_grom,2),

Err_L2_grom=zeros(time_step,1);
Err_L2_lrom=zeros(time_step,1);
for i=1:time_step
    [~,Err_L2_grom(i),~,~]=...
        POD_error_1D(Mass,Stiff,u_fem(:,i),u_grom(:,i));
    [~,Err_L2_lrom(i),~,~]=...
        POD_error_1D(Mass,Stiff,u_fem(:,i),u_lrom(:,i));
end

Err_L2_avg_grom=(sum(Err_L2_grom)/time_step)./Eng_avg_dns;
Err_L2_avg_lrom=(sum(Err_L2_lrom)/time_step)./Eng_avg_dns;
%Err_L2_avg=sqrt(sum(Err_L2.^2)/time_step),


end




