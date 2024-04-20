%-----------------------------------------------------
%  Project "Modeling and Control of CPS"
%           TASK 4 - Sparse observer
%                               
%                                Carlo Migliaccio
%                                12.04.2024
%-----------------------------------------------------
clear all
close all
clc 

load tracking_moving_targets.mat    %load A,D,Y

p=100;          %number of the cells
q=25;           %number of sensors
Ntarget=3;
Nattack=2;

lammba1=10; lambda2=20; 
lambda=[10*ones(p,1); 20*ones(q,1)];
G = [D eye(q)];                 %augmented sensing matrix
G = normalize(G); 
eps=1e-8; 
tau = (norm(G)^(-2))-eps; 

%------------------------SPARSE OBSERVER------------------------
Tmax=50;

%z_hat = [xtrue; zeros(q,1)];       %stato iniziale 

z_hat = [zeros(p,1); zeros(q,1)];       %stato iniziale 

mes_x = zeros(p,1); 
mes_a = zeros(q,1); 

z_hat_plus=zeros(p+q,1);
for k=0:(Tmax-1)
    %to do...
    arg=z_hat + tau*G'*(Y(:,k+1)-G*z_hat);                  %arg of STO
    %to do... Estimation (apply STO for each arg(i))
    for i=1:(p+q)
        z_hat_plus(i) = sto(arg(i), tau*lambda(i));        
    end
    %to do... Prediction
    z_hat=[A*z_hat_plus(1:p); z_hat_plus(p+1:end)];
    
    x_hat=A*z_hat_plus(1:p);
    a_hat=z_hat_plus(p+1:end);

    mes_x = [mes_x x_hat];
    mes_a = [mes_a a_hat];
end

%data cleaning
for j=0:(Tmax)
    %cleaning x_hat
    max_x_vec = maxk(mes_x(:,j+1),Ntarget);
    for i=1:p
        if(abs(mes_x(i,j+1))<max_x_vec(Ntarget))
            mes_x(i,j+1)=0; 
        end
    end

    %cleaning a_hat
    max_a_vec = maxk(mes_a(:,j+1),Nattack);
    for i=1:q
        if(abs(mes_a(i,j+1))<max_a_vec(Nattack))
            mes_a(i,j+1)=0; 
        end
    end
end

room(mes_x,mes_a,1,Tmax);



