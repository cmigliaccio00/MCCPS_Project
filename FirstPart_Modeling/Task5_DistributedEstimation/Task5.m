clear all
close all
clc

load distributed_localization_data.mat

n=100; 
q=25;

Q = Q_4; 
K=0;
delta=1e-8;
Z_prev=zeros(n+q,q);    
Z_next=zeros(n+q,q);   
tau=4e-7; 
lambda_1=10;    lambda_2=0.1; 
lambda = [lambda_1*ones(n,1); lambda_2*ones(q,1)]; 
G=[D eye(q)];

%Per ogni istante di tempo da 1 a Tmax
while 1
    Z_prev=Z_next; 
    %per ogni sensore da 1 a q
    for i=1:q
        %Local mean 
        A=0;
        for j=1:q
            A=A+Q(i,j)*Z_prev(:,j);
        end
        B=tau * G(i,:)' * ( y(i)-G(i,:)*Z_prev(:,i) );
        arg_sto=A+B;
        %applico l'operatore
        for j=1:(n+q)
            Z_next(j,i)=sto(arg_sto(j), tau*lambda(j));
        end
    end
    K=K+1; 
    %Controllare la condizione di terminazione
    Somma=0; 
    for i=1:q
        Somma=Somma+norm(Z_next(:,i)-Z_prev(:,i))^2;
    end
    
    if (Somma<delta)
        break
    end 
    
end                     

%-----------------PULIZIA DI X-----------------
Ntarget=2; 
x_clean=Z_next(1:n,1);
max_x_vec = maxk(x_clean,Ntarget);
for i=1:n
    if(abs(x_clean(i))<max_x_vec(Ntarget))
        x_clean(i)=0; 
    end
end
Supp_x=find(x_clean)'

A = Z_next(101:125, :);
tol=0.002; 
for i=1:q
    for j=1:q
        if(abs(A(i,j))<tol)
            A(i,j)=0;
        end
    end
end

Supp_a=find(A(:,1))'



