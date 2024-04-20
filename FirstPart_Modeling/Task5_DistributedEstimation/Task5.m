clear
close all
clc

load distributed_localization_data.mat

n=100; 
q=25;
Q_vec=[Q_4,Q_8,Q_12,Q_18];
Q_name = ["Q_4","Q_8","Q_{12}","Q_{18}"];
Q = zeros(q,q);
for m=1:4
    Q = Q_vec(:,((m-1)*q)+1:m*q);
    K=0;
    delta=1e-8;
    Z_prev=zeros(n+q,q);    
    Z_next=zeros(n+q,q);   
    tau=4e-7; 
    lambda_1=10;    lambda_2=0.1; 
    lambda = [lambda_1*ones(n,1); lambda_2*ones(q,1)]; 
    G=[D eye(q)];

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
    [max_x_vec,Supp_x] = maxk(x_clean,Ntarget);
    for i=1:n
        if(abs(x_clean(i))<max_x_vec(Ntarget))
            x_clean(i)=0; 
        end
    end

    A = Z_next(n+1:n+q, :);
    tol=0.002; 
    for i=1:q
        for j=1:q
            if(abs(A(i,j))<tol)
                A(i,j)=0;
            end
        end
    end

    Supp_a=find(A(:,1))';
    room(Supp_x,Supp_a,K,Q_name(m));
end
figure()
subplot(2,2,1)
plot(digraph(Q_4'))
title(Q_name(1)) 
subplot(2,2,2)
plot(digraph(Q_8'))
title(Q_name(2)) 
subplot(2,2,3) 
plot(digraph(Q_12'))
title(Q_name(3)) 
subplot(2,2,4)
plot(digraph(Q_18'))
title(Q_name(4)) 
sgt = sgtitle('Graphs','Color','red');
sgt.FontSize = 20;

%To do...
%   -> Rappresentazione della 'stanza' (file prof)
%   -> Cambiare Q (rappresentare tramite grafo)
%   -> Analizzare il comportamento in  base a esr(Q)
%   -> Condizione di terminazione sul singolo sensore  
%       --> Stop Gradiente
%       --> Info distribuita pari a quella dell'ultimo passo
%   -> Notare differenze rispetto al caso centralizzato: non ho G!
%           --> Ogni sensore ne conserva una parte
%           --> tau non posso calcolarla a priori
%           --> i lambda vanno ricalibrati
%           --> non posso fare normalize(G)
%   -> Numero di iterazioni!




