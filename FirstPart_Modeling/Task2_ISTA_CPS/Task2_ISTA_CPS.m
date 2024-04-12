%-----------------------------------------------------
%  Project "Modeling and Control of CPS"
%           TASK 2 - Implementation of ISTA for CPS
%                               
%                                Carlo Migliaccio
%                                28.03.2024
%-----------------------------------------------------
clear all
close all
clc
    
Attack_Type=1;          %1->UNAWARE        2->AWARE

%   hyperparameters
n=10;                           
q=25;                                    %number of sensors
h=2;                                     %sparsity 

Num=40;                                  %number of experiments
NCorrect=0;                              %number of correct support
delta=1e-12;                             
sigma=1e-2; 
K_max = zeros(Num, 1);                   %Time max of the algorithm
dist_x=zeros(Num,1);                     %distance of the solution 
                             
tol=0.1;                                 %are really non-zero components?

eps=1e-8;  
for i=1:1:Num
    C=randn(q,n);                            %Sensing matrix
    x_vero=randn(n,1);                       %State variables

    tau= (norm(C)^(-2))-eps;                    %step size
    lambda=[zeros(n,1)                          %non-sparse 
        ones(q,1)*1/(100*tau)];                 %sparse

    %where are the non-zero components?
    Supp_vero=randperm(q,h);                
    a_vero = zeros(q, 1);     
    %generate in the chosen position
    
    if Attack_Type==2
        %AWARE ATTACK
        y=C*x_vero;
        a_vero(Supp_vero) = y(Supp_vero)*(1/2); 
    else
        %UNAWARE ATTACK
        a_vero(Supp_vero) = sign(randn(1))*unifrnd(1,2, h, 1);
    end

    noise = sigma*randn(q, 1);              
    
    %Rearrangement of the variables for CPS problem
    G=[C eye(q)]; 
    z_vero=[x_vero; a_vero]; 
    y=G*z_vero+noise;

    %------------------------------------------------------
    % Iterative Shrinkage and Thresholding Algorithm (ISTA)
    % Solution of the problem (sparse optimization)
    

    z_calc=zeros((n+q), 1);                         
    %----------------------------------------------------------
    %From here start the ISTA----------------------------------
    while 1
        z_prev=z_calc;
        
        %.........passo k+1 dell'algoritmo
        qi=z_prev+tau*G'*(y-G*z_prev);
        %For each component, apply the S.T. operator
        for j=1:1:(n+q)
            %nuovo elemento di x  
            z_calc(j,1) = sto(qi(j), lambda(j)*tau);
        end
        %.................................................
        K_max(i) = K_max(i)+1;          %Add the step
        if(norm(z_calc-z_prev)<delta)
            break
        end
    end 
    

    %-----------Shrink to 0 elements of a_calc with a certain tol---------
    x_calc = z_calc(1:n);
    a_calc = z_calc(n+1:end);
    
    % ---------------Data cleaning
    for j=(1:q)
        if(abs(a_calc(j))<=tol) 
            a_calc(j)=0; 
        end
    end
    
    Supp_vero=sort(Supp_vero);
    Supp_calc = find(a_calc)';       %Support of the found solution
    dist_x(i)=norm(x_vero-x_calc)^2;       %Distance of x_calc
    %support of the solution found by the ISTA
     if(isequal(Supp_calc, Supp_vero))   NCorrect = NCorrect+1;  end
end

%Results of the computation-------------------------------
Rate=(NCorrect/Num)*100     %Rate of attack detection    
K_max
mean_dist = mean(dist_x)

