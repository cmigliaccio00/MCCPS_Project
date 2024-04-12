%-----------------------------------------------------
%  Project "Modeling and Control of CPS"
%           TASK 1 - Implementation of ISTA
%                               
%                                Carlo Migliaccio
%                                22.03.2024
%-----------------------------------------------------
clear all
close all
clc

%   hyperparameters
p=20;                                    %number of state variables  
q=25;                                    %number of sensors
k=2;                                     %sparsity 
Num=20;                                  %number of experiments
NCorrect=0;                              %number of correct support
delta=1e-12;                             
sigma=1e-2; 
                             
tol=0.1;                                 %are really non zero components?
Iterations=zeros(p, 1);                  %number of Iterations required to 
                                         %converge (in the norm-sense)
eps=1e-8;  
for i=1:1:Num
    C=randn(q,p);                            %Sensing matrix
    tau= (norm(C)^(-2))-eps;                 %step size
    lambda=1/(100*tau);                      %hyperparameter for the LASSO
    
    Supp_vero=randperm(p,k);                 %where are the non-zero components?
    x_vero = zeros(p, 1);                    
    x_vero(Supp_vero) = unifrnd(1,2, k, 1);  %generate in the chosen position  
                                             %the (random) state variables
    noise = sigma*randn(q, 1);              
    y=C*x_vero + noise;                      
    %------------------------------------------------------
    % Iterative Shrinkage and Thresholding Algorithm (ISTA)
    % Solution of the problem (sparse optimization)
    
    x_calc=zeros(p, 1);                        %solution of the LASSO 
    
    %----------------------------------------------------------
    %From here start the ISTA----------------------------------
    while 1
        x_prev=x_calc;
        qi=x_prev+tau*C'*(y-C*x_prev);
        %For each component, apply the S.T. operator
        for j=1:1:p
            %nuovo elemento di x  
            x_calc(j,1) = sto(qi(j), lambda*tau);
        end
        if(norm(x_calc-x_prev)<delta)
            break
        end
        Iterations(i)=Iterations(i)+1; 
    end 
    
    [x_calc x_vero]

    %-----------Shrink to 0 elements of x_calc with a certain tol---------
    for j=1:p
        if(abs(x_calc(j))<=tol) 
            x_calc(j)=0; 
        end
    end
    
    [x_calc x_vero]
    
    %support of the solution found by the ISTA
    Supp_calc = find(x_calc)';
    
    if(isequal(Supp_calc, Supp_vero))   NCorrect = NCorrect+1;  end
end

Rate = (NCorrect/Num)*100
Stat = [min(Iterations) max(Iterations) mean(Iterations)]
plot((1:20), Iterations)
