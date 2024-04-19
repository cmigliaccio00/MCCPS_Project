%-----------------------------------------------------
%  Project "Modeling and Control of CPS"
%           TASK 3 - Localization under sparse attacks
%                               
%                                Carlo Migliaccio
%                                28.03.2024
% 15.04.2024
%-----------------------------------------------------
clear
close all
clc

%load the provided data D, y
load localization.mat   

%   hyperparameters
p=100;                                   %number of cells of the grid
q=25;                                    %number of sensors

delta=1e-15;                             %used in the stop condition                           
                             
tol=0.1;                                 %are really non-zero components?

lambda1=10;   lambda2=20;               %weights of the lasso
eps=1e-8; 
lambda=[lambda1*ones(p,1); lambda2*ones(q,1)];  %sparse
G=[D eye(q)];                                   %augmented sensing matrix
G=normalize(G);
tau= (norm(G)^(-2))-eps;                        %step size

z_calc=zeros((p+q), 1);   %z_0                      
%----------------------------------------------------------
%From here start the ISTA----------------------------------
while 1
    z_prev=z_calc;
    %.........passo k+1 dell'algoritmo
    qi=z_prev+tau*G'*(y-G*z_prev);
    %For each component, apply the S.T. operator
    for j=1:1:(p+q)
        %nuovo elemento di x  
        z_calc(j,1) = sto(qi(j), lambda(j)*tau);
    end
    %.................................................
    if(norm(z_calc-z_prev)<delta)
        break
    end
end 

tol=4;
for i=1:(p+q)
    if(abs(z_calc(i))<tol)    
        z_calc(i)=0;
    end
end

%Separate x and a
x_calc = z_calc(1:p);
a_calc = z_calc(p+1:end);
Supp_x = find(x_calc)'
Supp_a = find(a_calc)'


   


