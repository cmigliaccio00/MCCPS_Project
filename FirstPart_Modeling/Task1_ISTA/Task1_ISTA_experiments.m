
%% First Graph: SUPPORT OF RECOVERY RATE  VARYING q
clear all
close all
clc

q_test = 5:1:20;
q_test=q_test';
q_dim = size(q_test);
lim_max = q_dim(1);

% Rate of attack detection varying q

colors = ["red", "blue", "green", "magenta"];
p=20;                                    %number of state variables  
%q=25;  
i_col=1; 
%number of sensors

coef_l=[30 100 1000];

for kk=1:size(coef_l,2)
    k=2;                                     %sparsity 
    Num=20;                                  %number of experiments
                          %number of correct support
    delta=1e-12;                             
    sigma=1e-2; 
                                 
    tol=0.1;                                 %are really non zero components?
    Iterations=zeros(lim_max, 1);                  %number of Iterations required to 
                                             %converge (in the norm-sense)
    eps=1e-8;  
    
    Supp_vero=randperm(p,k);                 %where are the non-zero components?
    x_vero = zeros(p, 1);                    
    x_vero(Supp_vero) = unifrnd(1,2, k, 1);  %generate in the chosen position  
                                         %the (random) state variables
    
    Rates_q = zeros(lim_max,1);
    %-----------------Variation of q---------------------------
    for ii=1:lim_max
        NCorrect=0;  
        
        q = q_test(ii);
        %Faccio per ogni numero di sensori un certo numero di esperimenti
        for jj=1:20
            %   hyperparameters
            C=randn(q,p); 
            noise = sigma*randn(q, 1);    
            tau= (norm(C)^(-2))-eps;                 %step size
            lambda=1/(coef_l(kk)*tau);                      %hyperparameter for the LASSO          
            y=C*x_vero + noise; 
                  
            %----------------QUI CONGELO I DATI----------------------------------------
            
            
            %------------------------------------------------------
            % Iterative Shrinkage and Thresholding Algorithm (ISTA)
            % Solution of the problem (sparse optimization)
            
            x_calc=zeros(p, 1);                      %initial condition
            
            
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
                Iterations(ii)=Iterations(ii)+1; 
            end 
            
            %[x_calc x_vero];
            
            %-----------Shrink to 0 elements of x_calc with a certain tol---------
            for j=1:p
                if(abs(x_calc(j))<=tol) 
                    x_calc(j)=0; 
                end
            end
            
            %[x_calc x_vero];
            
            %support of the solution found by the ISTA
            Supp_calc = sort(find(x_calc))';
            Supp_vero = sort(Supp_vero);
            
            if(isequal(Supp_calc, Supp_vero))   NCorrect = NCorrect+1;  end
            
        end
        NCorrect;
        Rate = (NCorrect/20);
    
        %Tasso di successo per q_test(ii) sensori
        Rates_q(ii) = Rate;
        info_exp = [q Rate];
        %Stat = [min(Iterations) max(Iterations) mean(Iterations)];
        %plot((1:20), Iterations)
    end
    
    %Costruzione del grafico
    scatter(q_test, Rates_q, "square", 'MarkerEdgeColor',colors(i_col))
    line(q_test, Rates_q,'Color', colors(i_col))
    i_col = i_col+1;
    %title("\textbf{Support recovery rate varying sensors number} ", 'Interpreter','latex', 'FontSize',15)
    xlabel('Number of sensors','FontSize',13, 'Interpreter', 'latex');
    xticks(q_test)
    ylabel("Recovery rate", 'FontSize',13, 'Interpreter', 'latex');
    hold on
    
end

legend('', '$\lambda_1=\frac{1}{30\tau}$', '', '$\lambda_2=\frac{1}{100\tau}$', '', '$\lambda_3=\frac{1}{1000\tau}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 13)
leg = legend; 
leg.Position(4)=0.25; 

%% Second Graph: Convergency time varying lambda

clear all
close all
clc

p=20;                                    %number of state variables  
%q=25;  
%number of sensors
k=2;                                     %sparsity 
Num=20;                                  %number of experiments
                      %number of correct support
delta=1e-12;                             
sigma=1e-2; 
                             
tol=0.1;                                 %are really non zero components?
q = 20;

Supp_vero=randperm(p,k);                 %where are the non-zero components?
x_vero = zeros(p, 1);                    
x_vero(Supp_vero) = unifrnd(1,2, k, 1);  %generate in the chosen position  
                                         %the (random) state variables
C=randn(q,p); 

eps=[1e-3; 1e-8];  
col=["red", "blue"];
noise = sigma*randn(q, 1);              
y=C*x_vero + noise;

for ee=1:size(eps,1)
    tau= (norm(C)^(-2))-eps(ee);                            %step size
    
    lambda_val = sort(1./(linspace(40,1000,10).*tau));
    lim_max = size(lambda_val);
    lim_max = lim_max(2);

    
    Iterations=zeros(lim_max, 1);                  %number of Iterations required to 
                                             %converge (in the norm-sense)
    for ii=1:lim_max
        % Faccio un certo numero di esperimenti
        %------------------------------------------------------
        % Iterative Shrinkage and Thresholding Algorithm (ISTA)
        % Solution of the problem (sparse optimization)
        
        x_calc=zeros(p, 1);                      %initial condition
        lambda=lambda_val(ii);
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
            Iterations(ii)=Iterations(ii)+1; 
        end 
    
        %-----------Shrink to 0 elements of x_calc with a certain tol---------
        for j=1:p
            if(abs(x_calc(j))<=tol) 
                x_calc(j)=0; 
            end
        end
    end
    
    hold on
    semilogx(lambda_val, Iterations, 'Color', col(ee));
    scatter(lambda_val, Iterations, 'd', 'Color',col(ee));
    [lambda_val' Iterations]

    %ylim([0 500]
    ylabel('Convergency time (\#Iterations)','Interpreter', 'latex', 'FontSize',13);
    xlabel('$\lambda$','Interpreter', 'latex', 'FontSize',14);
    %title('\textbf{Convergence time varying} $\lambda$', 'Interpreter', 'latex', 'FontSize',14)
end



legend(['$\tau_1=\Vert C \Vert_2^{-2}-10^{-3}$'], '', ...
    '$\tau_2=\Vert C \Vert_2^{-2}-10^{-8}$','', 'Interpreter','latex','FontSize',12)

leg = legend; 
leg.Position(4)=0.20; 

%% Third Graph: Convergency time varying tau

clear all
close all
clc

p=20;                                    %number of state variables  
%q=25;  
%number of sensors
k=2;                                     %sparsity 
Num=20;                                  %number of experiments
                      %number of correct support
delta=1e-12;                             
sigma=1e-2; 
                             
tol=0.1;                                 %are really non zero components?
q = 20;


Supp_vero=randperm(p,k);                 %where are the non-zero components?
x_vero = zeros(p, 1);                    
x_vero(Supp_vero) = unifrnd(1,2, k, 1);  %generate in the chosen position  
                                         %the (random) state variables
C=randn(q,p); 

eps=[1e-3; 1e-8];  
col=["red", "blue"];
noise = sigma*randn(q, 1);              
y=C*x_vero + noise;

tau0=(norm(C))^(-2)-1e-8;
lambda0 = 1/(100*tau0);
const = tau0*lambda0;

tau_val = linspace(0.01, 0.02, 30);
lim_max = size(tau_val);
lim_max = lim_max(2);

Iterations=zeros(lim_max, 1);                  %number of Iterations required to 
                                         %converge (in the norm-sense)
for ii=1:lim_max
    % Faccio un certo numero di esperimenti
    %------------------------------------------------------
    % Iterative Shrinkage and Thresholding Algorithm (ISTA)
    % Solution of the problem (sparse optimization)
    
    x_calc=zeros(p, 1);                      %initial condition
    tau = tau_val(ii);
    lambda=const/tau;
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
        Iterations(ii)=Iterations(ii)+1; 
    end 

    %-----------Shrink to 0 elements of x_calc with a certain tol---------
    for j=1:p
        if(abs(x_calc(j))<=tol) 
            x_calc(j)=0; 
        end
    end
end


line(tau_val, Iterations,'Color', 'blue');
hold on
scatter(tau_val, Iterations, 'square', 'MarkerEdgeColor','blue');
[tau_val' Iterations]

%ylim([0 500]
ylabel('Convergency time (\#Iterations)','Interpreter', 'latex', 'FontSize',13);
xlabel('$\tau$','Interpreter', 'latex', 'FontSize',14);
%title('\textbf{Convergence time varying} $\tau$', 'Interpreter', 'latex', 'FontSize',14)
