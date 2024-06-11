clear
close all
clc 

%tracking_lost_x = zeros(100,1);
%tracking_lost_a = zeros(100,2);
for repetition=1:100
%% INITIALIZATION
load tracking_moving_targets.mat    %load A,D,Y

p=100;          %number of the cells
q=25;           %number of sensors

lambda1=10; lambda2=20; 
lambda=[lambda1*ones(p,1); lambda2*ones(q,1)];
eps=1e-8;                            
G = [D eye(q)];                     % augmented sensing matrix
G=normalize(G);                     
tau= (norm(G)^(-2))-eps;            % step size
mes_x = zeros(p,1);                 % save inside the matrix the estimate of x at each step
mes_a = zeros(q,1);                 % save inside the matrix the estimate of a at each step
z_hat = [zeros(p,1); zeros(q,1)];   % Initial state observer

x_true=zeros(p,1);
a_true=zeros(q,1);
%% OPTIONS
Ntarget = 3;
Nattack = 25;
aware = 1;
change_sensors =0;
Tmax = 100;
%% CHECKS
if change_sensors==1 && aware==0
    return
end

if change_sensors && mod(Tmax,25)~=0
    return
elseif ~aware && Tmax~=50
    return
end
if ~aware && Ntarget~=3
    'The number of target must be 3'
    return
end
if ~aware && Nattack~=2
    'The number of sensors under attack must be 3'
    return
end

%% ALGORITHM
if aware    
    %generate inital condition on the position
    support_x_true = randperm(p);
    support_x_true = support_x_true(1:Ntarget); % I consider 3 targets
    x_true(support_x_true) = 1;


    %generate attack support
    tmp_supp_a = randperm(q);
    tmp_supp_a = tmp_supp_a(1:Nattack);
    
    sigma=1e-2; 
    noise = sigma*randn(q, 1);
    
    if change_sensors
        rep = Tmax/25;
        support_a_true=zeros(rep,Nattack);
        for i=1:rep        
            for j=((i-1)*25)+1:25*i %portion of Y to be filled
                x_true=A*x_true;
                Y(:,j) = D*x_true + a_true + noise;
                a_true(tmp_supp_a)=0.5*Y(tmp_supp_a,j);
            end
            %generate a new attack support
            a_true=zeros(q,1);
            support_a_true(i,:)=tmp_supp_a;
            tmp_supp_a = randperm(q);
            tmp_supp_a = tmp_supp_a(1:Nattack);
        end
    else
        for i=1:Tmax
            x_true=A*x_true;
            Y(:,i) = D*x_true + a_true + noise;
            a_true(tmp_supp_a)=0.5*Y(tmp_supp_a,i); 
        end
        support_a_true(1,:)=tmp_supp_a;
    end
else
    %The definition of initial condition on x_true and a_true are given 
    % only for graphical reason. In real cases we don't have the initial 
    % condition on the target and also we have to estimate which are the 
    % sensors under attack
    x_true = zeros(p,1);
    support_x_true = [87,23,36];
    x_true(support_x_true) = 1;
    
    a_true = zeros(q,1);
    support_a_true = [12,16];
    a_true(support_a_true) = 1;
    
end
  

%------------------------SPARSE OBSERVER------------------------
z_hat_plus=zeros(p+q,1);
for k=1:Tmax
  
    arg=z_hat + tau*G'*(Y(:,k)-G*z_hat);                  %arg of STO
   
    for i=1:(p+q)
        z_hat_plus(i) = sto(arg(i), tau*lambda(i));        
    end
   
    z_hat=[A*z_hat_plus(1:p); z_hat_plus(p+1:end)];
    
    x_hat=A*z_hat_plus(1:p);
    a_hat=z_hat_plus(p+1:end);

    mes_x = [mes_x x_hat];
    mes_a = [mes_a a_hat];
end
mes_x = mes_x(:,2:end); % we delete the first colum that is all 0
mes_a = mes_a(:,2:end); % we delete the first colum that is all 0

%-------------------------Data cleaning--------------------------
for j=1:Tmax
    %cleaning x_hat
    max_x_vec = maxk(abs(mes_x(:,j)),Ntarget);
    for i=1:p 
        if(abs(mes_x(i,j))<max_x_vec(end))
            mes_x(i,j)=0; 
        end
    end

    %cleaning a_hat
    max_a = max(abs(mes_a(:,j)));
    for i=1:q
        if(abs(mes_a(i,j))<max_a*50/100)
            mes_a(i,j)=0; 
        end
    end
end
%Plot the room
[x_correct_estimate,a_correct_estimate, total_estimate] = room(mes_x,mes_a,1,Tmax,support_x_true,support_a_true,change_sensors);
%%PLOT THE CONVERGENCE OF THE OBSERVER
time_interval = zeros(1,Tmax);

flag_x = 0;
flag_a = 0;
for k=1:Tmax-1
    %time_interval(k) = k;
    if change_sensors
        index = fix(k/25)+1;
        if a_correct_estimate(k+1) == 1 && a_correct_estimate(k) == 0
            conv_a(repetition,index) = k+1;
        elseif a_correct_estimate(k+1) == 0 && a_correct_estimate(k) == 1 && k ~=24
            conv_a(repetition,index) = 0;
        end
        if a_correct_estimate(k+1) == 1 && a_correct_estimate(k) == 0 && flag_a
            tracking_lost_a(repetition,index) = 1;
            flag_a=0;
        end
        if a_correct_estimate(k+1) == 0 && a_correct_estimate(k)==1 && k~=24
            flag_a=1;
        end
    else
        if a_correct_estimate(k+1) == 1 && a_correct_estimate(k) == 0
            conv_a(repetition,1) = k+1;
        elseif a_correct_estimate(k+1) == 0 && a_correct_estimate(k) == 1
            conv_a(repetition,1) = 0;
        end
       if a_correct_estimate(k+1) == 1 && a_correct_estimate(k) == 0 && flag_a
            tracking_lost_a(repetition,1) = 1;
            flag_a=0;
        end
        if a_correct_estimate(k+1) == 0 && a_correct_estimate(k)==1 && k~=24
            flag_a=1;
        end
    end
    if x_correct_estimate(k+1) == 1 && x_correct_estimate(k) == 0
        conv_x(repetition) = k+1;
    elseif x_correct_estimate(k+1) == 0 && x_correct_estimate(k) == 1
        conv_x(repetition) = 0; 

    end
    if x_correct_estimate(k+1) == 1 && x_correct_estimate(k) == 0 && flag_x
        tracking_lost_x(repetition) = 1;
        flag_x=0;
    end
    if x_correct_estimate(k+1) == 0 && x_correct_estimate(k)==1
        flag_x=1;
    end
    
    if total_estimate(k+1) == 1 && total_estimate(k) == 0
        conv_obs(repetition) = k+1;
    elseif total_estimate(k+1) == 0 && total_estimate(k) == 1
        conv_obs(repetition) = 0;
    end
end
% 
% figure
% grid on
% subplot(3,1,1)
% plot(time_interval,a_correct_estimate,'LineWidth',2,'Color','r');
% ylim([-0.5,1.5]);
% xlim([1,Tmax]);
% title("Convergence of the attack $a$","Interpreter","latex","FontSize",15,"Color",'r');
% xlabel("iteration","Interpreter","latex");
% ylabel("prediction","Interpreter","latex");
% yticklabels({'','Uncorrect', '','Correct'});
% subplot(3,1,2)
% grid on
% plot(time_interval,x_correct_estimate,'LineWidth',2,'Color','b');
% ylim([-0.5,1.5]);
% xlim([1,Tmax]);
% title("Convergence of the state $x$","Interpreter","latex","FontSize",15,"Color",'b');
% xlabel("iteration","Interpreter","latex");
% ylabel("prediction","Interpreter","latex");
% yticklabels({'','Uncorrect', '','Correct'});
% subplot(3,1,3)
% grid on
% plot(time_interval,total_estimate,'LineWidth',2,'Color','g');
% ylim([-0.5,1.5]);
% xlim([1,Tmax]);
% title("Convergence of the Sparse Observer","Interpreter","latex","FontSize",15,"Color",'g');
% xlabel("iteration","Interpreter","latex");
% ylabel("prediction","Interpreter","latex");
% yticklabels({'','Uncorrect', '','Correct'});
repetition
close all
end