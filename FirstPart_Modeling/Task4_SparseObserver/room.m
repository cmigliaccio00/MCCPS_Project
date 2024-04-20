%% April 11, 2024 - S. M. Fosson
%% Dynamic CPS: 3 targets moving in a room
%% The dyanmics is given by matrix A
function room(x_estimated,a_estimated,Tstart,Tmax)
    load("tracking_moving_targets.mat");
  
    n=size(x_estimated);
    n=n(1);
    q=size(a_estimated);
    q=q(1);
    H = 10; %height of the grid (# cells)
    L = 10; %length of the grid (# cells)
    W = 100; %width of a square cell (cm)

    room_grid = zeros(2,n);
    for i=1:n
	    room_grid(1,i) = floor(W/2)+ mod(i-1,L)*W; 
	    room_grid(2,i) = floor(W/2)+ floor((i-1)/L)*W;
    end
    count=0;
    %% sensors localization
    sensors=zeros(q,3); %vector with the position of sensors (repeated)
    numbers=1:q;
    numbers = num2str(numbers');
    for i=1:q
      [~,cell]=maxk(D(i,:),1);
      sensors(i,1)=cell;
    end
    for i=1:q
        index = find(sensors==sensors(i,1));
        switch size(index,1)
            case 1
                sensors(i,2) = room_grid(1,sensors(i,1));
                sensors(i,3) = room_grid(2,sensors(i,1));
            case 2
                sensors(index(1),2) = room_grid(1,sensors(index(1),1))-25;
                sensors(index(1),3) = room_grid(2,sensors(index(1),1));
                sensors(i,2) = room_grid(1,sensors(i,1))+25;
                sensors(i,3) = room_grid(2,sensors(i,1));
            case 3
                %TODO
            case 4
                %TODO
        end
    end
    
    %The definition of initial condition on x_true and a_true are given only for graphical
    %reason. In the reality we don't have the initial condition on the
    %target and also we have to estimate which are the sensors under
    %attack
    x_true = zeros(n,1);
    support_x = [87,23,36];
    x_true(support_x) = 1;
    
    a_true = zeros(q,1);
    support_a = [12,16];
    a_true(support_a) = 1;

    for move = Tstart:Tmax
        x_true=A*x_true;
        target_real = find(x_true);
        target_estimated = find(x_estimated(:,move));

        sensors_under_attack = sensors(support_a);
        support_attack = find(a_estimated(:,move));
        estimated_attack=zeros(size(support_attack,1),2);
        for i = 1:size(support_attack,1)
            estimated_attack(i,:) = sensors(support_attack(i),2:end);
        end
        
        plot(room_grid(1,target_real), room_grid(2,target_real),'square','MarkerSize',9, 'MarkerEdgeColor',1/255*[40 208 220],'MarkerFaceColor',1/255*[40 208 220]);
        hold on
        plot(room_grid(1,target_estimated), room_grid(2,target_estimated),'*','MarkerSize',9, 'MarkerEdgeColor',1/255*[28 55 189]);
        plot(sensors(:,2), sensors(:,3),'o','MarkerSize',9, 'MarkerEdgeColor',1/255*[247 176 240],'MarkerFaceColor',1/255*[247 176 240]);
        plot(room_grid(1,sensors_under_attack), room_grid(2,sensors_under_attack),'o','MarkerSize',9, 'MarkerEdgeColor',1/255*[255 0 0]);
        plot(estimated_attack(:,1), estimated_attack(:,2),'*','MarkerSize',9, 'MarkerEdgeColor',1/255*[255 0 0]);
        
        text(sensors(:,2)-25,sensors(:,3)-50,numbers);
        grid on,        
        legend('Targets','Estimated targets','Sensors','Sensors under attack','Estimated attacks','Location','eastoutside')

        xticks(100:100:1000)
        yticks(100:100:1000)
        xlabel('(cm)')
        ylabel('(cm)')
        axis([0 1000 0 1000])
        axis square
        str = sprintf(' Time = %d', move);
        text(1100,900,str);
        pause(1)
        hold off
    end

end