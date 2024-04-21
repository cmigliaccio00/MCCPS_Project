function room(x_estimated,a_estimated,Tstart,Tmax,support_x,support_a,change_sensor)
    load tracking_moving_targets.mat;
  
    n=size(x_estimated,1);
    q=size(a_estimated,1);
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
                sensors(index(1),2) = room_grid(1,sensors(index(1),1))-25;
                sensors(index(1),3) = room_grid(2,sensors(index(1),1))-20;
                sensors(index(2),2) = room_grid(1,sensors(index(1),1))+25;
                sensors(index(2),3) = room_grid(2,sensors(index(1),1))-20;
                sensors(i,2) = room_grid(1,sensors(i,1));
                sensors(i,3) = room_grid(2,sensors(i,1))+20;
            case 4
                sensors(index(1),2) = room_grid(1,sensors(index(1),1))-25;
                sensors(index(1),3) = room_grid(2,sensors(index(1),1))-25;
                sensors(index(2),2) = room_grid(1,sensors(index(1),1))+25;
                sensors(index(2),3) = room_grid(2,sensors(index(1),1))-25;
                sensors(index(3),2) = room_grid(1,sensors(index(1),1))+25;
                sensors(index(3),3) = room_grid(2,sensors(index(1),1))+25;
                sensors(i,2) = room_grid(1,sensors(i,1))-25;
                sensors(i,3) = room_grid(2,sensors(i,1))+25;    
        end
    end
    
    x_true=zeros(n,1);
    x_true(support_x) = 1;
    row=1;
    for move = Tstart:Tmax
        x_true=A*x_true;
        target_real = find(x_true);
        target_estimated = find(x_estimated(:,move));
        if change_sensor
            sensors_under_attack=support_a(row,:);
            if mod(move,25)==0
                row=row+1;
            end
        else
            sensors_under_attack = sensors(support_a(1,:));
        end        
        sensors_under_attack=sensors(sensors_under_attack,:);
        support_attack = find(a_estimated(:,move));      
        estimated_attack=zeros(size(support_attack,1),2);
        for i = 1:size(support_attack,1)
            estimated_attack(i,:) = sensors(support_attack(i),2:end);
        end
        
        
        plot(sensors(:,2), sensors(:,3),'o','MarkerSize',10, ...
            'MarkerEdgeColor',1/255*[247 176 240], ...
            'MarkerFaceColor',1/255*[247 176 240]);
        hold on
        plot(sensors_under_attack(:,2), sensors_under_attack(:,3), ...
            'o','MarkerSize',15, ...
            'MarkerEdgeColor',1/255*[255 0 0]);
        plot(estimated_attack(:,1), estimated_attack(:,2), ...
            '*','MarkerSize',10, ...
            'MarkerEdgeColor',1/255*[255 0 0]);
        text(sensors(:,2)-20,sensors(:,3)-30,numbers, ...
            'FontSize',8,'FontWeight','bold');
        plot(room_grid(1,target_real), room_grid(2,target_real), ...
            'square','MarkerSize',10, ...
            'MarkerEdgeColor',1/255*[40 208 220], ...
            'MarkerFaceColor',1/255*[40 208 220]);
        plot(room_grid(1,target_estimated), room_grid(2,target_estimated), ...
            '*','MarkerSize',10, ...
            'MarkerEdgeColor',1/255*[28 55 189]);
        
        grid on,        
        legend( 'Sensors', ...
                'Sensors under attack', ...
                'Estimated attacks', ...
                'Targets', ...
                'Estimated targets', ...
                'Location', ...
                'eastoutside')

        xticks(100:100:1000)
        yticks(100:100:1000)
        xlabel('(cm)')
        ylabel('(cm)')
        axis([0 1000 0 1000])
        axis square
        str = sprintf(' Time = %d', move);
        text(1100,900,str);
        pause(0.5)
        hold off
    end

end