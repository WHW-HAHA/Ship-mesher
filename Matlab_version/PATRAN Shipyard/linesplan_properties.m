function [ properties ] = linesplan_properties( plan )
% Computes basic properties of a linesplan

    wls = 10;

    y = linspace(0.01,plan.T,wls);
    for o=1:plan.ords % for each cross-section
        clear points
        x(o) = plan.coordinates{1,o};
        points(1,:) = plan.coordinates{2,o};  % z
        points(2,:) = plan.coordinates{3,o};  % y
        for i=2:length(points(2,:))
            if plan.coordinates{3,o}(i) == plan.coordinates{3,o}(i-1)
                points(2,i) = points(2,i) + 0.001*i;
            end
        end
        frame{o} = pchip(points(2,:),points(1,:));
        coords{o} = vertcat(y,fnval(frame{o},y));          % segment coordinates per section
    end


    for o=1:plan.ords % for each cross-section
        for w=1:length(y)
            if y(w)>plan.coordinates{3,o}(1)
                waterline{w}(o) = coords{o}(2,w);
            else 
                waterline{w}(o) = 0;
            end 
        end
    end

    for w=1:length(y)
       A_w(w) = 2 * simps(x,waterline{w}(1,:));
       I_t(w) = 2/3 * simps(x,waterline{w}(1,:).^3);
    end
    
    properties.Disp = simps(y,A_w);
    properties.CB = properties.Disp/(plan.L*plan.B*plan.T);
    properties.LB = plan.L/plan.B;
    properties.BT = plan.B/plan.T;
end

