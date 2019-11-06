function [ info ] = linesplan_panelsize(plan,KG,k_xx,k_yy,waterdepth)
%Computes panelsize and other information about a linesplan
% 
% 
% Revised:
%   Yijun 01052016: GM_l value
%                   kzz

    wls = 10; 
    g = 9.81;
    

    y = linspace(0.01,plan.T,wls);
    for o=1:plan.ords % for each cross-section
        clear points
        x(o) = plan.coordinates{1,o};
        points(1,:) = plan.coordinates{2,o};  % y
        points(2,:) = plan.coordinates{3,o};  % z
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

    for w=1:length(y) % y 就是 y
       A_w(w) = 2 * simps(x,waterline{w}(1,:)); % x 就是 x
       Lcf(w) = 2 * simps(x,waterline{w}(1,:).*x)/A_w(w); % LCF_x
       I_t(w) = 2/3 * simps(x,waterline{w}(1,:).^3);
    end
    
    displacement = simps(y,A_w);
    CB = displacement/(plan.L*plan.B*plan.T);
    BM = I_t(end)/displacement;
    KB = simps(y,A_w.*y)/displacement;
    LCB = simps(y,A_w.*(Lcf-plan.L/2))/displacement; %Yijun
    GM = KB + BM - KG;
    fprintf(['LCB =  ' num2str(LCB) '\n']);%Yijun
    fprintf(['KB =  ' num2str(KB) '\n']);  %Yijun

    
    omega_heave = sqrt(A_w(end)*g/displacement);
    omega_roll  = sqrt(GM*g/(k_xx^2));
    
    I_l  = 2*simps(x,waterline{end}(1,:).*(x.^2)) - (Lcf(end))^2*A_w(end);
    BM_l = I_l/displacement;
    GM_l = KB + BM_l - KG;

    omega_pitch = sqrt(GM_l*g/(k_yy^2));

    om = [omega_pitch omega_roll omega_heave];
    wavelength = dispersion_relation(max(om)*2,waterdepth);
    pansize = min(wavelength/4,plan.T/4);
    info.panelsize = round(pansize,1);
    info.Aw = A_w(end);
    info.Lcb = LCB;
    info.GM_l = GM_l;
    info.GM_t = GM;
    info.displacement = displacement;
    info.I_t = I_t(end);
    info.I_l = I_l;   
    info.wn_3 = omega_heave;
    info.wn_4 = omega_roll;
    info.wn_5 = omega_pitch;
end

