function [ scaledplan ] = scale_linesplan( plan, scalefactors )
%% Done...
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    scaledplan = plan;
    scaledplan.L = plan.L*scalefactors.L;
    scaledplan.B = plan.B*scalefactors.B;
    scaledplan.D = plan.D*scalefactors.D;
    scaledplan.T = plan.T*scalefactors.T;
    
    for o=plan.ords:-1:1
        scaledplan.coordinates{1,o} = plan.coordinates{1,o}*scalefactors.L;
        scaledplan.coordinates{2,o} = plan.coordinates{2,o}*scalefactors.B;
        scaledplan.coordinates{3,o} = plan.coordinates{3,o}*scalefactors.D;
    end

end

