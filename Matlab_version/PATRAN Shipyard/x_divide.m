function [ points ] = x_divide( curve, spacing, evaluations )
%Evaluates a spline curve at a number of x positions
    t = linspace(curve.breaks(1),curve.breaks(end),evaluations);
    x = fnval(curve,t);
    t2 = interp1(x(1,:),t,spacing);
    points = fnval(curve,t2);
    points(1,:) = spacing;
end

