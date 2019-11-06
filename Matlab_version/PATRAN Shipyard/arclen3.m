function [ len ] = arclen3( x1, x2, x3 )
%%
%ARCLEN3 returns the cumulative arclength of a 3D curve

len = cumsum( sqrt(diff(x1).^2 + diff(x2).^2 + diff(x3).^2) );

end

