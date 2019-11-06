function [ len ] = arclen( x1, x2 )
%% 
%ARCLEN returns the cumulative arclength of a 2D curve
%   
%   Usage: arclen(x1, x2) Where x1 and x2 are vectors of equal length
%   The function returns a vector of the same length as x1 and x2
%   containing the arclength

len = cumsum( sqrt(diff(x1).^2 + diff(x2).^2) );

end

