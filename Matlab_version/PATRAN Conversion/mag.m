function [ magnitude ] = mag( v )
%Computes magnitude of vector v
    magnitude = sqrt(sum(v.*v));

end

