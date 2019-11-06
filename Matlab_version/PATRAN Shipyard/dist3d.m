function [ distance ] = dist3d(point1,point2)
%%
%Distance between 2 points in 3D
    distance = sqrt( (point2(1)-point1(1))^2 + (point2(2)-point1(2))^2 + (point2(3)-point1(3))^2 );

end

