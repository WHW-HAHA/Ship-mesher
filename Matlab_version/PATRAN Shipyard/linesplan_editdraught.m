function [ plan ] = linesplan_editdraught( plan,pansize )
% Returns a new linesplan structure at a different draught

L = plan.coordinates(1,end); % The most precise formulation of L is needed to avoid losing the last ordinate in the interpolation 
L = L{:};   % matlab cell datatype is a bit odd... 
%% Spline the cross-sections
coords = cell(plan.ords,1);
frame = cell(1,plan.ords);
divs = 0;
for o=1:round(plan.ords/2) % for each cross-section
    points(1,:) = plan.coordinates{2,o};  % y
    points(2,:) = plan.coordinates{3,o};  % z
    frame{o} = cscvn(points);
    
    len = max(arclen(points(1,:),points(2,:)),frame{o}.breaks(end));
    divs = max(divs,round(len(end)/pansize)+1);        % number of segments along spline
    
    s = linspace(frame{o}.breaks(1),frame{o}.breaks(end),divs);
    coords{o} = vertcat(ones(1,length(s))*plan.coordinates{1,o},fnval(frame{o},s)); % segment coordinates per section
    
    clear points
end

%% Define the keel
keelline = zeros(2,plan.ords);
for o=1:plan.ords   % last x and z coordinate are assumed to be on the keel
    keelline(1,o) = plan.coordinates{1,o};
    keelline(2,o) = plan.coordinates{3,o}(1);
end
keelx = keelline(1,:);                   % keel x vector for interpolation
keelz = keelline(2,:);                   % keel z vector for interpolation
keelspline = pchip(keelx,keelz);         % used pchip because splines give wrong results

%% Collect lengthwise coordinates
% We want to draw a spline through every top coordinate per section, 
% and the one underneath and so on. Not every section has the same amount
% of points. To prevent the splines from stopping somewhere on a
% cross-section halway through the ship, extra points will be added on the
% keel. This makes sure the splines start and end either on the keel or on
% the last cross-section. 
coords2 = coords; % make a copy of the old coordinate distribution
for o=1:plan.ords-1
    extra = abs(length(coords2{o+1})-length(coords2{o}))-1; % determine the number of extra points needed
    if extra>0 
        dx = (plan.coordinates{1,o+1}-plan.coordinates{1,o})/(extra+1); % distance between the extra points
    end
    for i=1:extra % interpolate the extra points on the keel and add them to the coordinate list
        if o<plan.ords/2   ///////////// plan.ords = num_section
            x = plan.coordinates{1,o}+i*dx;
            z = fnval(keelspline,x); 
            coords{o} = horzcat([x 0 z]', coords{o});
        else
            x = plan.coordinates{1,o+1}-i*dx;
            z = fnval(keelspline,x);
            coords{o+1} = horzcat([x 0 z]', coords{o+1});
        end
    end
end

% now the coordinates + extra coordinates are listed per cross-section in
% 'coords'. They will be reorganised into 'points' where they form the
% basis for the lenghtwise splines.

for o=1:plan.ords % get the number of lengthwise splines needed
   m(o) = length(coords{o});
end

for i=0:max(m)-1
    c=1;
    clear xyz
    for o=1:plan.ords
        if size(coords{o},2)>i % if there are points to add
            xyz(:,c) = coords{o}(:,end-i); % add the point to xyz
            c=c+1;
        end
    end
%     plot3(xyz(1,:),xyz(2,:),xyz(3,:),'d','Color',[64 104 225]/255)
%     hold on
    if size(xyz,2)>1
	
        points{i+1} = xyz; %put xyz in points if it has a value
        strook{i+1} = cscvn(points{i+1});
    end
end

%% New cross sections
cuts = cosspace(plan.coordinates{1,1},L,30,0.5*pi);%keelx;% linspace(plan.coordinates{1,1},L,30);
keelpoints = zeros(3,length(cuts));
keelpoints(1,:) = cuts;
keelpoints(3,:) = fnval(keelspline,cuts);

for i=1:length(points)
   coords{i} = x_divide(strook{i},cuts,2*length(points));
end

for c=1:length(cuts)
    vec = zeros(3,length(points));
    for i=1:length(points)
        vec(:,i) = coords{i}(:,c);
    end
    vec(:,any(isnan(vec))) = [];
    vec(2,vec(2,:)<0) = 0;
    if not(isempty(vec)) && dist3d(vec(:,end),keelpoints(:,c))<(pansize/2)
        vec(:,end) = keelpoints(:,c);
    else
        vec(:,end+1) = keelpoints(:,c);
    end
    cross{c} = vec;
end

%% Cut off draft
c = 1;
while c<=length(cross)
   if cross{c}(3,end) > plan.T
        cross(c) = [];
        c = c - 1;
   else
      for i=size(cross{c},2):-1:1
         if cross{c}(3,i) > plan.T
            breadth_at_draft = interp1([cross{c}(3,i+1) cross{c}(3,i)],[cross{c}(2,i+1) cross{c}(2,i)],plan.T);
            for j=i-1:-1:1
                cross{c}(:,j) = [];
            end
            cross{c}(:,1) = [cross{c}(1,2) breadth_at_draft plan.T];
            break
         end
      end
   end
   c = c + 1;
end

%% Store
for c=1:length(cross)
    x = cross{c}(1,1);
    coordinates{1,c} = x;
    coordinates{2,c} = cross{c}(2,end:-1:1);
    coordinates{3,c} = cross{c}(3,end:-1:1);
end

plan.coordinates = coordinates;
plan.ords = length(coordinates);

end

