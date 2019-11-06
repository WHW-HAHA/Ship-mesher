function mesh_linesplan2(plan,pansize,pandist)
%% Prepare by cutting off draft
plan = linesplan_editdraught(plan,pansize);
assignin('base','plan',plan)
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

divs = 0;
for o=plan.ords:-1:round(plan.ords/2) % for each cross-section
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
        if o<plan.ords/2
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
        points{i+1} = xyz; % put xyz in points if it has a value
        strook{i+1} = cscvn(points{i+1});
    end

end
 
%% Plot Lengthsplines
% figure(6)
% hold on
% for i=1:length(points)
%     curve = cscvn(points{i});
%     pp = fnplt(curve);
%     plot3(pp(1,:),pp(2,:),pp(3,:),'Color',[154 205 49]/255,'LineWidth',1.5);
% end
% title('Bulkcarrier, 187.00 x 29.00 x 11.95 (13.00) meter.')

%% New cross sections
clear coords
figure
hold on
cuts = cosspace(plan.coordinates{1,1},L,ceil(L/pansize),pandist);
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
    vec(3,vec(3,:)<0) = 0;
    cross{c} = vec;
end

%% Extra points on keel
% now that there are new interpolated cross sections extra points have to
% be specified again. But now such that there is enough points to make
% quadrilateral panels
cross2 = cross;
for c=1:length(cuts)-1
    extra = abs(size(cross2{c+1},2)-size(cross2{c},2));
    if extra>0 
        dx = (cuts(c+1)-cuts(c))/(extra+1);
    end
    t=1;
    s=length(cross{o});
    for i=1:extra
        if c<round(length(cuts)/2)
            x = cuts(c)+t*dx;
            z = fnval(keelspline,x); 
            cross{c} = horzcat(cross{c},[x 0 z]');
        elseif c>round(length(cuts)/2)
            x = cuts(c+1)-t*dx;
            z = fnval(keelspline,x);
            cross{c+1} = horzcat(cross{c+1},[x 0 z]');
        end
        t=t+1;
    end
end

%% Plot hull points
% figure(6)
% hold on
% for i=1:length(cross)
%    plot3(cross{i}(1,:),cross{i}(2,:),cross{i}(3,:),'d','Color',[64 104 225]/255,'MarkerSize',2)
% end
% view(3)
% axis equal
% hold off

%% Make PATRAN structure
% write list of coordinates
% while bookkeeping their number
panelnumber = 1;
patr.crd = [];

% add hull
for c=1:length(cross)
    patr.crd = vertcat(patr.crd,(cross{c}'));
    cnums{c} = panelnumber:panelnumber+size(cross{c},2)-1;
    panelnumber = panelnumber+size(cross{c},2);    
end

 patr.pan = [];
for c=1:length(cross)
    if c<round(length(cross)/2)
        for i=1:length(cnums{c})-1
            patr.pan = vertcat(patr.pan,[cnums{c}(i) cnums{c+1}(i) cnums{c+1}(i+1) cnums{c}(i+1)]);
        end
    elseif c>round(length(cross)/2)
        for i=1:length(cnums{c})-1
            patr.pan = vertcat(patr.pan,[cnums{c}(i+1) cnums{c-1}(i+1) cnums{c-1}(i) cnums{c}(i)]);
        end
    end
end

%% Lid 
d = coords{1,1}(3);
for i=1:length(coords{1,1}(1,:))-1 %for each row
    %determine x and y coordinates (x scalar, y vector)
    x = coords{1,1}(1,i);
    xplus = coords{1,1}(1,i+1);
    b = coords{1,1}(2,i); %breadth
    bplus = coords{1,1}(2,i+1);
    if isnan(b);b=0;end
    if isnan(bplus);bplus=0;end
    p = ceil(b/pansize)+1;
    y = linspace(0,b,p);
    yplus = linspace(0,bplus,p);

    if b>eps(b) && bplus>eps(bplus) %square panels can be drawn
        for j=1:p-1
             l = length(patr.crd);
             patr.crd = vertcat(patr.crd,[x y(j) d],[x y(j+1) d],[xplus yplus(j+1) d],[xplus yplus(j) d]);
             patr.pan = vertcat(patr.pan,[l+1 l+2 l+3 l+4]);
        end
    end

    if b<eps(b) && bplus>eps(bplus) %triangle starting at 0 breadth
        l = length(patr.crd);
        patr.crd = vertcat(patr.crd,[x 0 d],[x+0.5*pansize 0 d],[xplus 0 d],[xplus yplus(end) d]);
        patr.pan = vertcat(patr.pan,[l+1 l+2 l+3 l+4]);
    end
    if b>eps(b) && bplus<eps(bplus) %triangle ending at 0 breadth
        l = length(patr.crd);
        patr.crd = vertcat(patr.crd,[x 0 d],[x+0.5*pansize 0 d],[xplus 0 d],[x y(end) d]);
        patr.pan = vertcat(patr.pan,[l+1 l+2 l+3 l+4]);
    end
end

%% Bow 
div = ceil((cross{1,1}(3,1)-cross{1,1}(3,end))/(pansize/2))+1;
spacing = (cross{1,1}(3,1)-cross{1,1}(3,end))/div;
x = cross{1,1}(1,1);
if div>1        
    for i=1:div %for each row
        %determine z and y coordinates (z scalar, y vector)
        z = cross{1,1}(3,end) + spacing*(i-1);
        zplus = cross{1,1}(3,end) + spacing*(i);
        b = interp1(cross{1,1}(3,:),cross{1,1}(2,:),z); %breadth
        bplus = interp1(cross{1,1}(3,:),cross{1,1}(2,:),zplus);
        if isnan(b);b=0;end
        if isnan(bplus);bplus=0;end
        p = ceil(b/(pansize/2))+1;
        y = linspace(0,b,p);
        yplus = linspace(0,bplus,p);

        if b>eps(b) && bplus>eps(bplus) %square panels can be drawn
            for j=1:p-1
                 l = length(patr.crd);
                 patr.crd = vertcat(patr.crd,[x y(j) z],[x y(j+1) z],[x yplus(j+1) zplus],[x yplus(j) zplus]);
                 patr.pan = vertcat(patr.pan,[l+4 l+3 l+2 l+1]);
            end
        end

        if b<eps(b) && bplus>eps(bplus) %triangle starting at 0 breadth
            l = length(patr.crd);
            patr.crd = vertcat(patr.crd,[x 0 z],[x 0 z+0.5*(zplus-z)],[x 0 zplus],[x yplus(end) zplus]);
            patr.pan = vertcat(patr.pan,[l+4 l+3 l+2 l+1]);
        end
        if b>eps(b) && bplus<eps(bplus) %triangle ending at 0 breadth
            l = length(patr.crd);
            patr.crd = vertcat(patr.crd,[x 0 z],[x 0 z+0.5*(zplus-z)],[x 0 zplus],[x y(end) z]);
            patr.pan = vertcat(patr.pan,[l+4 l+3 l+2 l+1]);
        end
    end
end

%% Stern
div = ceil((cross{1,end}(3,1)-cross{1,end}(3,end))/(pansize/2))+1;
spacing = (cross{1,end}(3,1)-cross{1,end}(3,end))/div;
x = cross{1,end}(1,1);
if div>1        
    for i=1:div %for each row
        %determine z and y coordinates (z scalar, y vector)
        z = cross{1,end}(3,end) + spacing*(i-1);
        zplus = cross{1,end}(3,end) + spacing*(i);
        b = interp1(cross{1,end}(3,:),cross{1,end}(2,:),z); %breadth
        bplus = interp1(cross{1,end}(3,:),cross{1,end}(2,:),zplus);
        if isnan(b);b=0;end
        if isnan(bplus);bplus=0;end
        p = ceil(b/(pansize/2))+1;
        y = linspace(0,b,p);
        yplus = linspace(0,bplus,p);

        if b>eps(b) && bplus>eps(bplus) %square panels can be drawn
            for j=1:p-1
                 l = length(patr.crd);
                 patr.crd = vertcat(patr.crd,[x y(j) z],[x y(j+1) z],[x yplus(j+1) zplus],[x yplus(j) zplus]);
                 patr.pan = vertcat(patr.pan,[l+1 l+2 l+3 l+4]);
            end
        end

        if b<eps(b) && bplus>eps(bplus) %triangle starting at 0 breadth
            l = length(patr.crd);
            patr.crd = vertcat(patr.crd,[x 0 z],[x 0 z+0.5*(zplus-z)],[x 0 zplus],[x yplus(end) zplus]);
            patr.pan = vertcat(patr.pan,[l+1 l+2 l+3 l+4]);
        end
        if b>eps(b) && bplus<eps(bplus) %triangle ending at 0 breadth
            l = length(patr.crd);
            patr.crd = vertcat(patr.crd,[x 0 z],[x 0 z+0.5*(zplus-z)],[x 0 zplus],[x y(end) z]);
            patr.pan = vertcat(patr.pan,[l+1 l+2 l+3 l+4]);
        end
    end
end

%finish structure
patr.nbody = 1;
patr.bpan = length(patr.pan);
patr.npan = length(patr.pan);
patr.ncrd = length(patr.crd);
patr.info = [date ' RHDHV ' plan.info];

% Remove broken panel check
for p=1:length(patr.pan)
    vertices = patr.pan(p,:);
    v1 = patr.crd(vertices(1),:) - patr.crd(vertices(2),:);
    base = mag(v1);
    
    v2 = patr.crd(vertices(2),:) - patr.crd(vertices(3),:);
    proj = dot(v2,v1/mag(v1));
    
    height = sqrt(mag(v2)^2 - proj^2);
    
    onehalf = 0.5*base*height;
    
    v3 = patr.crd(vertices(3),:) - patr.crd(vertices(4),:);
    base = mag(v3);
    
    v4 = patr.crd(vertices(4),:) - patr.crd(vertices(1),:);
    proj = dot(v4,v3/mag(v3));
    
    height = sqrt(mag(v4)^2 - proj^2);
    
    otherhalf = 0.5*base*height;
    area(p) = onehalf + otherhalf;
end
patr.pan(isnan(area'),:) = [];

patr.bpan = length(patr.pan);
patr.npan = length(patr.pan);
patr.ncrd = length(patr.crd);

assignin('base','info',patr.info);
assignin('base','patr',patr);

oldpath = pwd;
path = [pwd '\PATRAN\'];
cd(path)
patran_write(patr,path,'hull.pat');
patran_mirror('hull.pat',path);
patran_combine('hull.pat','hull_mir.pat',path)
patran_offset('combined.pat',path,[-plan.L/2 0 -d]')
cd(oldpath)