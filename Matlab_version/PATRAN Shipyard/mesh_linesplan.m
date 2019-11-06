function mesh_linesplan(plan,pansize,pandist)
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
% figure()
% hold on
% for i=1:length(cross)
%    plot3(cross{i}(1,:),cross{i}(2,:),cross{i}(3,:),'d','Color',[64 104 225]/255,'MarkerSize',2)
% end
% view(3)
% axis equal
% hold off
%% Bow and stern
%stern
for i=1:size(coords2{1},2)
    p = ceil(coords2{1}(2,i)/pansize)+1;
    vec = [repmat(coords2{1}(1,i),1,p); linspace(coords2{1}(2,i),0,p); repmat(coords2{1}(3,i),1,p)];
    stern{i} = vec;
end

for i=1:length(stern)-1
    extra = abs(size(stern{i+1},2)-size(stern{i},2));
    if extra>0 
        dz = (coords2{1}(3,i+1)-coords2{1}(3,i))/(extra+1);
        for e=1:extra
            stern{i} = horzcat(stern{i},stern{i}(:,end)+[0 0 dz]');
        end
    end
    %plot3(stern{i}(1,:),stern{i}(2,:),stern{i}(3,:),'.r-','MarkerSize',5)
end
%plot3(stern{i+1}(1,:),stern{i+1}(2,:),stern{i+1}(3,:),'.r-','MarkerSize',5)

for i=2:size(cross{end},2)-1
    if cross{end}(2,i)<0.05
        cross{end}(2,i)=0.05;
    end 
end

%bow
for i=1:size(coords2{end},2)
    if coords2{end}(2,i)<=0.05
        coords2{end}(2,i)=0.05;
    end
    p = ceil(coords2{end}(2,i)/pansize)+1;
    vec = [repmat(coords2{end}(1,i),1,p); linspace(coords2{end}(2,i),0,p); repmat(coords2{end}(3,i),1,p)];
    bow{i} = vec;
end

for i=1:length(bow)-1
    extra = size(bow{i+1},2)-size(bow{i},2);
    if extra>0 
        dz = (coords2{end}(3,i+1)-coords2{end}(3,i))/(extra+1);
        for e=1:extra
            bow{i} = horzcat(bow{i},bow{i}(:,end)+[0 0 dz]');
        end
    end
    if extra<0 
        dz = -(coords2{end}(3,i+1)-coords2{end}(3,i))/(abs(extra)+1);
        for e=1:abs(extra)
            bow{i+1} = horzcat(bow{i+1},bow{i+1}(:,end)+[0 0 dz]');
        end
    end
    %plot3(bow{i}(1,:),bow{i}(2,:),bow{i}(3,:),'.k-','MarkerSize',5)
end
%plot3(bow{i+1}(1,:),bow{i+1}(2,:),bow{i+1}(3,:),'.k-','MarkerSize',5)

%% Lid
for i=1:length(coords{1})
    p = ceil(coords{1}(2,i)/pansize)+1;
    vec = [repmat(coords{1}(1,i),1,p); linspace(coords{1}(2,i),0,p); repmat(coords{1}(3,i),1,p)];
    lid{i} = vec;
end

lid2 = lid;
for i=1:length(lid)-1
    extra = abs(size(lid2{i+1},2)-size(lid2{i},2));
    if extra>0 
        dx = (coords{1}(1,i+1)-coords{1}(1,i))/(extra+1);
        for e=1:extra
            if i<round(length(lid)/2)
                lid{i} = horzcat(lid{i},lid{i}(:,end)+[dx 0 0]');
            elseif i>round(length(lid)/2)
                lid{i+1} = horzcat(lid{i+1},lid{i+1}(:,end)+[-dx 0 0]');
            end
        end
    end
    plot3(lid{i}(1,:),lid{i}(2,:),lid{i}(3,:),'.b-','MarkerSize',5)
end
plot3(lid{end}(1,:),lid{end}(2,:),lid{end}(3,:),'.g-','MarkerSize',5)
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

% add stern
for i=1:length(stern)
    scnums{i} = [cnums{1}(size(coords2{1},2)+1-i)]; 
    for c=2:size(stern{i},2)
        patr.crd = vertcat(patr.crd,stern{i}(:,c)');
        scnums{i} = vertcat(scnums{i},panelnumber); 
        panelnumber = panelnumber+1;
    end
end

for c=1:length(stern)-1
    for i=1:length(scnums{c})-1
        patr.pan = vertcat(patr.pan,[scnums{c}(i+1) scnums{c+1}(i+1) scnums{c+1}(i) scnums{c}(i)]);
    end
end

% add bow
for i=1:length(bow)
    bcnums{i} = [cnums{end}(size(coords2{end},2)+1-i)]; 
    for c=2:size(bow{i},2)
        patr.crd = vertcat(patr.crd,bow{i}(:,c)');
        bcnums{i} = vertcat(bcnums{i},panelnumber); 
        panelnumber = panelnumber+1;
    end
end

for c=1:length(bow)-1
    for i=1:length(bcnums{c})-1
        patr.pan = vertcat(patr.pan,[bcnums{c+1}(i) bcnums{c+1}(i+1) bcnums{c}(i+1) bcnums{c}(i)]);
    end
end

% add lid
%lcnums{1} = scnums{end}(1:ceil(coords{1}(2,1)/pansize)+1);
for i=1:length(lid)
    lcnums{i} = [cnums{i}(1)]; 
    for c=2:size(lid{i},2)
        patr.crd = vertcat(patr.crd,lid{i}(:,c)');
        lcnums{i} = vertcat(lcnums{i},panelnumber); 
        panelnumber = panelnumber+1;
    end
end
%lcnums{i+1} = bcnums{end}(1:ceil(coords{1}(2,end)/pansize)+1);

for c=1:length(lid)
    for i=1:length(lcnums{c})-1
        if c<round(length(lid)/2)
            patr.pan = vertcat(patr.pan,[lcnums{c+1}(i) lcnums{c+1}(i+1) lcnums{c}(i+1) lcnums{c}(i)]);
        elseif c>round(length(lid)/2)
            patr.pan = vertcat(patr.pan,[lcnums{c}(i+1) lcnums{c-1}(i+1) lcnums{c-1}(i) lcnums{c}(i)]);
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
patr.npan = patr.npan - length(patr.pan(isnan(area'),:));
patr.bpan = patr.bpan - length(patr.pan(isnan(area'),:));
patr.pan(isnan(area'),:) = [];


assignin('base','info',patr.info);
assignin('base','patr',patr);

oldpath = pwd;
path = [pwd '\PATRAN\'];
cd(path)
d= coords{1,1}(3);
patran_write(patr,path,'hull.pat');
patran_mirror('hull.pat',path);
patran_combine('hull.pat','hull_mir.pat',path)
patran_offset('combined.pat',path,[-plan.L/2 0 -d]')
cd(oldpath)