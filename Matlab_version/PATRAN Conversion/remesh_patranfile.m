function [ output_args ] = remesh_patranfile( file, sf,split )
   
   
    %% Read PATRAN file and apply scalefacors
    pat = patran_read_pat(file);
    pat.crd(:,1) = pat.crd(:,1)*sf.L;
    pat.crd(:,2) = pat.crd(:,2)*sf.B;
    pat.crd(:,3) = pat.crd(:,3)*sf.D;
    reduc = sf.reduc;
    
    %% Split panels if needed
    if split==1
        %% Make the new panels
        pan = []; %new panel list
        for p=1:length(pat.pan)
            %to divide one panel into 4 matlab is going to find the midpoints on
            %each face of the panel and add an extra point in the middle. Between
            %these points and the 4 existing points (a total of 9) new panels are
            %defined.

            %get all corner points
            vertices = pat.pan(p,:);

            %midpoint from 1 to 2
            v1 = pat.crd(vertices(2),:) - pat.crd(vertices(1),:);
            mid1 = pat.crd(vertices(1),:) + 0.5*v1;

            %midpoint from 2 to 3
            v2 = pat.crd(vertices(3),:) - pat.crd(vertices(2),:);
            mid2 = pat.crd(vertices(2),:) + 0.5*v2;

            %midpoint from 3 to 4
            v3 = pat.crd(vertices(4),:) - pat.crd(vertices(3),:);
            mid3 = pat.crd(vertices(3),:) + 0.5*v3;

            %midpoint from 4 to 1
            v4 = pat.crd(vertices(1),:) - pat.crd(vertices(4),:);
            mid4 = pat.crd(vertices(4),:) + 0.5*v4;

            %mid of panel, taken as the halfway point between two opposing corners
            v5 = mid3 - mid1;
            midmid = mid1 + 0.5*v5; 

            %there are now 4 panels instead of the old one. The old one is removed
            %and the new ones are added
            % - vertices1 -> mid1 -> midmid -> mid4
            % - mid1 -> vertices2 -> mid2 -> midmid
            % - mid4 -> midmid -> mid3 -> vertices4
            % - midmid -> mid2 -> vertices3 -> mid3
            pan1 = [vertices(1) length(pat.crd)+1 length(pat.crd)+5 length(pat.crd)+4];
            pan2 = [length(pat.crd)+1 vertices(2) length(pat.crd)+2 length(pat.crd)+5];
            pan3 = [length(pat.crd)+4 length(pat.crd)+5 length(pat.crd)+3 vertices(4)];
            pan4 = [length(pat.crd)+5 length(pat.crd)+2 vertices(3) length(pat.crd)+3];

            %add panels to the panel list
            pan = vertcat(pan,pan1,pan2,pan3,pan4);

            %add the newly found midpoints to the coordinate list
            pat.crd = vertcat(pat.crd,mid1,mid2,mid3,mid4,midmid);

        end

        pat.pan = pan; %overwrite old panel list


        %% Remove unused coordinates
        %now there are a lot of double points to be found in the coordinate list
        %because panels share faces, and each neighbour made it's own midpoints.
        %Here we will remove those double entries.
        pat.crd = int64(ceil(pat.crd.*1e8)); %convert to integers to avoid floating point errors

        c=1;while c<length(pat.crd)
            index = sort(find(all(bsxfun(@eq, pat.crd, pat.crd(c,:)), 2))); %find indexes of duplicate coordinates
            original = index(1); %store the index of the original coordinate
            for i=2:length(index) %for each duplicate, 
                pat.pan(pat.pan==index(i)) = original; %overwrite references to the copy
                pat.crd(index(i),:) = []; %remove the copy itself
                pat.pan(pat.pan>index(i)) = pat.pan(pat.pan>index(i))-1; %shift all the references to coordinates down
                index = index - 1;
            end
            c=c+1;
        end

        %% Write PATRAN file
        pat.crd = double(pat.crd)/1e8; %convert the integers back to doubles
        pat.npan = length(pat.pan);
        pat.ncrd = length(pat.crd);
        pat.bpan = pat.npan;
    end
    %% Remove waterline panels
    % CODE FROM MATPAT BY MARIN
    % Find waterline panels
    nwl=0;b=1;wl(1)=0;
    for i=1:pat.npan
        % Check for waterline panels
        if sum(abs(pat.crd(pat.pan(i,:),3)))==0
            nwl=nwl+1;
            panwl(nwl)=i;
        end
        % Check number of bodypanels left
        if i==sum(pat.bpan(1:b))
            if nwl~=0
                wl(b+1)=length(panwl)-sum(wl(1:b));
                bpan(b)=pat.bpan(b)-wl(b+1);
            else
                disp(['No waterline elementes found in body ' num2str(b) ' in ' path file '.']);
            end
            b=b+1;
        end
    end
    if exist('panwl','var')
        if nwl==0;disp(['No waterline elementes found in ' path file '.']);return;end
        pat.pan(panwl,:)=[];
        pat.bpan=bpan;
        pat.npan=size(pat.pan,1);
        % Find unused coordinates
        nuc=0;
        for i=1:pat.ncrd
            if isempty(find(pat.pan==i,1))
                nuc=nuc+1;
                uc(nuc)=i;
            end
        end

        % Delete unused coordinates and order panels
        if exist('uc','var')
            uc=sort(uc,'descend');
            for i=1:length(uc)
                pat.pan(find(pat.pan>uc(i)))=pat.pan(find(pat.pan>uc(i)))-1;
            end
            pat.crd(uc,:)=[];
            pat.ncrd=size(pat.crd,1);
        end
    end
    % END OF CODE FROM MATPAT BY MARIN
    
    %% Waterline
    V = pat.crd;
    V(V(:,2)<0,:) = [];
    xpoints = linspace(min(V(:,1)),max(V(:,1)),100);
    waterline{1}(:,2) = griddata(V(:,1),V(:,3),V(:,2),xpoints,-reduc);
    waterline{1}(:,1) = xpoints;
    waterline{1}(isnan(waterline{1}(:,2)),2) = 0; 

    %% Put the PATRAN file at its new draught
    delete = zeros(size(pat.pan));
    for p=1:length(pat.pan)
        %get all corner points
        vertices = pat.pan(p,:);
        c = [];
        %check how many of them are above the new draft
        for v=1:4
           if pat.crd(vertices(v),3)>-(reduc) 
               c = [c v];
           end
        end
        if length(c)==4
            %if four points are above the new draft the entire panel can simply
            %be removed
            delete(p) = 1;
        end
        if length(c)==1
            %if 3 points are below the new draft it is important that the one 
            %above the draft gets put on the new waterline. After this the
            %point closest to the waterline is chosen and also put on the
            %waterline.
            pat.crd(vertices(c(1)),3) = -(reduc);
            pat.crd(vertices(c(1)),2) = sign(pat.crd(vertices(c(1)),2))*interp1(waterline{1}(:,1),waterline{1}(:,2),pat.crd(vertices(c(1)),1));
            b=1:4;b=b(b~=c(1));
            [mi,i] = min(abs(-(reduc)-pat.crd(vertices(b),3)));
            pat.crd(vertices(b(i)),3) = -(reduc);
            pat.crd(vertices(b(i)),2) = sign(pat.crd(vertices(b(i)),2))*interp1(waterline{1}(:,1),waterline{1}(:,2),pat.crd(vertices(b(i)),1));
        end
        if length(c)==2
           %when 2 points are above the new draft it is first determined if the
           %bottom or the top part of the panel is closer to the waterline. The
           %closer one of those ends is then put on the waterline.
            [mi,i] = min(abs(-(reduc)-pat.crd(vertices(1:4),3)));

            if numel(c(c==i))
                pat.crd(vertices(c(1)),3) = -(reduc);
                pat.crd(vertices(c(1)),2) = sign(pat.crd(vertices(c(1)),2))*interp1(waterline{1}(:,1),waterline{1}(:,2),pat.crd(vertices(c(1)),1));
                pat.crd(vertices(c(2)),3) = -(reduc);
                pat.crd(vertices(c(2)),2) = sign(pat.crd(vertices(c(2)),2))*interp1(waterline{1}(:,1),waterline{1}(:,2),pat.crd(vertices(c(2)),1));
            else
                b=1:4;b=b(b~=c(1)&b~=c(2));
                pat.crd(vertices(b(1)),3) = -(reduc);
                pat.crd(vertices(b(1)),2) = sign(pat.crd(vertices(b(1)),2))*interp1(waterline{1}(:,1),waterline{1}(:,2),pat.crd(vertices(b(1)),1));
                pat.crd(vertices(b(2)),3) = -(reduc); 
                pat.crd(vertices(b(2)),2) = sign(pat.crd(vertices(b(2)),2))*interp1(waterline{1}(:,1),waterline{1}(:,2),pat.crd(vertices(b(2)),1));
            end
        end
    end

    %remove the panels that are now above the draught because they have been
    %altered
    for p=1:length(pat.pan)
        %get all corner points
        vertices = pat.pan(p,:);
        c = [];
        %check how many of them are above the new draft
        for v=1:4
           if pat.crd(vertices(v),3)>=-(reduc) 
               c = [c v];
           end
        end
        if length(c)==4
            %if four points are above the new draft the entire panel can simply
            %be removed
            delete(p) = 1;
        end
    end
    
    pat.crd(:,3) = pat.crd(:,3) + reduc; %make sure the new waterline lies at y=0 again
    pat.pan(delete==1,:) = []; %delete the panels that need to be deleted
    
    %% Create Waterplane
    %the waterplane is divided into a number (div) of rows for which the
    %breadth is determined. The breadth is divided in a number of panels that
    %are close to square.
    div = 60;
    xmin =  waterline{1}(1,1);xmax = max(pat.crd(:,1));
    spacing = (xmax - xmin) / div;
    for i=1:div %for each row
        %determine x and y coordinates (x scalar, y vector)
        x = xmin+(i-1)*spacing;
        xplus = xmin+(i)*spacing;
        b = interp1(waterline{1}(:,1),waterline{1}(:,2),x); %breadth
        bplus = interp1(waterline{1}(:,1),waterline{1}(:,2),xplus);
        if isnan(b);b=0;end
        if isnan(bplus);bplus=0;end
        p = ceil(b/spacing)+1;
        y = linspace(0,b,p);
        yplus = linspace(0,bplus,p);

        if b>eps(b) && bplus>eps(bplus) %square panels can be drawn
            for j=1:p-1
                 l = length(pat.crd);
                 pat.crd = vertcat(pat.crd,[x y(j) 0],[x y(j+1) 0],[xplus yplus(j+1) 0],[xplus yplus(j) 0]);
                 pat.pan = vertcat(pat.pan,[l+1 l+2 l+3 l+4]);
                 % mirrored
                 pat.crd = vertcat(pat.crd,[x -y(j) 0],[xplus -yplus(j) 0],[x -y(j+1) 0],[xplus -yplus(j+1) 0]);
                 pat.pan = vertcat(pat.pan,[l+5 l+6 l+8 l+7]);
            end
        end

        if b<eps(b) && bplus>eps(bplus) %triangle starting at 0 breadth
            l = length(pat.crd);
            pat.crd = vertcat(pat.crd,[x 0 0],[x+0.5*spacing 0 0],[xplus 0 0],[xplus yplus(end) 0]);
            pat.pan = vertcat(pat.pan,[l+1 l+2 l+3 l+4]);
            % mirrored
            pat.crd = vertcat(pat.crd,[x 0 0],[x+0.5*spacing 0 0],[xplus 0 0],[xplus -yplus(end) 0]);
            pat.pan = vertcat(pat.pan,[l+5 l+6 l+7 l+8]);
        end
        if b>eps(b) && bplus<eps(bplus) %triangle ending at 0 breadth
            l = length(pat.crd);
            pat.crd = vertcat(pat.crd,[x 0 0],[x+0.5*spacing 0 0],[xplus 0 0],[x y(end) 0]);
            pat.pan = vertcat(pat.pan,[l+1 l+2 l+3 l+4]);
            % mirrored
            pat.crd = vertcat(pat.crd,[x 0 0],[x+0.5*spacing 0 0],[xplus 0 0],[x -y(end) 0]);
            pat.pan = vertcat(pat.pan,[l+5 l+6 l+7 l+8]);
        end
    end
    
    
    pat.npan = length(pat.pan);
    pat.bpan = length(pat.pan);

    figure
    patran_open(pat)
    [pathstr, name] = fileparts(file); 
    patran_rem_unusedcrd(pat,[pwd '\'],[name '_scaled.pat']) %write and plot PATRAN file


end

