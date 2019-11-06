function [ output_args ] = scale_patranfile(file,sf)


    kL = 0.25; %k_xx = kL * Lpp
    kB = 0.30; %k_yy = kB * B
    KGT = 0.9; %KG = KGT * T
    
    % used for scaling
    length_scalefactor = sf.L; % *100%
    breadth_scalefactor = sf.B; % *100%
    depth_scalefactor = sf.D; % *100%
    draught_reduction = sf.To - sf.T; % in units, after scaling

    %% Read PATRAN file and apply scalefacors
    pat = patran_read_pat(file);
    pat.crd(:,1) = pat.crd(:,1)*length_scalefactor;
    pat.crd(:,2) = pat.crd(:,2)*breadth_scalefactor;
    pat.crd(:,3) = pat.crd(:,3)*depth_scalefactor;

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
    % END OF CODE FROM MATPAT BY MARIN

    %% Data
    V = pat.crd; %short notation for vertices

    xmin = min(V(:,1));
    xmax = max(V(:,1));
    
    % used for estimating required panel size
    k_xx = kL*(xmax-xmin); %radius of gyration about the x-axis
    k_yy = kB*(max(V(:,2)) - min(V(:,2))); %radius of gyration about the y-axis
    KG = KGT*max(-V(:,3)); %distance in height from baseline to CoG
    g = 9.81;
    
    wls = 15;
    newdraft = max(-V(:,3))-draught_reduction;
    wl = newdraft/wls;

    %% Get projected length of waterlines
    xkeel = linspace(xmin,xmax,1000);
    keelline = griddata(V(:,1),V(:,2),V(:,3),xkeel,0,'linear');

    figure;
    plot(xkeel,keelline);
    axis equal
    title('Interpolated keelline')

    reduc = draught_reduction;
    %because direct interpolation did not work, the keelline is traversed until
    %the next waterline depth is reached. The minimum and maximum x coordinate
    %are stored so that the hull can be interpolated between those lengths
    %later
    for i=1:wls  
        j=1;while keelline(j)>(-reduc-(i-1)*wl)
            j=j+1;
            end
        xmi(i) = xkeel(j);    
        j=length(xkeel);while keelline(j)<(-reduc-(i-1)*wl)
            j=j-1;
            end
        xma(i) = xkeel(j);   
    end
    xmi(1) = xmin;
    xma(1) = xmax;

    %% Compute waterline areas and area moments
    V(V(:,2)<0,:) = []; %remove other half of the ship to make interpolation possible
    waterline = cell(wls,2);
    A_w = zeros(1,wls);
    I_t = zeros(1,wls);
    for i=1:wls
        %interpolate waterlines
        xpoints = linspace(xmi(i),xma(i),100);
        waterline{i}(:,2) = griddata(V(:,1),V(:,3),V(:,2),xpoints,-reduc-(i-1)*wl);
        waterline{i}(:,1) = xpoints;
        waterline{i}(isnan(waterline{i}(:,2)),2) = 0; 

        % Integrate waterline area and moment
        A_w(i) = 2 * simps(waterline{i}(:,1),waterline{i}(:,2));
        I_t(i) = 2/3 * simps(waterline{i}(:,1),waterline{i}(:,2).^3);
    end
    y = linspace(0,newdraft,wls+1);
    A_w = fliplr([A_w 0]);

    figure
    plot(waterline{1}(:,1),waterline{1}(:,2))
    title('Interpolated waterline')
    axis equal

    %% Compute eigenfrequencies and minimal panel size
    displacement = simps(y,A_w);
    BM = I_t(1)/displacement;
    KB = simps(y,A_w.*y)/displacement;
    GM = KB + BM - KG;

    omega_heave = sqrt(A_w(end)*g/displacement);
    omega_roll = sqrt(GM*g/(k_xx^2));

    I_l = 2*simps(waterline{1}(:,1),waterline{1}(:,2).*waterline{1}(:,1).^2); 
    BM_l = I_l/displacement;
    GM_l = KB + BM_l - KG;

    omega_pitch = sqrt(GM_l*g/(k_yy^2));

    om = [omega_pitch omega_roll omega_heave];
    wavelength = ((2*pi*g)/(max(om)*2)^2);

    pansize = wavelength/4;
    recommended_panel_size = round(pansize,1);

    %% Establish panel size in PATRAN file
    % taken as the median of half of the faces of all panels
    mags = zeros(length(pat.pan)*2,1);
    for p=1:length(pat.pan)
        vertices = pat.pan(p,:);
        v1 = pat.crd(vertices(1),:) - pat.crd(vertices(2),:);
        v2 = pat.crd(vertices(2),:) - pat.crd(vertices(3),:);
        mags(2*p) = mag(v1);
        mags(2*p-1) = mag(v2);
    end

    psize = median(round(mags,2,'significant'));
    figure
    hist(mags)
    title('Panel sizes in PATRAN file')

    disp(['The approximate panel size in the PATRAN file is ' num2str(psize)]);
    disp(['The minimum calculated panel size is ' num2str(recommended_panel_size)]);

    clear bpan g i j mags om p v1 v2 wavelength wl wls xma xmi xmax xmin

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
    xmin = min(pat.crd(:,1));xmax = max(pat.crd(:,1));
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
    patran_rem_unusedcrd(pat,[pwd '\'],outputfilename) %write and plot PATRAN file
end

