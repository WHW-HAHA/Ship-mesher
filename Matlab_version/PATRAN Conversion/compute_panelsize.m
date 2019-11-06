function [ info ] = compute_panelsize(file,sf,k_xx,k_yy,KG,waterdepth)
   
    length_scalefactor = sf.L; % *100%
    breadth_scalefactor = sf.B; % *100%
    depth_scalefactor = sf.D; % *100%

    g = 9.81;

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
    if nwl==0;disp(['No waterline elementes found in ' path file '.']);end
    if exist('panwl','var')
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
    clear b uc nuc nwl panwl 

    %% Data
    reduc = sf.reduc;
    V = pat.crd; %short notation for vertices

    xmin = min(V(:,1));
    xmax = max(V(:,1));
    L = (xmax-xmin);

    wls = 15;
    newdraft = max(-V(:,3))-reduc;
    wl = newdraft/wls;

    %% Get projected length of waterlines
    xkeel = linspace(xmin,xmax,1000);
    keelline = griddata(V(:,1),V(:,2),V(:,3),xkeel,0,'linear');

    figure;
    plot(xkeel,keelline);
    axis equal
    title('Interpolated keelline')
   
    %because direct interpolation did not work, the keelline is traversed until
    %the next waterline depth is reached. The minimum and maximum x coordinate
    %are stored so that the hull can be interpolated between those lengths
    %later
    for i=1:wls  
        j=1;while keelline(j)>(-reduc-(i-1)*wl) && j<length(keelline)
            j=j+1;
        end
        xmi(i) = xkeel(j);    
        j=length(xkeel);while keelline(j)<(-reduc-(i-1)*wl) && j>1
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
    wavelength = dispersion_relation(max(om)*2,waterdepth);
    pansize = wavelength/4;
    info.needpsize = round(pansize,1);
    info.Aw = A_w(end);
    info.GM_l = GM_l;
    info.GM_t = GM;
    info.displacement = displacement;
    info.I_t = I_t(end);
    info.I_l = I_l;   
    info.wn_3 = omega_heave;
    info.wn_4 = omega_roll;
    info.wn_5 = omega_pitch;

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
    
    info.nowpsize = round(psize,1);

end

