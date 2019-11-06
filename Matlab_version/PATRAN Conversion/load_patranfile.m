function [ ship ] = load_patranfile( file, sf )

    % used for scaling
    length_scalefactor = sf.L; % *100%
    breadth_scalefactor = sf.B; % *100%
    depth_scalefactor = sf.D; % *100%

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

    %% Data
    reduc = sf.reduc;
    V = pat.crd; %short notation for vertices

    xmin = min(V(:,1));
    xmax = max(V(:,1));
    ship.L = (xmax-xmin);
    ship.B = max(V(:,2)) - min(V(:,2));
    ship.D = max(-V(:,3));
    ship.T = max(-V(:,3))-reduc;

    wls = 15;
    wl = ship.T/wls;

    %% Get projected length of waterlines
    xkeel = linspace(xmin,xmax,1000);
    keelline = griddata(V(:,1),V(:,2),V(:,3),xkeel,0,'linear');

    
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
    y = linspace(0,ship.T,wls+1);
    A_w = fliplr([A_w 0]);

    %% Compute eigenfrequencies and minimal panel size
    ship.displacement = simps(y,A_w);

end

