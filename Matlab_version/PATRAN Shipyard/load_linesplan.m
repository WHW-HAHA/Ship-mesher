function [ linesplan ] = load_linesplan(filepath)
    %   LOAD_LINESPLAN  Loads a TU Delft linesplan from file
	%   
    %   Usage: load_linesplan(filepath)
	%   The function returns a structure containing the following fields;
	%   coordinates: a 3x#ords cell containing y and z coordinates per x
	%   L: the length of the ship
	%   B: the breath of the ship
    %   T: the draft of the ship
    %   filepath: the path to the file
    %   type: the ship type as specified in the linesplan
    
    try
        fid = fopen(filepath);
        fgetl(fid); %skip first line
    catch
        error = 'The file could not be opened.'
    end
    name = fgetl(fid); %second line contains ship details
    fileinfo = textscan(name,'%s %[^,] %s %n %s %n %s %n %[^)] %[^.] $c');

    line = fgetl(fid); %third line contains scaling information
    scaleinfo = textscan(line,'%n %n %n %n');

    line = fgetl(fid); %forth line are the number of ordinates
    ords = textscan(line,'%n');
    ords = ords{1};

    spacing = [];
    for l=1:ceil(ords/6) %gets the ordinate spacing information
        line = fgetl(fid);
        n = textscan(line,'%n');
        spacing = vertcat(spacing,n{1});
    end
    x = [0 ; cumsum(spacing)];

    fgetl(fid); %skip a line

    coordinates = cell(2,ords+1);
    for o=1:ords+1 %gets y and z coordinates per ordinate
        clear z y
        line = fgetl(fid);
        firstline = textscan(line,'%n %n %n');
        coords = firstline{2} + 1;

        zy = [];
        for i=1:ceil((coords*2)/6)
            line = fgetl(fid);
            n = textscan(line,'%n');
            zy = vertcat(zy,n{1});
        end
        c=1;
        for i=1:2:length(zy)-1
            z(c) = zy(i);
            y(c) = zy(i+1);
            c = c + 1;
        end

        coordinates{1,o} = z;
        coordinates{2,o} = y;
    end

    line = fgetl(fid); %last line lists the main particulars
    main_particulars = textscan(line,'%n %n %n');
    L = main_particulars{1};
    B = main_particulars{2};
    T = main_particulars{3};
    d = textscan(fileinfo{9}{:},'%c %n');
    D = d{2};

    x = x*L; %apply scaling and store coordinates in structure
    for o=1:ords+1 
       misspacing = scaleinfo{1} - coordinates{2,o}(end);
       linesplan.coordinates{2,o} = coordinates{1,o}*B/scaleinfo{3};
       linesplan.coordinates{3,o} = (coordinates{2,o}+misspacing)*T;
       linesplan.coordinates{1,o} = x(o);
    end
    
    linesplan.L = L; %attach more information to the structure
    linesplan.B = B;
    linesplan.T = T;
    linesplan.D = D;
    linesplan.filepath = filepath;
    linesplan.type = fileinfo{2}{:};
    linesplan.ords = o;
    linesplan.info = name;
    
    fclose(fid);
    
