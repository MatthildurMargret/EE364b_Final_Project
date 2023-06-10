
%Upload point clouds from txt files
d=dir('*.txt');         % all .txt files in working directory
N=length(d);            
sizeA = [4 938];
formatSpec = '%f';  

%% Compute ground-truth volume over all frames

volume_track = zeros(27,1);

for i=1:27
    % Open each set of RV vertices
    fid = fopen(d(i).name,'r');
    A = fscanf(fid,formatSpec,sizeA);
    % Retrieve the relevant data and shift to include origin
    points = A(2:end,:)';
    x = points(:,1);
    y = points(:,2)+40;
    z = points(:,3)-100;
    % Get the corresponding connvectivity matrix
    allLines = readlines(d(i).name);
    allLines = allLines(939:2810);
    cells = num2cell(allLines);
    matrix = [cells{:}];
    T = zeros(1872,3);
    for j=1:1872
        c = regexp(matrix(j),' ', 'split');
        T(j,1) = str2double(c(4));
        T(j,2) = str2double(c(6));
        T(j,3) = str2double(c(8));
    end
    T = T+1;
    % Compute the volume and save 
    [v,a] = triangulationVolume(T,x,y,z);
    volume_track(i) = v;
    fclose(fid);
end

%% Sweeping horizontal sensor

% Sweep horizontal plane over z axis, from 20% to 80% of height
shifts = [0.2 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80] ;
errors = zeros(13,1);
best_error = 100;
best_height = 0;

for shift = 1:13
    
    z_shift = shifts(shift);

    % To measure the length of the sensor
    sensor_track = zeros(27,1);

    %Loop over each frame
    for i=1:27
        % Retrieve pointcloud and mesh corresponding to frame i
        fid = fopen(d(i).name,'r');
        A = fscanf(fid,formatSpec,sizeA);
        points = A(2:end,:)';
        allLines = readlines(d(i).name);
        allLines = allLines(939:2810);
        cells = num2cell(allLines);
        matrix = [cells{:}];
        T = zeros(1872,3);
        for j=1:1872
            c = regexp(matrix(j),' ', 'split');
            T(j,1) = str2double(c(4));
            T(j,2) = str2double(c(6));
            T(j,3) = str2double(c(8));
        end
        T = T+1;
    
        % Define the surface vertices and faces for the SurfaceIntersection
        % function
        Surface.vertices = points;
        Surface.faces = T;

        % A necessary rotation to match the cutting plane
        Surface.vertices = Surface.vertices * rotz(100);
    
        % Construct horizontal plane at given height
        xmax = max(Surface.vertices(:,1));
        ymax = max(Surface.vertices(:,2));
        zmax = max(Surface.vertices(:,3));
        xmin = min(Surface.vertices(:,1));
        ymin = min(Surface.vertices(:,2));
        zmin = min(Surface.vertices(:,3));
        x_mid = (xmax+xmin)/2;
        y_mid = (ymax+ymin)/2;

        height = abs(zmax-zmin);
        b = zmin + z_shift*height;

        % Define a large enough plane to intersect with entire surface
        radius_u = max(abs(xmax),abs(ymax));
        radius_l = max(abs(xmin),abs(ymin));
        radius = 3*max(radius_u,radius_l);
        [xpc,ypc] = pol2cart((0:2)'*pi/2,radius);
        z_coords = [b;b;b];
        Plane_sensor.vertices = [xpc ypc z_coords; xpc -ypc z_coords];
        Plane_sensor.faces = [1:3;4:6];

        % Rotate sensor plane by multiplying this with roty(degree)
        Plane_sensor.vertices = Plane_sensor.vertices; 

        % Compute intersection with RV
        [intersect, Surf] = SurfaceIntersection(Surface, Plane_sensor);

        % Restrict the intersection to feasible region
        temp = Surf.vertices;
        mask = temp(:,2) > 31.2208;
        sensor = Surf.vertices(mask,:);

        % Take the norm and save 
        sensor_track(i) = norm(sensor,2);
    
    end
    
    % Normalize volume to compare with strain
    volume_n = normalize((volume_track-min(volume_track))/max(volume_track));
    sensor_n = normalize((sensor_track-min(sensor_track))/max(sensor_track));
    
    % Compute the error for given coordinates
    errors(shift) = rmse(volume_n,sensor_n);

    % Keep track of best results
    if rmse(volume_n,sensor_n)^2 < best_error
        best_error = rmse(volume_n,sensor_n)^2;
        best_height = z_shift;
        best_strain = sensor_track;
    end

end

best_strain_horizontal = best_strain;

%% See the error curve 

plot(shifts*100,errors,'LineWidth',2)
xlabel("Shift from midpoint along z axis, %",'FontSize',14)
ylabel("RMS error",'FontSize',14)

%% Visualize the results 

plot(normalize(best_strain_horizontal),'LineWidth',2)
hold on
plot(normalize(volume_track),'LineWidth',2)
xlabel('Frame','FontSize',14)
title('Change in strain vs change in volume','FontSize',14)
legend('Sensor strain','Volume changes','FontSize',14)

%% Sweep over a1 rotations

% We rotate the sensor plane in one degree increments 
angle = -45:45;
errors_p = zeros(91,1);
best_error = 100;
best_angle = 0;

for rot = 1:91
    
    Rx = roty(angle(rot));

    % To measure the length of the sensor
    sensor_track = zeros(27,1);
    
    for i=1:27
        % Open each RV pointcloud and corresponding connectivity matrix
        fid = fopen(d(i).name,'r');
        A = fscanf(fid,formatSpec,sizeA);
        points = A(2:end,:)';
        allLines = readlines(d(i).name);
        allLines = allLines(939:2810);
        cells = num2cell(allLines);
        matrix = [cells{:}];
        T = zeros(1872,3);
        for j=1:1872
            c = regexp(matrix(j),' ', 'split');
            T(j,1) = str2double(c(4));
            T(j,2) = str2double(c(6));
            T(j,3) = str2double(c(8));
        end
        T = T+1;
        
        % Define Surface vertices and faces for intersection method
        Surface.vertices = points;
        Surface.vertices = Surface.vertices * rotz(100);
        Surface.faces = T;
    
        %Construct horizontal plane at middle
        xmax = max(Surface.vertices(:,1));
        ymax = max(Surface.vertices(:,2));
        xmin = min(Surface.vertices(:,1));
        ymin = min(Surface.vertices(:,2));
        x_mid = (xmax+xmin)/2;
        y_mid = (ymax+ymin)/2;
        
        zmax = max(Surface.vertices(:,3));
        zmin = min(Surface.vertices(:,3));
        z_mid = (zmax+zmin)/2;
        radius_u = max(abs(xmax),abs(ymax));
        radius_l = max(abs(xmin),abs(ymin));
        radius = 3*max(radius_u,radius_l);
        [xpc,ypc] = pol2cart((0:2)'*pi/2,radius);
        z_coords = [z_mid-20;z_mid-20;z_mid-20];
        Plane_sensor.vertices = [xpc ypc z_coords; xpc -ypc z_coords];
        Plane_sensor.faces = [1:3;4:6];

        %Rotate sensor plane
        Plane_sensor.vertices = Plane_sensor.vertices*Rx;

        %Compute intersection with RV
        [intersect, Surf] = SurfaceIntersection(Surface, Plane_sensor);

        %Cut the sensor
        mask = Surf.vertices(:,2) > 31.2208;
        sensor = Surf.vertices(mask,:);

        sensor_track(i) = norm(sensor,2);
    
    end
    
    %Normalize volume to compare with strain
    volume_n = (volume_track-min(volume_track))/max(volume_track);
    strain_n = (sensor_track-min(sensor_track))/max(sensor_track);
    
    %Compute the error for given coordinates
    errors_p(rot) = rmse(volume_n,strain_n)^2;
    if rmse(volume_n,sensor_n)^2 < best_error
        best_error = rmse(volume_n,sensor_n)^2;
        best_angle = angle(rot);
        best_strain = sensor_track;
    end

end

best_strain_rot = best_strain;

%% See the error curve 

plot(angle,errors_p,'LineWidth',2)
xlabel("Rotation around x-axis",'FontSize',14)
ylabel("RMS error",'FontSize',14)
title('Error over varying orientations of sensor plane','FontSize',14)

%% Visualize the results 

plot(normalize(best_strain),'LineWidth',2)
hold on
plot(normalize(volume_track),'LineWidth',2)
xlabel('Frame','FontSize',14)
title('Change in strain vs change in volume','FontSize',14)
legend('Sensor strain','Volume changes','FontSize',14)

%% Sweep over vertical planes

angle = -50:50;
errors = zeros(101,1);
best_error = 100;
best_angle_y = 0;

for rot = 1:101
    
    Ry = roty(angle(rot));

    %To measure the length of the entire intersection
    strain_track = zeros(27,1);

    %To measure the length of the sensor
    sensor_track = zeros(27,1);

    Rz_heart = rotz(100);
    
    for i=1:27
        fid = fopen(d(i).name,'r');
        A = fscanf(fid,formatSpec,sizeA);
        points = A(2:end,:)';
        allLines = readlines(d(i).name);
        allLines = allLines(939:2810);
        cells = num2cell(allLines);
        matrix = [cells{:}];
        T = zeros(1872,3);
        for j=1:1872
            c = regexp(matrix(j),' ', 'split');
            T(j,1) = str2double(c(4));
            T(j,2) = str2double(c(6));
            T(j,3) = str2double(c(8));
        end
        T = T+1;
    
        Surface.vertices = points;
        Surface.vertices = Surface.vertices * Rz_heart;
        Surface.faces = T;
    
        %Construct horizontal plane at middle
        xmax = max(Surface.vertices(:,1));
        ymax = max(Surface.vertices(:,2));
        zmax = max(Surface.vertices(:,3));
        xmin = min(Surface.vertices(:,1));
        ymin = min(Surface.vertices(:,2));
        zmin = min(Surface.vertices(:,3));
        x_mid = (xmax+xmin)/2;
        y_mid = (ymax+ymin)/2;
        z_mid = (zmax+zmin)/2;
        radius_u = max(abs(xmax),abs(ymax));
        radius_l = max(abs(xmin),abs(ymin));
        radius = 3*max(radius_u,radius_l);
        [xpc,ypc] = pol2cart((0:2)'*pi/2,radius);
        z_coords = [z_mid-10;z_mid-10;z_mid-10];
        Plane_custom.vertices = [xpc ypc z_coords; xpc -ypc z_coords];
        Plane_custom.faces = [1:3;4:6];

        %Rotate sensor plane
        Plane_custom.vertices = Plane_custom.vertices*Ry;

        %Compute intersection with RV
        [intersect, Surf] = SurfaceIntersection(Surface, Plane_custom);

        %Cut the sensor
        mask = Surf.vertices(:,2) > 31.2208;
        sensor = Surf.vertices(mask,:);

        strain_track(i) = norm(Surf.vertices,2);
        sensor_track(i) = norm(sensor,2);
    
    end
    
    %Normalize volume to compare with strain
    volume_n = normalize(volume_track);
    strain_n = normalize(strain_track);
    sensor_n = normalize(sensor_track);
    
    %Compute the error for given coordinates
    errors(rot) = rmse(volume_n,sensor_n)^2;
    if rmse(volume_n,sensor_n) < best_error
        best_error = rmse(volume_n,sensor_n);
        best_angle = angle(rot);
        best_strain = sensor_track;
    end

end

best_strain_roty = best_strain;

%% See the error curve 

plot(angle,errors)
xlabel("Rotation around y-axis")
ylabel("RMS error")

%% Visualize the results 

plot(normalize(best_strain))
hold on
plot(normalize(volume_track))
xlabel('Frame')
title('Change in strain vs change in volume')
legend('Sensor strain','Volume changes')