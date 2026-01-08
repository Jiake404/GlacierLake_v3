function [DEM,fill_flag,New_r,New_c,Remaining_dist_matrix,external_water,s_flag,total_external_input]=flow_speed2(DEM,s_rol,ro_water,s_mwe,...
    t_step,runoff,tt,s_tnum,New_r,New_c,Remaining_dist_matrix,fill_flag,DEM_ini,external_water,s_flag,tr_a,cell_size,spilloverRow,spilloverCol)

gs = 0.001; %mean grain size,m
mu = 1.763*10^(-3); % viscosity of water,(Pa s)= 1.76*10^-3 Leeson et al.,2012;=1.8*10^-3 Arnold 1998 
n = 0.05; %Leeson et al.,2012
dx = cell_size; %m
Rh = 0.035;
s_thro = 0.1;%mwe;

[rows, cols] = size(DEM);
flowdir = zeros(rows, cols); % flow direction
flowdir(isnan(DEM))=nan;
slope = zeros(rows, cols); % slope
slope(isnan(DEM))=nan;

%%Calculate the flow direction
% Define eight directions
drow = [-1, -1, 0, 1, 1, 1, 0, -1]; % 
dcol = [0, 1, 1, 1, 0, -1, -1, -1]; % 

% Define the distances of the eight directions ( the grid size or the diagonal distance)
distances = [dx, dx*sqrt(2), dx, dx*sqrt(2), dx, dx*sqrt(2), dx, dx*sqrt(2)];
D8_codes = [1,2,3,4,5,6,7,8];

% for each grid
for i = 1:rows
    for j = 1:cols
        if ~isnan(DEM(i,j))
            min_elev = DEM(i,j);
            dir = 0;
            % for each directions
            for k = 1:8
                % Calculate the row and column indices of the target grid
                target_i = i + drow(k);
                target_j = j + dcol(k);
                if target_i >= 1 && target_i <= rows && target_j >= 1 && target_j <= cols
                    if DEM(target_i, target_j) <= min_elev
                        min_elev = DEM(target_i, target_j);
                        best_dir = k;
                    end
                end
            end
            flowdir(i, j) = best_dir;
            % calculate slope
            slope(i,j) = (DEM(i,j)-DEM(i+drow(flowdir(i,j)),j+dcol(flowdir(i,j)))) / distances(flowdir(i,j));
        end

    end
end
slope(slope<0)=0;

% Find all the lowest points (sinks)
min_elev = min(DEM(:));
[sink_rows, sink_cols] = find(DEM == min_elev); 

%% Calculate the flow rate
if tr_a == 2  %snow
    k = 0.0778*gs^2*exp(-7.8*s_rol/ro_water);%water permeability,Shimizu's equation,1969
    v = k.*slope./mu;
else
    v = Rh^(2/3).*slope.^0.5.*1/n;
end
%% Calculate the distance from each grid to the sink.
    % Initialize the step matrix, -1 indicates that the location has not been visited.
    distance = -1 * ones(rows, cols);
    queue = [sink_rows, sink_cols];  % BFS 队列
    for i = 1:length(sink_rows)
        distance(sink_rows(i), sink_cols(i)) = 0; 
    end
    while ~isempty(queue)
        r = queue(1,1);
        c = queue(1,2);
        queue(1, :) = []; 

        for d = 1:8
            nr = r + drow(d);
            nc = c + dcol(d);
            if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols
                expected_dir = D8_codes(mod(d+3, 8) + 1);
                if flowdir(nr, nc) == expected_dir && distance(nr, nc) == -1
                    step_dist = dx; 
                    if mod(d,2) == 0 
                        step_dist = dx * sqrt(2);
                    end

                    distance(nr, nc) = distance(r, c) + step_dist;
                    queue = [queue; nr, nc]; 
                end
            end
        end
    end

%% Calculate the total amount entering the sink
W_move=zeros([size(DEM),tt]);
Distance_new=zeros([size(DEM),tt]);
water_depth=0;
for layer = s_tnum:tt
    new_r = squeeze(New_r(:,:,layer));
    new_c = squeeze(New_c(:,:,layer));
    remaining_dist_matrix = squeeze(Remaining_dist_matrix(:,:,layer));
    for r = 1:rows
        for c = 1:cols
            w_move(r,c) = v(new_r(r,c),new_c(r,c)).*3600.*t_step;
            distance_new(r,c) = distance(new_r(r,c),new_c(r,c))-remaining_dist_matrix(r,c);
        end
    end
    grid_num = sum((w_move(:) > distance_new(:)) & (distance_new(:)>0));
    w_move(w_move>distance_new)=distance_new(w_move > distance_new);
    water_depth_one = grid_num*runoff(layer);
    water_depth = water_depth_one+water_depth;
    W_move(:,:,layer)=w_move;
    Distance_new(:,:,layer)=distance_new;
end

%% filling
while water_depth > 0
    neighbors = [];
    [sink_rows, sink_cols] = find(DEM == min(DEM(:))); 
    for i = 1:length(sink_rows)
        for j = 1:8
            nr = sink_rows(i,1) + drow(j);
            nc = sink_cols(i,1) + dcol(j);
            if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols && fill_flag(nr, nc) == 0
                neighbors = [neighbors; DEM(nr, nc), nr, nc];
               
                [~, unique_idx] = unique(neighbors(:, 2:3), 'rows', 'stable');
                
                neighbors_unique = neighbors(unique_idx, :);
            end
        end
    end
    if isempty(neighbors_unique)
        continue;
    end
    % Remove the row and column numbers of the sink itself from the "neighbors_unique" array.
    sink_coords = [sink_rows, sink_cols];  
    neighbors_coords = neighbors_unique(:, 2:3);  

    % Check if the "sink" itself is included in the "neighbors_unique" list
    sink_idx = ismember(neighbors_coords, sink_coords, 'rows');

    % remove sink
    neighbors_unique(sink_idx, :) = [];

    neighbors_unique = sortrows(neighbors_unique, 1);  
    next_height = neighbors_unique(1, 1); 

    current_height(1) = DEM(sink_rows(1),sink_cols(1));
    required_water_depth = sum(next_height - current_height)*length(sink_rows);
     % 5. Check whether the water volume is sufficient
    if water_depth >= required_water_depth
        % The water volume is sufficient and will be filled up to next_height.
        DEM(sub2ind([rows, cols], sink_rows, sink_cols)) = next_height;
        water_depth = water_depth - required_water_depth;
        
    else
        % Insufficient water volume, average filling
        final_height = current_height + water_depth / length(sink_rows);
        DEM(sub2ind([rows, cols], sink_rows, sink_cols)) = final_height;
        water_depth = 0;  
        
    end
    for i = 1:length(sink_rows)
        fill_flag(sink_rows(i,1),sink_cols(i,1)) = 1;
    end
end
%% Track the position of the water in each grid after a step, and calculate for each layer.
for layer = s_tnum:tt
    w_move = squeeze(W_move(:,:,layer));
    remaining_dist_matrix = squeeze(Remaining_dist_matrix(:,:,layer));
    remaining_dist_matrix = remaining_dist_matrix + w_move;
    new_r = squeeze(New_r(:,:,layer));
    new_c = squeeze(New_c(:,:,layer));
    for r = 1:rows
        for c = 1:cols
            if w_move(r, c) > 0  
                remaining_dist = remaining_dist_matrix(r, c);
                curr_r = new_r(r,c);
                curr_c = new_c(r,c);

                % Track the movement of water
                while remaining_dist > 0
                    flow_d = flowdir(curr_r, curr_c); % The flow direction of the current grid
                    d_idx = find(D8_codes == flow_d); % Obtain the corresponding D8 index

                    if isempty(d_idx)  
                        break;
                    end

                    % Calculate the position of the next grid.
                    next_r = curr_r + drow(d_idx);
                    next_c = curr_c + dcol(d_idx);
                
                    % Determine whether it is within the catchment
                    if next_r < 1 || next_r > rows || next_c < 1 || next_c > cols
                        break;
                    end

                    % 
                    step_dist = dx; 
                    if mod(d_idx, 2) == 0  
                        step_dist = dx * sqrt(2);
                    end

                    % Determine whether there is still remaining distance to cover
                    if remaining_dist < step_dist
                        remaining_dist_matrix(r, c) = remaining_dist;
                        break;
                    else
                        % Move to the next square
                        curr_r = next_r;
                        curr_c = next_c;
                        remaining_dist = remaining_dist - step_dist;
                    end
                end
                % record the new position
                new_r(r, c) = curr_r;
                new_c(r, c) = curr_c;
                remaining_dist_matrix(r, c) = remaining_dist; % 记Record the remaining distance that has not been covered
            end
        end
    end
    Remaining_dist_matrix(:,:,layer) = remaining_dist_matrix;
    New_r(:,:,layer)=new_r;
    New_c(:,:,layer)=new_c;
end
%% Calculate the inflow of external meltwater
dem_diff = DEM-DEM_ini;
dem_diff(dem_diff<0)=0;
total_external_input = floor(dem_diff./0.1);
% nansum(total_external_input(:))
s_flag(total_external_input>0)=1;
external_input = total_external_input-nansum(external_water(:,:,1:tt-1),3);
if any(external_input(:)<0)
    disp('The lake depth has experienced a negative growth!!!')
end
external_water(:,:,tt)=external_input;
if fill_flag(spilloverRow,spilloverCol)==1
    disp('Overflow!!!')
    disp(tt)
end









