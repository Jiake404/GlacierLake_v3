function [DEM,fill_flag,external_water,s_flag,total_external_input]=flow_speed0(DEM,tt,DEM_ini,external_water,s_flag,water_depth,fill_flag,...
    spilloverRow,spilloverCol)
%% filling
drow = [-1, -1, 0, 1, 1, 1, 0, -1]; 
dcol = [0, 1, 1, 1, 0, -1, -1, -1]; 
[rows, cols] = size(DEM);
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
%     for i = 1:length(sink_rows)
%         current_height(i,1) = DEM(sink_rows(i,1),sink_cols(i,1));
%     end
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

%% Calculate the inflow of external meltwater
dem_diff = DEM-DEM_ini;
dem_diff(dem_diff<0)=0;
total_external_input = floor(dem_diff./0.1);
% nansum(total_external_input(:))
s_flag(total_external_input>0)=1;
external_input = total_external_input-sum(external_water(:,:,1:tt-1),3);
if any(external_input(:)<0)
    disp('The lake depth has experienced a negative growth!!!')
end
external_water(:,:,tt)=external_input;
if fill_flag(spilloverRow,spilloverCol)==1
    disp('Overflow!!!')
    disp(tt)
end









