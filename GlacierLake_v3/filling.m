function [DEM]=filling(DEM)
[rows, cols] = size(DEM);
flowdir = zeros(rows, cols); %flow direction
flowdir(isnan(DEM))=nan;
slope = nan(rows, cols); % slope
% slope(isnan(DEM))=nan;


%%calculate flow direction
% Define eight directions
drow = [-1, -1, 0, 1, 1, 1, 0, -1]; 
dcol = [0, 1, 1, 1, 0, -1, -1, -1]; 
dx=100;
% Define the distances of the eight directions ( the grid size or the diagonal distance)
distances = [dx, dx*sqrt(2), dx, dx*sqrt(2), dx, dx*sqrt(2), dx, dx*sqrt(2)];
D8_codes = [1,2,3,4,5,6,7,8];

% for each grid
for i = 1:rows
    for j = 1:cols
        if ~isnan(DEM(i,j))
            min_elev = DEM(i,j);
            best_dir = 0;
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

            if flowdir(i,j)==0
                slope(i,j)=0;
            else
                slope(i,j) = (DEM(i,j)-DEM(i+drow(flowdir(i,j)),j+dcol(flowdir(i,j)))) / distances(flowdir(i,j));
            end
        end
    end
end
%% filling
[min_val, ~] = min(DEM(:));  %the lowest value
iteration_count = 0;  % Record the number of iterations

while any(slope(:) == 0 & DEM(:) > min_val)

    [row_idx, col_idx] = find(slope == 0);  
    if isempty(row_idx) 
        break;
    end
    
    for i = 1:length(row_idx)
        if DEM(row_idx(i), col_idx(i)) ~= min_val
            prep_DEM = [];
            for k = 1:8
                target_i = row_idx(i) + drow(k);
                target_j = col_idx(i) + dcol(k);
                if target_i >= 1 && target_i <= rows && target_j >= 1 && target_j <= cols
                    prep_DEM(k, 1) = DEM(target_i, target_j);
                end
            end
            [min_DEM, ~] = min(prep_DEM); 
            if ~isnan(min_DEM)
                DEM(row_idx(i), col_idx(i)) = min_DEM + 0.1;  % 填洼
            else
                disp('There is insufficient valid data, which makes it impossible to fill the gap.--filling.m')
            end
        end
    end

    % Re-calculate the flow direction and slope
    for i = 1:rows
        for j = 1:cols
            if ~isnan(DEM(i, j))
                min_elev = DEM(i, j);
                best_dir = 0;
                % for each direction
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
                if flowdir(i,j)==0
                    slope(i,j)=0;
                else
                    slope(i,j) = (DEM(i,j)-DEM(i+drow(flowdir(i,j)),j+dcol(flowdir(i,j)))) / distances(flowdir(i,j));
                end
            end
        end
    end
    
    % Check if the slope of the lowest point is less than 0.
    if sum(slope(:) < 0) == 1  
        slope(slope < 0) = 0;
        disp('Filling successfully √ ----filling.m');
        break;
    end
    
    iteration_count = iteration_count + 1;
    if iteration_count > 10000  
        disp('Filling failed. The number of iterations exceeded the limit.----filling.m');
        break;
    end
end


% Find all the lowest points (sinks)
min_elev = min(DEM(:));
[sink_rows, sink_cols] = find(DEM == min_elev); 
% Initialize the marking matrix
reaches_sink = false(rows, cols);
queue = [sink_rows, sink_cols]; 

% Mark all "sink" entries as "true"
for k = 1:size(queue, 1)
    reaches_sink(queue(k, 1), queue(k, 2)) = true;
end

% 
for i = 1:rows
    for j = 1:cols
        if reaches_sink(i, j) || flowdir(i, j) == 0
            continue; %
        end
        
        % Track the flow path
        path = []; % Record the grids on the current path
        [cur_i, cur_j] = deal(i, j);
        while true
            path = [path; cur_i, cur_j]; % record the path
            direction = flowdir(cur_i, cur_j);
            if direction == 0
                break; 
            end
            
            % Find the next flow direction pixel
            d = find(D8_codes == direction, 1);
            if isempty(d)
                break;
            end
            cur_i = cur_i + drow(d);
            cur_j = cur_j + dcol(d);

            % Is it within the catchment
            if cur_i < 1 || cur_i > rows || cur_j < 1 || cur_j > cols
                break;
            end

            % If one has reached the grid at the lowest point of the flow direction
            if reaches_sink(cur_i, cur_j)
                % All the grids on the path are marked as true
                for p = 1:size(path, 1)
                    reaches_sink(path(p, 1), path(p, 2)) = true;
                end
                break;
            end
        end
    end
end

if any(~reaches_sink(~isnan(DEM)))
    disp('Not all the grids flow towards the sink!----filling.m')
else
    disp('All the grids flow towards the sink. √----filling.m')
end


