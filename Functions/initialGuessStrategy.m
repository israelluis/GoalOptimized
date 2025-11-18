function initialX = initialGuessStrategy(SelectedMethod,numInitialPoints,to_gc, min_gc,lb_list,ub_list,varNames,constraintFcn)
% Method 1: lhsdesign with constraints
% Method 2: lhsdesign iteratively to generate spread control parameters
% Method 3: lhsdesign iteratively to generate diverse torque profiles
rng('shuffle');  % randomize
nNodes=3;
nVars  =length(lb_list);
% constraintFcn = @(x)(x.t_p + x.t_f <= to_gc) & (x.t_p - x.t_r >= min_gc);

% Naive
% initialPoints    = lhsdesign(numInitialPoints, nVars);
% initialX         = array2table(initialPoints .* (ub_list' - lb_list') + lb_list', 'VariableNames', varNames);
if SelectedMethod==0
    initialPoints    = lhsdesign(numInitialPoints, nVars);
    initialX         = array2table(initialPoints .* (ub_list' - lb_list') + lb_list', 'VariableNames', varNames);
elseif SelectedMethod==1
    A = [1, 0, 1, 0;  % for constraint 1: 1*x1 + 0*x2 + 1*x3 + 0*x4 <= to_gc
        -1, 1, 0, 0]; % for constraint 2:(-1)*x1 + 1*x2 + 0*x3 + 0*x4 <= -15
    b = [to_gc; -min_gc];
    lb = lb_list'; % Lower bounds as a row vector
    ub = ub_list'; % Upper bounds as a row vector
    initialPoints = lhsdesigncon(numInitialPoints, nVars, lb, ub, [false false],A, b);
    initialX      = array2table(initialPoints, 'VariableNames', varNames);

elseif SelectedMethod==2
    numCandidates   = numInitialPoints * 10; % Generate 10x more points
    candidatePoints = lhsdesign(numCandidates, nVars);
    candidateX      = array2table(candidatePoints .* (ub_list' - lb_list') + lb_list', 'VariableNames', varNames);

    % Filter feasible points
    feasibleMask = false(height(candidateX), 1);
    for i = 1:height(candidateX)
        feasibleMask(i) = constraintFcn(candidateX(i,:));
    end
    feasibleX = candidateX(feasibleMask, :);

    % Select the most spread points using maximin distance
    if height(feasibleX) >= numInitialPoints
        initialX = selectSpreadPoints(feasibleX, numInitialPoints);
    else
        warning('Only %d feasible points found, using all', height(feasibleX));
        initialX = feasibleX;
    end
elseif SelectedMethod==3
    numCandidates   = numInitialPoints * 10;
    candidatePoints = lhsdesign(numCandidates, nVars);
    candidateX      = array2table(candidatePoints .* (ub_list' - lb_list') + lb_list', 'VariableNames', varNames);

    % Filter feasible points
    feasibleMask = false(height(candidateX), 1);
    for i = 1:height(candidateX)
        feasibleMask(i) = constraintFcn(candidateX(i,:));
    end
    feasibleX = candidateX(feasibleMask, :);

    % Generate moment_spline for all feasible candidates
    fprintf('Generating moment_spline profiles for %d feasible candidates...\n', height(feasibleX));
    moment_profiles = zeros(100, height(feasibleX)); % Assuming 101 time points

    for i = 1:height(feasibleX)
        point = feasibleX(i,:);
        [moment_spline, ~, ~, ~] = generateSplineBasedAssistance(nNodes, point{1,1:3}, point{1,4}, 0);
        moment_profiles(:, i) = moment_spline;
    end

    % Select points that maximize moment_spline diversity
    if height(feasibleX) >= numInitialPoints
        initialX = selectDiverseProfiles(feasibleX, moment_profiles, numInitialPoints);
    else
        warning('Only %d feasible points found, using all', height(feasibleX));
        initialX = feasibleX;
    end

elseif SelectedMethod == 4

    % Hybrid: Stratified torque sampling with profile diversity within each segment
    numPerSegment = round(numInitialPoints / 3);
    
    % Define torque range segments for var4
    torque_min = lb_list(4);
    torque_max = ub_list(4);
    torque_ranges = [torque_min, torque_min + (torque_max-torque_min)/3; 
                    torque_min + (torque_max-torque_min)/3, torque_min + 2*(torque_max-torque_min)/3;
                    torque_min + 2*(torque_max-torque_min)/3, torque_max];
    
    initialX_segments = cell(1,3);
    
    for segment = 1:3
        fprintf('Generating diverse profiles for torque segment %d...\n', segment);
        
        % Generate candidates for this torque segment
        numCandidates = numPerSegment * 10;
        candidatePoints = lhsdesign(numCandidates, nVars);
        candidateX = array2table(candidatePoints .* (ub_list' - lb_list') + lb_list', 'VariableNames', varNames);
        
        % Filter by constraints AND torque range
        feasibleMask = false(height(candidateX), 1);
        moment_profiles = zeros(100, height(candidateX));
        
        for i = 1:height(candidateX)
            point = candidateX(i,:);
            torque_value = point.(varNames{4});
            feasibleMask(i) = constraintFcn(point) && ...
                             (torque_value >= torque_ranges(segment,1)) && ...
                             (torque_value <= torque_ranges(segment,2));
            
            % Generate torque profile for diversity selection
            if feasibleMask(i)
                [moment_spline, ~, ~, ~] = generateSplineBasedAssistance(nNodes, point{1,1:3}, point{1,4}, 0);
                moment_profiles(:, i) = moment_spline;
            end
        end
        
        feasibleX = candidateX(feasibleMask, :);
        feasible_profiles = moment_profiles(:, feasibleMask);
        
        % Select diverse profiles within this torque segment
        if height(feasibleX) >= numPerSegment
            initialX_segments{segment} = selectDiverseProfiles(feasibleX, feasible_profiles, numPerSegment);
        else
            warning('Segment %d: Only %d feasible points found', segment, height(feasibleX));
            initialX_segments{segment} = feasibleX;
        end
    end
    
    % Combine segments
    initialX = vertcat(initialX_segments{:});
    
    % Handle size adjustments (same as before)
    if height(initialX) > numInitialPoints
        rng(42);
        randomIndices = randperm(height(initialX), numInitialPoints);
        initialX = initialX(randomIndices, :);
    end
    
    if height(initialX) < numInitialPoints
        additionalNeeded = numInitialPoints - height(initialX);
        additionalX = generateAdditionalPoints(additionalNeeded, to_gc, min_gc, lb_list, ub_list, varNames, nVars, constraintFcn);
        initialX = [initialX; additionalX];
    end

elseif SelectedMethod == 5
    % Strategy 5: Discretize both peak torque (M_p) and peak time (t_p) into 3 segments each
    % Creates 3x3 = 9 distinct regions with equal points per region
    
    fprintf('Strategy 5: Equal points across 9 torque-time regions\n');
    
    % Define segments for peak torque (var4) and peak time (var2)
    torque_min = lb_list(4);
    torque_max = ub_list(4);
    time_min = lb_list(1);
    time_max = ub_list(1);
    
    % Create 3x3 grid of regions
    torque_ranges = [torque_min, torque_min + (torque_max-torque_min)/3; 
                    torque_min + (torque_max-torque_min)/3, torque_min + 2*(torque_max-torque_min)/3;
                    torque_min + 2*(torque_max-torque_min)/3, torque_max];
    
    time_ranges = [time_min, time_min + (time_max-time_min)/3;
                  time_min + (time_max-time_min)/3, time_min + 2*(time_max-time_min)/3;
                  time_min + 2*(time_max-time_min)/3, time_max];
    
    % Calculate maximum equal points per region
    points_per_region = floor(numInitialPoints / 9);
    total_points = points_per_region * 9;
    
    if total_points < numInitialPoints
        fprintf('Using %d points (%d per region × 9 regions)\n', total_points, points_per_region);
    else
        fprintf('Using %d points (%d per region × 9 regions)\n', numInitialPoints, points_per_region);
    end
    
    initialX_regions = cell(3, 3); % 3x3 grid for torque x time
    region_success = false(3, 3);
    
    fprintf('Generating %d points for each of 9 regions:\n', points_per_region);
    
    for torque_seg = 1:3
        for time_seg = 1:3
            fprintf('  Region %d-%d (Tq[%.1f-%.1f], Tm[%.1f-%.1f]): ', ...
                torque_seg, time_seg, ...
                torque_ranges(torque_seg,1), torque_ranges(torque_seg,2), ...
                time_ranges(time_seg,1), time_ranges(time_seg,2));
            
            region_points = [];
            attempts = 0;
            max_attempts = 50;
            
            % Keep trying until we get the required points for this region
            while height(region_points) < points_per_region && attempts < max_attempts
                attempts = attempts + 1;
                
                % Generate fresh candidates for this specific region
                needed = (points_per_region - height(region_points)) * 10;
                candidatePoints = lhsdesign(needed, nVars);
                candidateX = array2table(candidatePoints .* (ub_list' - lb_list') + lb_list', 'VariableNames', varNames);
                
                % Filter by constraints AND region boundaries
                feasibleMask = false(height(candidateX), 1);
                for i = 1:height(candidateX)
                    point = candidateX(i,:);
                    torque_value = point.(varNames{4});
                    time_value = point.(varNames{1});
                    
                    feasibleMask(i) = constraintFcn(point) && ...
                                     (torque_value >= torque_ranges(torque_seg,1)) && ...
                                     (torque_value <= torque_ranges(torque_seg,2)) && ...
                                     (time_value >= time_ranges(time_seg,1)) && ...
                                     (time_value <= time_ranges(time_seg,2));
                end
                
                new_feasible = candidateX(feasibleMask, :);
                
                if ~isempty(new_feasible)
                    if isempty(region_points)
                        region_points = new_feasible;
                    else
                        region_points = [region_points; new_feasible];
                    end
                    % Remove duplicates
                    region_points = unique(region_points, 'rows');
                end
                
                fprintf('.');
            end
            
            % Select exactly points_per_region with maximum torque diversity
            if height(region_points) >= points_per_region
                % Maximize torque diversity within this region
                initialX_regions{torque_seg, time_seg} = selectMaxTorqueDiversity(region_points, points_per_region, varNames);
                region_success(torque_seg, time_seg) = true;
                fprintf(' %d points ✓\n', points_per_region);
            else
                % Use all we could find
                initialX_regions{torque_seg, time_seg} = region_points;
                fprintf(' %d points (wanted %d)\n', height(region_points), points_per_region);
            end
        end
    end
    
    % Combine all regions into single table
    initialX = [];
    for torque_seg = 1:3
        for time_seg = 1:3
            if ~isempty(initialX_regions{torque_seg, time_seg})
                initialX = [initialX; initialX_regions{torque_seg, time_seg}];
            end
        end
    end
    
    % Final distribution analysis
    torque_vals = initialX.(varNames{4});
    time_vals = initialX.(varNames{1});
    
    fprintf('\n=== FINAL DISTRIBUTION ===\n');
    fprintf('Total points: %d\n', height(initialX));
    fprintf('Points per region (ideal: %d):\n', points_per_region);
    
    for torque_seg = 1:3
        for time_seg = 1:3
            count = sum(torque_vals >= torque_ranges(torque_seg,1) & ...
                       torque_vals <= torque_ranges(torque_seg,2) & ...
                       time_vals >= time_ranges(time_seg,1) & ...
                       time_vals <= time_ranges(time_seg,2));
            if region_success(torque_seg, time_seg)
            status = '✓'; %region_success(torque_seg, time_seg) 
            else 
            status =  '⚠';
            end
            fprintf('  Torque%d-Time%d: %d points %s\n', torque_seg, time_seg, count, status);
        end
    end
    
    % Calculate torque diversity metric
    torque_diversity = calculateTorqueDiversity(initialX, varNames);
    fprintf('Overall torque diversity: %.3f (1.0 = perfect)\n', torque_diversity);
end
end

% Supplementary functions
function selected = selectSpreadPoints(points, nSelect)
    % Select points that maximize minimum distance between them
    nTotal = height(points);
    pointMatrix = table2array(points);
    
    % Start with a random point
    selectedIndices = randi(nTotal);
    
    for i = 2:nSelect
        % Calculate minimum distance from each point to already selected points
        minDistances = zeros(nTotal, 1);
        for j = 1:nTotal
            if ~ismember(j, selectedIndices)
                dists = sqrt(sum((pointMatrix(j,:) - pointMatrix(selectedIndices,:)).^2, 2));
                minDistances(j) = min(dists);
            else
                minDistances(j) = -inf;
            end
        end
        
        % Select point with maximum minimum distance
        [~, nextIdx] = max(minDistances);
        selectedIndices = [selectedIndices; nextIdx];
    end
    
    selected = points(selectedIndices, :);
end

function selected = selectDiverseProfiles(points, profiles, nSelect)
    % Select points that maximize diversity of moment_spline profiles
    % points: table of parameter combinations
    % profiles: matrix where each column is a moment_spline profile
    % nSelect: number of points to select
    
    nTotal = height(points);
    selectedIndices = [];
    
    % Start with the profile that's most different from mean profile
    mean_profile = mean(profiles, 2);
    distances_to_mean = zeros(nTotal, 1);
    for i = 1:nTotal
        distances_to_mean(i) = norm(profiles(:, i) - mean_profile);
    end
    [~, firstIdx] = max(distances_to_mean);
    selectedIndices = firstIdx;
    
    % Greedy selection: at each step, pick the profile that is most different
    % from all already selected profiles
    for step = 2:nSelect
        max_min_distance = -inf;
        best_candidate = -1;
        
        for i = 1:nTotal
            if ~ismember(i, selectedIndices)
                % Calculate minimum distance to any already selected profile
                min_dist = inf;
                for j = 1:length(selectedIndices)
                    dist = norm(profiles(:, i) - profiles(:, selectedIndices(j)));
                    if dist < min_dist
                        min_dist = dist;
                    end
                end
                
                if min_dist > max_min_distance
                    max_min_distance = min_dist;
                    best_candidate = i;
                end
            end
        end
        
        if best_candidate > 0
            selectedIndices = [selectedIndices; best_candidate];
        else
            break; % No more candidates
        end
    end
    
    selected = points(selectedIndices, :);
    
    fprintf('Selected %d points with diverse assistance profiles\n', length(selectedIndices));
end

function additionalX = generateAdditionalPoints(numNeeded, to_gc, min_gc, lb_list, ub_list, varNames, nVars, constraintFcn)
    % Generate additional points when segments don't provide enough
    numCandidates = numNeeded * 20;
    candidatePoints = lhsdesign(numCandidates, nVars);
    candidateX = array2table(candidatePoints .* (ub_list' - lb_list') + lb_list', 'VariableNames', varNames);
    
    % Filter feasible points
    feasibleMask = false(height(candidateX), 1);
    for i = 1:height(candidateX)
        feasibleMask(i) = constraintFcn(candidateX(i,:));
    end
    feasibleX = candidateX(feasibleMask, :);
    
    if height(feasibleX) >= numNeeded
        additionalX = selectSpreadPoints(feasibleX, numNeeded);
    else
        warning('Could only generate %d additional points (needed %d)', height(feasibleX), numNeeded);
        additionalX = feasibleX;
    end
end

function selected = selectMaxTorqueDiversity(points, nSelect, varNames)
    % Select points that maximize torque diversity within the region
    % This ensures we get the full range of torque values in each time segment
    
    nTotal = height(points);
    if nTotal <= nSelect
        selected = points;
        return;
    end
    
    % Extract torque values
    torque_vals = points.(varNames{4});
    
    % Strategy: Select points that are evenly spaced across the torque range
    sorted_torque = sort(torque_vals);
    
    % Create target torque values that are evenly distributed
    target_torques = linspace(min(torque_vals), max(torque_vals), nSelect);
    
    selectedIndices = [];
    
    for i = 1:nSelect
        % Find point closest to target torque value
        [~, closest_idx] = min(abs(torque_vals - target_torques(i)));
        
        % Make sure we don't select duplicates
        if ~ismember(closest_idx, selectedIndices)
            selectedIndices = [selectedIndices; closest_idx];
        else
            % If duplicate, find next closest unique point
            available_indices = setdiff(1:nTotal, selectedIndices);
            if ~isempty(available_indices)
                available_torques = torque_vals(available_indices);
                [~, rel_idx] = min(abs(available_torques - target_torques(i)));
                selectedIndices = [selectedIndices; available_indices(rel_idx)];
            end
        end
    end
    
    selected = points(selectedIndices, :);
end

function diversity = calculateTorqueDiversity(points, varNames)
    % Calculate how well torque values are distributed (0-1, where 1 is perfect)
    
    torque_vals = points.(varNames{4});
    
    if length(torque_vals) < 2
        diversity = 0;
        return;
    end
    
    % Calculate evenness of distribution using histogram evenness
    num_bins = min(10, length(unique(torque_vals)));
    [counts, ~] = histcounts(torque_vals, num_bins);
    
    % Perfect distribution would have equal counts in all bins
    expected_count = length(torque_vals) / num_bins;
    chi_squared = sum((counts - expected_count).^2 / expected_count);
    
    % Convert to 0-1 scale (0 = uneven, 1 = perfectly even)
    diversity = 1 / (1 + chi_squared / num_bins);
end