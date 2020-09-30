function [multi] = MultistaticProcessing_Triangular(scenario)
%MULTISTATICPROCESSING_TRIANGULAR Multistatic processing for triangular
%tests
%   Takes radar scenario object and returns list of target coordinates
%   estimated by multilateration.

%% Unpack Variables

multi = scenario.multi;

%% Perform 3-D Multilateration

% Pre-allocate data struction
multi.lat_points_tri = nan(3, multi.n_fr, multi.n_re/3);

% Loop through each triangle of three receivers
for tri_num = 1:(multi.n_re/3)
    
    for frame = 1:multi.n_fr
        
        % Perform multilateration
        multi.lat_points_tri(:,frame, tri_num) = Multilateration3D( ...
            multi.radar_pos(:, 3*(tri_num-1) + 1), ...
            multi.ranges(frame, 3*(tri_num-1) + 1), ...
            ...
            multi.radar_pos(:, 3*(tri_num-1) + 2), ...
            multi.ranges(frame, 3*(tri_num-1) + 2), ...
            ...
            multi.radar_pos(:, 3*(tri_num-1) + 3), ...
            multi.ranges(frame, 3*(tri_num-1) + 3));
    
    end
end

end

