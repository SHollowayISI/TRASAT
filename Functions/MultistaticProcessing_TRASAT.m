function [multi] = MultistaticProcessing_TRASAT(scenario)
%MULTISTATICPROCESSING_TRASAT Multistatic processing for TRASAT project
%   Takes radar scenario object and returns list of target coordinates
%   estimated by multilateration.

%% Unpack Variables

multi = scenario.multi;

%% Perform 3-D Multilateration

% Pre-allocate data struction
multi.lat_points = nan(3, multi.n_fr);

% Break if only one transceiver reporting data
if multi.n_re < 3
    return
end

% Determine indices of transceivers, depending on method
switch scenario.multi.method
    case 'SNR'
        
        % Remove non-detection values
        sort_SNR = multi.SNR;
        sort_SNR(multi.detect == 0) = -Inf;
        
        % Sort SNR values
        [~, sort_ind] = sort(sort_SNR, 2, 'descend');
        multi.lat_index = sort(sort_ind(:,1:3), 2);
        
    case 'Range'
        
        % Remove non-detection values
        sort_range = multi.ranges;
        sort_range(multi.detect == 0) = Inf;
        
        % Sort Range values
        [~, sort_ind] = sort(sort_range, 2, 'ascend');
        multi.lat_index = sort(sort_ind(:,1:3), 2);
        
    otherwise
        % Default to first two transceivers
        multi.lat_index = ones(multi.n_fr, [1 2 3]);
end


% Loop through each frame of simulation
for frame = 1:multi.n_fr
    
    % Check if 3 radar units detected targets
    if nnz(scenario.multi.detect(frame,:)) < 3
        % Return NaN if multilateration can not be completed
        multi.lat_points(:,frame) = nan(3,1);
        
    else
        % Calculate multilateration results
        multi.lat_points(:,frame) = Multilateration3D( ...
            multi.radar_pos(:, multi.lat_index(frame, 1)), ...
            multi.ranges(frame, multi.lat_index(frame,1)), ...
            multi.radar_pos(:, multi.lat_index(frame, 2)), ...
            multi.ranges(frame, multi.lat_index(frame,2)), ...
            multi.radar_pos(:, multi.lat_index(frame, 3)), ...
            multi.ranges(frame, multi.lat_index(frame,3)));
    end
    
end

end

