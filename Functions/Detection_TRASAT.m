function [detection] = Detection_TRASAT(scenario)
%DETECTION_TRASAT Performs target detection for TRASAT project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing information about detected targets.

%% Unpack Variables

simsetup = scenario.simsetup;
cube = scenario.cube;

%% Perform Detection

% Estimate noise power
noise_pow = pow2db(median(mean(scenario.cube.pow_cube,1)));

% Calculate threshold in absolute
abs_thresh = db2pow(simsetup.thresh + noise_pow);

% Perform detection
detection.detect_cube_CPI = (cube.pow_cube_CPI > abs_thresh);


%% Perform m-of-n Detection

% Determine bins with m-of-n detections
detection.detect_cube = (sum(double(detection.detect_cube_CPI),3) >= simsetup.det_m);

% Determine if any target is detected
detection.detect_logical = logical(nnz(detection.detect_cube) > 0);

% Determine location of maximum power bin
[max_pow, max_ind] = max(cube.pow_cube, [], 'all', 'linear');
[detection.max_r, detection.max_v] = ...
    ind2sub(size(cube.pow_cube), max_ind);

% Save max SNR value
detection.max_SNR = pow2db(max_pow) - noise_pow;

%% Estimate Target Location

if detection.detect_logical
    
    % OVERWRITE: Using old method of range/velocity estimation
%     [r_idx, d_idx] = find(detection.detect_cube);
%     
%     detection.target_range = cube.range_axis(ceil(median(r_idx)));
%     detection.target_vel = cube.vel_axis(ceil(median(d_idx)));
%     
    % Find objects in detection image
%     cc = bwconncomp(detection.detect_cube);
%     
%     % Find region with max SNR
%     for n = 1:cc.NumObjects
%         
%         if any(cc.PixelIdxList{n} == max_ind)
%             cent_idx = n;
%             break
%         end
%         
%     end
    
    % Save centroid of target with max SNR
%     rp = regionprops(cc, 'Centroid');
%     cent_loc = rp(cent_idx).Centroid;
%     
%     detection.target_range = interp1( ...
%         1:length(cube.range_axis), cube.range_axis, cent_loc(2));
%     
%     detection.target_vel = interp1( ...
%         1:length(cube.vel_axis), cube.vel_axis, cent_loc(1));
    
    % OVERWRITE: Save location of maximum
    detection.target_range = cube.range_axis(detection.max_r);
    detection.target_vel = cube.vel_axis(detection.max_v);
    
    % Alternative SNR calculations
%     sum_pow = sum(cube.pow_cube(cc.PixelIdxList{cent_idx}));
%     detection.sum_SNR = pow2db(sum_pow) - noise_pow;
%     detection.mean_SNR = pow2db(sum_pow/length(cc.PixelIdxList{cent_idx})) - ...
%         noise_pow;
    
else
    
    % Fill in blank values if no target is detected
    detection.target_range = [];
    detection.target_vel = [];
    detection.target_SNR = [];
    
end

end

