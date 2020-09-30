function [results] = ErrorEstimation_TRASAT(scenario)
%ERRORESTIMATION_TRASAT Estimate errors from TRASAT simulations
%   Takes radar scenario object as input, returns modified scenario.results
%   object as output.

%% Unpack Variables

traj = scenario.traj;
multi = scenario.multi;
simsetup = scenario.simsetup;

%% Unpack trajectory positions

% Calculate points in trajectory
traj_ind = round(simsetup.n_p * simsetup.cpi_fr * ...
    (0.5:(multi.n_fr-0.5)) + 1);
traj_points = traj.pos(:,traj_ind);

% Determine true range and angle from each receiver
true_ranges = zeros(length(traj_ind), multi.n_re);
true_angles = zeros(length(traj_ind), multi.n_re);
for fr = 1:length(traj_ind)
    for re = 1:multi.n_re
        [true_ranges(fr, re), ang] = rangeangle(traj_points(:,fr), multi.radar_pos(:,re));
        true_angles(fr, re) = ang(2);
    end
end

%% Calculate monostatic errors

% Determine range error for each frame
results.range_error = true_ranges - multi.ranges;
results.range_error_rms = squeeze(rms(results.range_error, 1));

% Determine angle error for each frame
results.angle_error = true_angles - multi.angles - multi.gimbal_angles;
results.angle_error_rms = squeeze(rms(results.angle_error, 1));

%% Calculate multistatic errors

% Determine trilateration error for each frame
results.trilat_error = multi.lat_points - traj_points;
results.trilat_error_rms = rms(results.trilat_error, 2);

% Determine trilateration error with tracking
results.track_error = multi.track_points - traj_points;
results.track_error_rms = rms(results.track_error, 2);

%% Calculate errors for triangle tests

if mod(multi.n_re, 3) == 0
    results.range_error_tri = rms(reshape(results.range_error_rms, [], 3),2);
    results.angle_error_tri = rms(reshape(results.angle_error_rms, [], 3),2);
    results.trilat_error_tri = squeeze(rms(multi.lat_points_tri - traj_points, 2));
end

end

