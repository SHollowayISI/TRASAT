function [multi] = Monopulse_TRASAT(scenario, method)
%MONOPULSE_TRASAT Performs monopulse processing for TRASAT project
%   Takes radar scenario object as input, passes scenario.multi as output,
%   containing angle data. 'method' input can be set as 'mean', 'median',
%   'mode', or 'histogram'.

%% Unpack Variables

simsetup = scenario.simsetup;
detection = scenario.detection;
multi = scenario.multi;
flags = scenario.flags;
rx_sig = scenario.rx_sig;

%% Perform Monopulse Processing

if detection.detect_logical
    
    % Find phase difference between two antennas
    phase_del = angle(squeeze(rx_sig(detection.max_r,:,1)./rx_sig(detection.max_r,:,2)));
    angle_del = asind(phase_del/(simsetup.d_ant*2*pi));
    
    switch method
        case 'mean'
            % Takes mean of angles
            mono_angle = mean(angle_del);
        case 'median'
            % Takes median of angles
            mono_angle = median(angle_del);
        case 'mode'
            % Takes mode of rounded angles
            mono_angle = mode(round(angle_del/simsetup.mono_res)*simsetup.mono_res);
        case 'histogram'
            % Finds peak in distribution of angles
            % Old code, used when distribution seemed to be bimodal, due to
            % coding error.
            
            % Determine limiting angle
            lim_ang = ceil(asind(1/(2*simsetup.d_ant))/simsetup.mono_res)*simsetup.mono_res;
            
            % Sample histogram of phase difference
            [hist_N, hist_edge] = histcounts(angle_del, -lim_ang:simsetup.mono_res:lim_ang);
            
            % Smooth histogram curve
            hist_N_ave = movmean(hist_N, simsetup.mono_smooth);
            
            % Estimate angle from histogram peak
            [~, mono_ind] = max(hist_N_ave);
            mono_angle = hist_edge(mono_ind);
            %}
    end
    
    % Save monopulse angle estimation
    multi.angles(flags.frame, flags.unit) = mono_angle;
    
    % Save elevation estimation
    multi.elevs(flags.frame, flags.unit) = ...
        detection.target_range * ...
        sind(multi.angles(flags.frame, flags.unit) + ...
        multi.gimbal_angles(flags.frame, flags.unit));
    
end

end

