% P = 64;
% G = db2pow(23);
% t_p = 20e-6;
% f_c = 9.525e9;
% lambda = physconst('LightSpeed')/f_c;
% rcs = db2pow(0);
% n_p = 512;
% k = physconst('Boltzmann');
% T = 290;
% NF = db2pow(2.7);
% L_r = db2pow(4);
% 
% R = 1;
% 
% SNR = (P*G*G*t_p*lambda*lambda*rcs*n_p)./(((4*pi)^3)*(R.^4)*k*T*NF*L_r);
% 
% Sth = 10;
% S = SNR/Sth;
% Rmax = S^(1/4)
% Rmax_nmi = Rmax/1852
% 
% bw = 50;
% 
% N = -log(2)/log(cosd(bw));
% 
% xmax = Rmax*((N/(N+2))^(N/4))*sqrt(2/(N+2))
% xmax_nmi = xmax/1852
% 
% y_at_xmax = Rmax*((N/(N+2))^((N+2)/4))
% y_at_xmax_nmi = y_at_xmax/1852
% 
% phi = -(pi/2):0.0001:(pi/2);
% rho = Rmax*(cos(phi).^(N/2));
% 
% polar(phi, rho)


%% Accuracy vs. Range Plot
%{
load('MAT Files\Scenario Objects\scenario_T2_RangeTri_1kmAlt_TRASAT_031520_1820.mat')

elevs_test = scenario.multi.elevs(:,1:18);
elevs_diff = elevs_test - 1000;
elevs_error_1km = rms(reshape(rms(elevs_diff,1),3,[]),1)'
tilt_error = sqrt(((1:6)'*1852).^2 + 1000^2)*sind(0.1);
elevs_error_1km_tilt = sqrt(tilt_error.^2 + elevs_error_1km.^2)

trilat_z_test = squeeze(scenario.multi.lat_points_tri(3,:,1:6));
trilat_z_diff = trilat_z_test - 1000;
trilat_z_error_1km = rms(trilat_z_diff,1)'

load('MAT Files\Scenario Objects\scenario_T2_RangeTri_6kmAlt_TRASAT_031520_2016.mat')

elevs_test = scenario.multi.elevs(:,1:15);
elevs_diff = elevs_test - 6000;
elevs_error_6km = rms(reshape(rms(elevs_diff,1),3,[]),1)'
tilt_error = sqrt(((1:5)'*1852).^2 + 6000^2)*sind(0.1);
elevs_error_6km_tilt = sqrt(tilt_error.^2 + elevs_error_6km.^2)

trilat_z_test = squeeze(scenario.multi.lat_points_tri(3,:,1:5));
trilat_z_diff = trilat_z_test - 6000;
trilat_z_error_6km = rms(trilat_z_diff,1)'

figure;
plot(1:6, elevs_error_1km_tilt, 1:6, trilat_z_error_1km, ...
     1:5, elevs_error_6km_tilt, 1:5, trilat_z_error_6km)
grid on
title('Elevation Error vs. Ground Range')
xlabel('Ground Range [nmi]', 'FontWeight', 'bold')
ylabel('RMS Elevation Error [m]', 'FontWeight', 'bold')
legend('Monopulse Error - 3nmi Ground Range', 'Trilateration Error - 3nmi Ground Range', ...
    'Monopulse Error - 5nmi Ground Range', 'Trilateration Error - 5nmi Ground Range');




%}

%% Accuracy vs. Elevation Plot
%{
load('MAT Files\Scenario Objects\scenario_T4_AltTri_3nmiR_TRASAT_031620_1526.mat')

real_elev = 200*(0:4);

elevs_test = scenario.multi.elevs;
elevs_diff = elevs_test - reshape(repmat(real_elev, 3, 1), 1, []);
elevs_error_3nmi = rms(reshape(rms(elevs_diff,1),3,[]),1)'
tilt_error = sqrt((3*1852).^2 + (200*(0:4)').^2)*sind(0.1);
elevs_error_3nmi_tilt = sqrt(tilt_error.^2 + elevs_error_3nmi.^2)

trilat_z_test = squeeze(scenario.multi.lat_points_tri(3,:,:));
trilat_z_diff = trilat_z_test;
trilat_z_error_3nmi = rms(trilat_z_diff,1)'


load('MAT Files\Scenario Objects\scenario_T4_AltTri_5nmiR_TRASAT_031620_1638.mat')

real_elev = 200*(0:4);

elevs_test = scenario.multi.elevs;
elevs_diff = elevs_test - reshape(repmat(real_elev, 3, 1), 1, []);
elevs_error_5nmi = rms(reshape(rms(elevs_diff,1),3,[]),1)'
tilt_error = sqrt((5*1852).^2 + (200*(0:4)').^2)*sind(0.1);
elevs_error_5nmi_tilt = sqrt(tilt_error.^2 + elevs_error_3nmi.^2)

trilat_z_test = squeeze(scenario.multi.lat_points_tri(3,:,:));
trilat_z_diff = trilat_z_test;
trilat_z_error_5nmi = rms(trilat_z_diff,1)'

figure;
plot(0:0.2:0.8, elevs_error_3nmi_tilt, 0:0.2:0.8, trilat_z_error_3nmi, ...
     0:0.2:0.8, elevs_error_5nmi_tilt, 0:0.2:0.8, trilat_z_error_5nmi)
grid on
title('Elevation Error vs. Target Elevation')
xlabel('Target Elevation [km]', 'FontWeight', 'bold')
ylabel('RMS Elevation Error [m]', 'FontWeight', 'bold')
legend('Monopulse Error - 3nmi Ground Range', 'Trilateration Error - 3nmi Ground Range', ...
    'Monopulse Error - 5nmi Ground Range', 'Trilateration Error - 5nmi Ground Range');
%}

%% Feasible Trajectory Plotting
%{
close all
figure;

for n = -1:1
    traj_a_x = [-1389, -1389, 1389, 1389, -1389] - n*1852*1.5;
    traj_a_y = [0, 8000, 7500, 1500, 0];
    fill(traj_a_x, traj_a_y, [0.5 0.75 1], 'EdgeColor', [0 0 0.5])
    hold on
    
    traj_b_x = -traj_a_x;
    traj_b_y = traj_a_y;
    fill(traj_b_x, traj_b_y, [0.5 0.75 1], 'EdgeColor', [0 0 0.5])
    hold on
end

for n = -1:1
    traj_both_x = [0, 1389, 1389, 0, -1389, -1389, 0] - n*1852*1.5;
    traj_both_y = [750, 1500, 7500, 7750, 7500, 1500, 750];
    fill(traj_both_x, traj_both_y, [0 0.4470 0.7410], 'EdgeColor', [0 0 0.5])
    hold on
end

radar_loc_x = [-1389, 1389];
radar_loc_y = [0, 0];
scatter(radar_loc_x, radar_loc_y, 'red', 'filled')
text(0, 4000, 'Multiple Detection Zone', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
text(1500, 800, 'Single Detection Zone', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
text(-1500, 800, 'Single Detection Zone', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')

grid on
xlim([-4000 4000])
ylim([0 8000])
title('SEMTA Viable Multilateration Area')
xlabel('Distance Along Track [m]', 'FontWeight', 'bold')
ylabel('Distance From Radar Line [m]', 'FontWeight', 'bold')
%}


%% Fitting for Tracking Parameters

close all;

% load('C:\Users\sholloway\Documents\MATLAB\TRASAT\MAT Files\Scenario Objects\scenario_T2_Excursion_2000m_6nmi_n10dB_Repeat_TRASAT_032420_1909.mat')
load('C:\Users\sholloway\Documents\MATLAB\TRASAT\MAT Files\Scenario Objects\scenario_T2_Excursion_2000m_Repeat_TRASAT_031820_0807.mat')

miss_ind = (abs(scenario.multi.lat_points(3,:)-2000) > 500);
scenario.multi.lat_points(:,miss_ind) = nan;

multi = scenario.multi;
simsetup = scenario.simsetup;
traj = scenario.traj;


% Interpolate through NaNs
index = 1:scenario.multi.n_fr;
det_index = index(~isnan(scenario.multi.lat_points(1,:)));
miss_index = index(isnan(scenario.multi.lat_points(1,:)));
diff_index = diff(det_index);

int_type = {'spline', 'spline', 'cubic'};
for dim = 1:3
    scenario.multi.lat_points(dim,miss_index) = ...
        interp1(det_index, scenario.multi.lat_points(dim,det_index), miss_index, int_type{dim});
    
    det_index = index;
    diff_index = diff(det_index);
end



traj_ind = round(scenario.simsetup.n_p * scenario.simsetup.cpi_fr * ...
    (0.5:(scenario.multi.n_fr-0.5)) + 1);
traj_points = scenario.traj.pos(:,traj_ind);

% Values from 3-19-20:
% scenario.multi.sigma_v = [5, 0.001, 150];
% scenario.multi.sigma_z = [1, 0.5, 1100];
% pre_wind = [1, 1, 1];
% post_wind = [5, 5, 31];

scenario.multi.sigma_v = [5, 0.001, 150];
scenario.multi.sigma_z = [1, 0.1, 1100];
pre_wind = [1, 1, 1];
post_wind = [5, 5, 31];

start_ind = 5;
stop_ind = 5;

traj_points = traj_points(:,start_ind:(end - stop_ind + 1));
scenario.multi.lat_points = scenario.multi.lat_points(:,start_ind:(end - stop_ind + 1));
scenario.multi.n_fr = size(scenario.multi.lat_points,2);

for n = 1:3
    scenario.multi.lat_points(n,:) = ...
        movmean(scenario.multi.lat_points(n,:), pre_wind(n));
end

scenario.multi = TargetTracking_TRASAT(scenario);

for n = 1:3
    scenario.multi.track_points(n,:) = ...
        movmean(scenario.multi.track_points(n,:), post_wind(n));
end

track_error = scenario.multi.track_points - traj_points;

post_start = 5;
post_end = 5;

track_error = track_error(:,post_start:(end-post_end+1));

track_error_rms = rms(track_error,2)
rms(track_error_rms)

[~, index] = sort(scenario.multi.SNR, 2, 'descend');
elevs = scenario.multi.elevs(index(:,1));
elevs = elevs((start_ind):(end - stop_ind + 1),:);
% elevs = movmean(elevs, 1);
% elevs_error = elevs - traj_points(3,:)';
for n = 1:scenario.multi.n_fr
    elevs(n) = scenario.multi.elevs(n,index(n,1));
    ranges(n) = scenario.multi.ranges(n,index(n,1));
end

tilt_error = ranges'.*sind(0.1*(rand(size(elevs))-0.5));

elevs = elevs + tilt_error;

elevs = movmean(elevs, 31);

elevs_error = elevs - traj_points(3,:)';

elevs_error = elevs_error((start_ind):(end - stop_ind + 1),:);

rms(elevs_error)
rms([track_error_rms(1:2); rms(elevs_error)])

figure;
plot(track_error')
hold on
plot(elevs_error)
grid on
xlim([0 100])
ylim([-100 100])
xlabel('Time [Radar Frame]', 'FontWeight', 'bold')
ylabel('Error [m]', 'FontWeight', 'bold')
legend('Cross-Track Direction', 'Along-Track Direction', 'Multilateration Altitude', 'Monopulse Altitude')

figure;
scatter3(scenario.multi.lat_points(1,:)/1000, scenario.multi.lat_points(2,:)/1000, scenario.multi.lat_points(3,:)/1000)
hold on
plot3(scenario.multi.track_points(1,:)/1000, scenario.multi.track_points(2,:)/1000, scenario.multi.track_points(3,:)/1000)
hold on
scatter3(scenario.multi.radar_pos(1,:)/1000, scenario.multi.radar_pos(2,:)/1000, scenario.multi.radar_pos(3,:)/1000, 'red', 'filled')
xlabel('Distance from Center [km]', 'FontWeight', 'bold')
ylabel('Distance Along Track [km]', 'FontWeight', 'bold')
zlabel('Altitude [km]', 'FontWeight', 'bold')
zlim([0 10])
legend('Multilateration Points', 'Tracking Result', 'Radar Unit Locations')
title('3D Scatter Plot - Multilateration')

figure;
scatter3(scenario.multi.lat_points(1,:)/1000, scenario.multi.lat_points(2,:)/1000, elevs/1000)
hold on
plot3(scenario.multi.track_points(1,:)/1000, scenario.multi.track_points(2,:)/1000, scenario.multi.track_points(3,:)/1000)
hold on
scatter3(scenario.multi.radar_pos(1,:)/1000, scenario.multi.radar_pos(2,:)/1000, scenario.multi.radar_pos(3,:)/1000, 'red', 'filled')
xlabel('Distance from Center [km]', 'FontWeight', 'bold')
ylabel('Distance Along Track [km]', 'FontWeight', 'bold')
zlabel('Altitude [km]', 'FontWeight', 'bold')
zlim([0 10])
legend('Multilateration Points', 'Tracking Result', 'Radar Unit Locations')
title('3D Scatter Plot - Multilateration w/ Monopulse')

%% Monopulse Error Investigation
%{
[~, index] = sort(scenario.multi.SNR, 2, 'descend');
elevs = scenario.multi.elevs(index(:,1));

traj_ind = round(scenario.simsetup.n_p * scenario.simsetup.cpi_fr * ...
    (0.5:(scenario.multi.n_fr-0.5)) + 1);
traj_points = scenario.traj.pos(:,traj_ind);

for n = 1:scenario.multi.n_fr
    [true_ranges(n), true_angles(:,n)] = rangeangle( ...
        traj_points(:,n), scenario.multi.radar_pos(:,index(n,1)), ...
        azelaxes( ...
        scenario.multi.radar_dir(1,index(n,1)), ...
        scenario.multi.radar_dir(2,index(n,1))));
end

true_elevs = true_ranges .* sind(true_angles(2,:));

for n = 1:scenario.multi.n_fr
    angle_error_chosen(n) = scenario.results.angle_error(n,index(n,1));
    elev_error_chosen(n) = scenario.multi.elevs(n,index(n,1)) - traj_points(3,n);
end

% plot(angle_error_chosen)
figure
plot(elev_error_chosen)

% elevs = elevs.*cosd(true_angles(1,:)');

% plot(elevs)

% mean(elevs)
%}





% 
% 
% hold on
% t = 0.512*(0:1:120);
% 
% x = 2*sin(0.1*t);
% y = 0.2*t;
% z = 6*ones(size(y));
% 
% plot3(x', y', z');
% legend('Multilateration Points', 'Tracking Result', ...
%     'Radar Unit Locations', 'Real Trajectory')


















