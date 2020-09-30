%% TRASAT Radar System - Example Scenario Initialization File
%{

    Sean Holloway
    TRASAT Init File
    
    This file specifies all parameters of a radar simulation test for the
    TRASAT project. 

    Use script 'FullSystem_AutomatedSimulation.m' to run scenarios.
    
%}

%% Housekeeping

% Initialize class
scenario = RadarScenario;
scenario.flags.sim_rate = sim_rate;

%% Multistatic Scenario Setup

% Number of frames and receivers to simulate
scenario.multi.n_fr = 120;
scenario.multi.n_re = 9;

% Allocate data structures
multiSetup(scenario);

% Locations of radar units
scenario.multi.radar_pos = ...
    [3 * nm * ones(1,scenario.multi.n_re)...
                .*[-1 1 -1 1 -1 1 -1 1 -1]; ...               % Alternating x location
    1.5 * nm * ((0:(scenario.multi.n_re-1))-1); ...       % Incremental y distance
    0 * nm * ones(1,scenario.multi.n_re)];            % Constant z elevation

% Normal direction of radar units
scenario.multi.radar_dir = [0 180 0 180 0 180 0 180 0; ...
                            0 0 0 0 0 0 0 0 0];

% Multistatic processing method
%   Method of deciding which units to use for multilateration each frame
scenario.multi.method = 'SNR';                  % 'SNR' or 'Range'
scenario.multi.track_method = 'Kalman Filter'; % 'Moving Average' or 'None'

% Tracking Parameters
scenario.multi.sigma_v = [5, 0.001, 150];
scenario.multi.sigma_z = [1, 0.5, 1100];
scenario.multi.interp = true;
scenario.multi.smooth = true;

%% Target RCS Setup

% Options setup
scenario.rcs = struct( ...
    ...
    ... % RCS options
    'rcs_model',    'model', ...        % Set 'model' or 'constant'
    'ave_rcs',      -10, ...            % Target RCS in dBm^2
    ...
    ... % Model options
    'dim',      [10; 10; 0], ...          % [x; y; z] size of target
    'n_sc',     50, ...                 % Number of point scatterers
    'res_a',    0.1, ...                % Angle resolution in degrees
    'freq',     9e9:10e6:10e9);    % Frequency range to model

% Run RCS model
scenario.rcs = TargetRCSModel(scenario.rcs);

% View RCS results
% viewRCSFreq(scenario);
% viewRCSAng(scenario);

%% Target Trajectory Setup
% 
% Options setup
scenario.traj = struct( ...
    ...
    ... % Trajectory options
    'alt',      2000, ...                % Altitude in meters
    'yvel',     200, ...                % Along track velocity in m/s
    'exc',      2000, ...                % Excursion distance in meters
    'per',      0.1, ...                % Excursion period (Nominally 0.05 to 0.2)
    ...
    'pos_st',   [100; 0; 0], ...          % Position input if 'static' is used
    'vel_st',   [200; 0; 0], ...          % Velocity input if 'static' is used
    ...
    ... % Model options
    'model',    'model', ...               % Set 'static' or 'model'
    'time',     0 : 100e-6 : 6.2*10.24);     % Time of simulation in seconds
% NOTE: Set time step to PRI

% Run Trajectory model
scenario.traj = TrajectoryModel(scenario.traj);

% View Trajectory
% viewTraj(scenario);

%% Simulation Setup

% Radar simulation and processing setup
scenario.simsetup = struct( ...
    ...
    ... % Waveform Properties
    'f_c',      9.525e9, ...            % Operating frequency in Hz
    'f_s',      50e6, ...               % ADC sample frequency in Hz
    't_p',      20e-6, ...             % Pulse duration in seconds
    'n_ch',     13, ...                 % Number of chips per pulse
    'prf',      10e3, ...               % Pulse repetition frequency in Hz
    'n_p',      512, ...                % Number of pulses per CPI
    'cpi_fr',   10, ...                 % Number of CPI per frame
    ...
    ... % Transceiver Properties
    'n_ant',        16, ...             % Number of elements in antenna array
    'd_ant',        5, ...              % Vertical spacing between antenna arrays
    'tx_pow',       64, ...             % Transmit power in Watts
    'tx_ant_gain',  0, ...             % Tx antenna gain in dBi (N/A?)
    'rx_sys_gain',  -4 + 3, ...              % Rx system gain in dB (N/A?)
    'rx_nf',        2.7, ...            % Rx noise figure in dB
    'rx_ant_gain',  0, ...             % Rx antenna gain in dBi (N/A?)
    ...
    ... % Processing Properties
    'win_type',     'blackmanharris', ...   % Type of window for doppler processing
    'mono_res',     0.01, ...               % Resolution of monopulse histogram
    'mono_smooth',  25, ...                 % Smoothing window size for monopulse histogram
    ...
    ... % Detection Properties
    'thresh',       9, ...             % Detection threshold above noise power in dB
    'det_m',        2, ...              % M for m-of-n processing
...
    ... % Tracking Properties
    'wind_size',    [3; 3; 3]);         % Size of averaging window in [x;y;z] direction 
    

% Set up Phased Array Toolbox system objects
scenario.sim = PhasedSetup_TRASAT(scenario);

