%% TRASAT Radar System
%{

    Sean Holloway
    TRASAT (Tracking Radar for Single Airborne Target) System
    MATLAB Simulation & Processing

    This shell file runs successive scripts and gauges progress.

    TODOs in AutomatedSimulation
    
%}


%% Housekeeping
clear variables
close all
addpath(genpath(pwd));
tic

profile on

%% Definitions

nm = 1852;                      % Nautical miles in meters
c = physconst('LightSpeed');    % Speed of light in m/s

%% User Options

% Filename for saving radar scenario object
filename = 'TrajectoryScenario_TRASAT';

% Format to save figure files
save_format = {'.png', '.fig'};

% Rate to divide up fast time x slow time processing rates
sim_rate = 2^9;
% Optimized at 2^9 for this project

%% Setup Radar Scenario

% Scenario file to run
ExampleScenario_TRASAT

%% Run Simulation

% Perform main loop of simulation
MainLoop_TRASAT


%% Save Figures and Data
%{
% Establish file name
save_name = [filename, '_', datestr(now, 'mmddyy_HHMM')];

mat_path = 'MAT Files\Scenario Objects\';
fig_path = ['Figures\', save_name, '\'];

% Save scenario object
SaveScenario(scenario, save_name, mat_path);

% Save open figures
for ftype = 1:length(save_format)
    SaveFigures(save_name, fig_path, save_format{ftype});
end

% Read elapsed time
toc
%}

profile viewer









