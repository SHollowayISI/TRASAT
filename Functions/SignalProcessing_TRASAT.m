function [cube] = SignalProcessing_TRASAT(scenario)
%SIGNALPROCESSING_TRASAT Performs signal processing for TRASAT
%   Takes scenario struct as input, retuns scenario.cube struct containing
%   processed Range-Doppler cube

%% Unpack Variables

sim = scenario.sim;
simsetup = scenario.simsetup;

%% Define Constants

c = physconst('LightSpeed');
lambda = c/simsetup.f_c;

%% Perform Range Matched Filter

% Apply matched filter to each channel
% cube.range_cube = sim.match(sum(scenario.rx_sig,3));
% cube.range_cube = permute(sim.match(scenario.rx_sig), [1 2 4 3]);
cube.range_cube = sim.match(scenario.rx_sig(:,:,1));

% Split into fast time x slow time x CPI x antenna
% cube.range_cube = reshape(cube.range_cube, ...
%     size(cube.range_cube,1), simsetup.n_p, simsetup.cpi_fr, 2);

cube.range_cube = reshape(cube.range_cube, ...
    size(cube.range_cube,1), simsetup.n_p, simsetup.cpi_fr);

%% Perform Doppler FFT

% Calculate FFT size
N_d = 2^ceil(log2(size(cube.range_cube,2)));

% Apply windowing
expression = '(size(cube.range_cube,2))).*cube.range_cube;';
expression = ['transpose(', simsetup.win_type, expression];
cube.rd_cube = eval(expression);

% FFT across slow time dimension
cube.rd_cube = fftshift(fft(cube.rd_cube, N_d, 2), 2);

% Wrap max negative frequency and positive frequency
cube.rd_cube = [cube.rd_cube, cube.rd_cube(:,1,:)];

%% Calculate Square Power Cube

% Take square magnitude of radar cube, integrated between antennas?
cube.pow_cube_CPI = abs(cube.rd_cube).^2;

% Average across CPI
cube.pow_cube = mean(cube.pow_cube_CPI, 3);

%% Derive Axes

% Derive Range axis
cube.range_offset = (c/2)*simsetup.t_p;
cube.range_res = (c/2)/simsetup.f_s;
cube.range_axis = (1:(simsetup.f_s/simsetup.prf))*cube.range_res - cube.range_offset;

% Derive Doppler axis
cube.vel_res = lambda*(simsetup.prf)/(2*simsetup.n_p);
cube.vel_axis = -cube.vel_res*((-N_d/2):(N_d/2));

end

