%clc; clear; close all;

%% Constants
muS = 132712440018; % Gravitational Parameter of the Sun [km^3/s^2]
AU = 149597885.092797; % Astronomical Unit in km
G = 6.67430e-11; % Gravitational constant [m^3 kg^-1 s^-2]

%% Orbital Parameters
% Earth
a1 = 1 * AU; % Semi-major axis of Earth's orbit in km
e1 = 0; % Eccentricity of Earth's orbit
i1 = 0; % Inclination
muE = 3.99e+5 ; % gravatational parameter of Earth
rE = 42157; % Orbit around Earth in Km (6371 [radius] + 35786 [GTO])
rE1 = 6571; % 200Km LEO Orbit for end of mission

% Asteroid: Apollo
a4 = 2.2e8; % Semi-major axis of Apollos orbit in km
e4 = 0.559871597; % Eccentricity of Apollo orbit
i4 = 6.352454347; % Inclination of Apollo orbit in degrees
Omega4 = 35.55648326; % Right Ascension of Apollo orbit in degrees
omega4 = 286.0329799; % Argument of periapsis of Apollo orbit in degrees
M4 = 41.059574; % Mean anomaly of Apollo orbit in degrees
rA = 5.775; % raduis of intial orbit around Apollo (Km)
muA = 1.33e-7; % Gravatational parameter of Apollo

%% Orbit Definition
% Orbital velocity for circular orbit
v_orbit = sqrt(2*muA*1000 / rA*1000);

% Orbital period
T = 2 * pi * sqrt(rA^3 / muA);

% Output results
fprintf('Altitude: %.2f Km\n', rA - 0.775);
fprintf('Orbital Velocity: %.2f m/s\n', v_orbit);
fprintf('Orbital Period: %.2f minutes', T / 60);
fprintf(' (%.2f hours)\n\n', T / 60 / 60);

%% Convert orbital elements position and velocity
% Cartesian --> elements
orb = [a4; e4; i4/180*pi; Omega4/180*pi; omega4/180*pi; M4/180*pi];
rv = E2C(orb, muS); % E2C converts elements to cartesian

% elements --> cartesian
orb_elements = C2E(rv, muS); % C2E converts cartesian to elements

%% Velocity increments for Hohmann transfer Return Leg
delta_v1r = abs(sqrt((2*muS/a4)-(2*muS/(a1+a4))) - sqrt(muS/a4)); % Calculating deltav for first Hohmann transfer manouver
deltavr_escape = abs(sqrt((2*muA/rA)+(delta_v1r^2)) - sqrt((2*muA/rA)-(muA/a4))); % Calculating deltav for Apollo escape
delta_v2r = abs(sqrt(muS/a1) - sqrt((2*muS/a1) - (2*muS/(a1+a4)))); % Calculating deltav for second Hohmann transfer manouver
deltavr_capture = abs(sqrt(((2*muE)/rE1) + (delta_v2r^2)) - sqrt(((2*muE)/rE1)-(muE/a1))); % Calculating deltav for Earth capture 
Totaldvr = abs(deltavr_escape + delta_v2r + deltavr_capture); % Total deltav of the transfer

% Display results
fprintf('Delta v1r (km/s): %.4f\n', delta_v1r);
fprintf('Delta v2r (km/s): %.4f\n', delta_v2r);
fprintf('Delta vr Total (km/s): %.4f\n', Totaldvr);

% Propellant mass 
Isp = 461; % Specific impulse of the rocket engine in seconds
g = 9.81; % Standard gravity in m/s^2
m2 = 2232; % Spacecraft mass + RL10 Engine - with cargo in Kg
ve = (Isp*g)/1000; % Exhaust Velocity in Km/s
mi1 = m2/(exp(-Totaldvr/ve)); % Total Mass of Spacecraft (spacecraft + fuel)
mfuel1 = mi1 - m2; % Mass of Fuel for Return Leg
fprintf('Propellant Mass Return (kg): %.4f\n\n', mfuel1);

%% Velocity increments for Hohmann transfer 
delta_v1 = abs(sqrt((2*muS/a1)-(2*muS/(a1+a4))) - sqrt(muS/a1)); % Calculating deltav for first Hohmann transfer manouver
deltav_escape = abs(sqrt((2*muE/rE)+(delta_v1^2)) - sqrt((muE/rE))); % Calculating deltav for Earth escape
delta_v2 = abs(sqrt(muS/a4) - sqrt((2*muS/a4) - (2*muS/(a1+a4)))); % Calculating deltav for second Hohmann transfer manouver
deltav_capture = abs(sqrt(((2*muA)/rA) + (delta_v2^2)) - sqrt(((2*muA)/rA)-(muA/a4))); % Calculating deltav for Apollo capture 
Totaldv = abs(deltav_escape + delta_v2 + deltav_capture); % Total deltav of the transfer

% Display results
fprintf('Delta v1 (km/s): %.4f\n', delta_v1);
fprintf('Delta v2 (km/s): %.4f\n', delta_v2);
fprintf('Delta v Total (km/s): %.4f\n', Totaldv);

% Propellant mass
m3 = 1732; % Spacecraft mass + RL10 Engine - no cargo in Kg
mi2 = mi1/(exp(-Totaldv/ve)); % Earth to Apollo Total Mass
mfuel2 = mi2 - m3; % Earth to Apollo Fuel Mass
fprintf('Propellant Mass Outgoing (kg): %.4f\n\n', mfuel2);

%Total Mass
mfuel_launch = mfuel1 + mfuel2;
mi_launch = mfuel_launch + m3;
fprintf('Total Launch Mass (kg): %.4f\n', mi_launch);
fprintf('Total Fuel Mass (kg): %.4f\n\n', mfuel_launch);

% Adjust periapsis distance for Apollo to a realistic value
close_approach_distance = 0.35 * AU; % Adjusted periapsis distance of Apollo during close approach

%% Hohmann Transfer Orbit for Close Approach

rp = close_approach_distance; % Periapsis distance of Apollo at close approach
a2 = (a1 + rp) / 2; % Semi-major axis for the Hohmann transfer orbit
e2 = (rp - a1) / (rp + a1); % Eccentricity of the transfer orbit
T_transfer = pi * sqrt(a2^3 / muS); % Time for half of the Hohmann transfer
fprintf('Hohmann Transfer Time: %.4f seconds', T_transfer);
fprintf(' (%.2f days)\n\n', T_transfer / 60 / 60 / 24);

%% Earth's position at launch
% Earth's mean motion
n_Earth = sqrt(muS / a1^3); % rad/s

% True anomaly of Earth at launch
theta_Earth_launch = n_Earth * T_transfer; % Angle by which Earth has moved during the transfer time

% X and Y coordinates of Earth at launch
x_Earth_launch = a1 * cos(theta_Earth_launch);
y_Earth_launch = a1 * sin(theta_Earth_launch);

%% Solving for True Anomaly at Close Approach for Apollo
n_Apollo = sqrt(muS / a4^3); % Mean motion of Apollo
M_Apollo_capture = M4 + n_Apollo * T_transfer; % Mean anomaly after transfer time

% Solve Kepler's equation for eccentric anomaly (E)
E_Apollo_capture = M_Apollo_capture; % Initial guess for E
for iter = 1:1000 % Iterate to solve for E
    E_Apollo_capture = M_Apollo_capture + e4 * sin(E_Apollo_capture); 
end

% True anomaly (nu) of Apollo at capture
nu_Apollo_capture = 2 * atan2(sqrt(1 + e4) * sin(E_Apollo_capture / 2), sqrt(1 - e4) * cos(E_Apollo_capture / 2));

%% Rotating the Transfer Orbit
% Create a theta array for plotting the transfer orbit
theta_transfer = linspace(0, pi, 1000);  % For half of the transfer orbit (periapsis to apoapsis)

% Compute the transfer orbit's X and Y coordinates
x_transfer = a2 * (cos(theta_transfer) - e2); % X-coordinates for transfer orbit (shifted for eccentricity)
y_transfer = a2 * sqrt(1 - e2^2) * sin(theta_transfer); % Y-coordinates for transfer orbit

% Rotate the transfer orbit to match Earth's position at launch
x_transfer_rotated = x_transfer * cos(theta_Earth_launch) - y_transfer * sin(theta_Earth_launch);
y_transfer_rotated = x_transfer * sin(theta_Earth_launch) + y_transfer * cos(theta_Earth_launch);

%% Plot
% Create a theta array for plotting orbits
theta = linspace(0, 2*pi, 1000);         % For full orbits (Earth and Apollo)

% Plot orbits with corrections
figure;
hold on;
grid on;
axis equal;

% Plot the Sun at the origin
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y'); % Yellow Sun at (0, 0)
legend_entries = {'Sun'}; % Add Sun to legend

% Earth's Orbit
x_Earth = a1 * cos(theta); % X-coordinates for Earth's orbit
y_Earth = a1 * sin(theta); % Y-coordinates for Earth's orbit
plot(x_Earth, y_Earth, 'b', 'LineWidth', 1.5); % Plot Earth's orbit
legend_entries{end+1} = 'Earth Orbit';

% Apollo Orbit
x_Apollo = a4 * (cos(theta) - e4); % X-coordinates for Apollo orbit (correcting for center offset)
y_Apollo = a4 * sqrt(1 - e4^2) * sin(theta); % Y-coordinates for Apollo orbit
plot(x_Apollo, y_Apollo, 'r', 'LineWidth', 1.5); % Plot Apollo orbit
legend_entries{end+1} = 'Apollo Orbit';

% Plot the rotated transfer orbit
plot(x_transfer_rotated, y_transfer_rotated, 'g--', 'LineWidth', 1.5); % Plot transfer orbit
legend_entries{end+1} = 'Transfer Orbit';

% Plot the positions of the bodies at launch and capture
plot(x_Earth_launch, y_Earth_launch, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Earth at launch
plot(a4 * (cos(nu_Apollo_capture) - e4), a4 * sqrt(1 - e4^2) * sin(nu_Apollo_capture), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Apollo at capture
legend_entries{end+1} = 'Earth at Launch';
legend_entries{end+1} = 'Apollo at Capture';

% Set axis limits to improve visibility and keep the plot centered
xlim([-1.75 * 10^8, 1.75 * 10^8]); % Set X-axis limits based on the orbit size
ylim([-1.25 * 10^8, 1.25 * 10^8]); % Set Y-axis limits for better proportion

% Labels and title
xlabel('X position (km)');
ylabel('Y position (km)');
title('Hohmann Transfer for Apollo Close Approach in 2039');
legend(legend_entries, 'Location', 'best');
hold off;

%% Communication Window Analysis
% Time vector for the entire mission 
t_max = 365.25 * 24 * 3600; % 1 year in seconds
time_step = 3600; % 1 hour in seconds
vehicle_time = 0:time_step:t_max; % Time array from 0 to t_max with 1-hour intervals

% Earth Position
x_Earth = a1 * cos(2 * pi / t_max * vehicle_time);
y_Earth = a1 * sin(2 * pi / t_max * vehicle_time);
z_Earth = zeros(size(vehicle_time)); % Earth is assumed to be in the equatorial plane

% Apollo Position
x_Apollo = a4 * (cos(2 * pi / t_max * vehicle_time) - e4);
y_Apollo = a4 * sqrt(1 - e4^2) * sin(2 * pi / t_max * vehicle_time);
z_Apollo = zeros(size(vehicle_time)); % Apollo is also assumed to be in the equatorial plane

% Spacecraft Position (assuming spacecraft is initially at Earth and moves towards Apollo)
x_spacecraft = linspace(a1, a4, length(vehicle_time)); % Linearly move from Earth's orbit to Apollo's orbit
y_spacecraft = linspace(0, 0, length(vehicle_time)); % Keep in equatorial plane for simplicity
z_spacecraft = zeros(size(vehicle_time)); % No inclination change

% Combine Spacecraft Position Vectors into State Vector
vehicle_state_vector = [x_spacecraft', y_spacecraft', z_spacecraft'];

% Define Constants for Communication Function
constants.r_earth = 6371; % Earth's radius in km

% Call Communication Window Function
planet_state_vector = []; % Not used in current comms_window implementation

% Call the communication window function
comms_status = comms_window(vehicle_state_vector, vehicle_time, planet_state_vector, constants);

% Plot the Communication Status
figure;
plot(vehicle_time / 3600, comms_status, 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Communication Window Status');
title('Communication Window between Spacecraft and Ground Station');
ylim([0, 1]);  % Ensure y-axis limits are between 0 and 1 for clarity
yticks([0 1]); % Add y-axis tick marks for 0 and 1 to make binary status clear
grid on;
legend('0 = No Communication, 1 = Communication Open');
