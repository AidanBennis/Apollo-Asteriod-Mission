clc; clear; close all;

%% Constants
mu = 132712440018; % km^3/s^2
AU = 149597885.092797; % Astronomical Unit in km
G = 6.67430e-11; % Gravitational constant [m^3 kg^-1 s^-2]


%% Orbital parameters
% Earth
a1 = 1 * AU; % Semi-major axis of Earth's orbit in km
e1 = 0; % Eccentricity of Earth's orbit
i1 = 0; % Inclination
r2 = 42157; % Orbit around Earth in Km (6371 [radius] + 35786 [GTO])
muT = 3.99e+5 ;% gravatational parameter of Earth

% Asteroid: Apollo
a3 = 2.2e8; % Semi-major axis of Apollos orbit in km
e3 = 0.559871597; % Eccentricity of Apollo orbit
i3 = 6.352454347; % Inclination of Apollo orbit in degrees
Omega3 = 35.55648326; % Right Ascension of Apollo orbit in degrees
omega3 = 286.0329799; % Argument of periapsis of Apollo orbit in degrees
M3 = 41.059574; % Mean anomaly of Apollo orbit in degrees
rL = 5.775; % raduis of intial orbit around launch body (including raduis of body) (km)
rT = 42157; % raduis of capture orbit around target body (including raduis of body) (km)
rT1 = 6571; % 200Km LEO Orbit for end of mission
muL = 1.33e-7;% gravatational parameter of Apollo

% Sun
muS = 1.33e+11; % gravatational parameter of sun

%% Orbit Definition
% Orbital velocity for circular orbit
v_orbit = sqrt(muL*1000 / rL*1000);

% Orbital period
T = 2 * pi * sqrt(rL^3 / muL);

% Output results
fprintf('Altitude: %.2f Km\n', rL - 0.775);
fprintf('Orbital Velocity: %.2f m/s\n', v_orbit);
fprintf('Orbital Period: %.2f minutes', T / 60);
fprintf(' (%.2f hours)\n\n', T / 60 / 60);
%% Convert orbital elements position and velocity
% Cartesian --> elements
orb = [a3; e3; i3/180*pi; Omega3/180*pi; omega3/180*pi; M3/180*pi];
rv = E2C(orb, mu); % E2C is a function that converts elements to cartesian

% elements --> cartesian
orb_elements = C2E(rv, mu); % C2E is a function that converts cartesian to elements

%% Hohmann transfer
% Calculating periapsis and apoapsis for Apollo
rp = a3 * (1 - e3); % Periapsis distance of Apollo orbit
ra = a3 * (1 + e3); % Apoapsis distance of Apollo orbit

% Semi-major axis of the Hohmann transfer orbit
a2 = (rp + ra) / 2;

%% Velocity increments for Hohmann transfer Return Leg
delta_v1r = abs(sqrt((2*muS/a3)-(2*muS/(a1+a3))) - sqrt(muS/a3)); % Calculating deltav for first Hohmann transfer manouver
deltavr_escape = abs(sqrt((2*muL/rL)+(delta_v1r^2)) - sqrt((2*muL/rL)-(muL/a3))); % Calculating deltav for launch body escape
delta_v2r = abs(sqrt(muS/a1) - sqrt((2*muS/a1) - (2*muS/(a1+a3)))); % Calculating deltav for second Hohmann transfer manouver
deltavr_capture = abs(sqrt(((2*muT)/rT1) + (delta_v2r^2)) - sqrt(((2*muT)/rT1)-(muT/a1))); % Calculating deltav for target body capture 
Totaldvr = abs(deltavr_escape + delta_v2r + deltavr_capture); % Total deltav of the transfer

% Display results
fprintf('Delta v1r (km/s): %.4f\n', delta_v1r);
fprintf('Delta v2r (km/s): %.4f\n', delta_v2r);
fprintf('Delta vr Total (km/s): %.4f\n', Totaldvr);

% Propellant mass 
Isp = 450; % Specific impulse of the rocket engine in seconds
g = 9.81; % Standard gravity in m/s^2
m2 = 2000; % Spacecraft mass - with cargo in Kg
ve = (Isp*g)/1000; % Exhaust Velocity in Km/s
mi1 = m2/(exp(-Totaldvr/ve)); % Total Mass of Spacecraft
mfuel1 = mi1 - m2;
fprintf('Propellant Mass Return (kg): %.4f\n\n', mfuel1);

%% Velocity increments for Hohmann transfer 
delta_v1 = abs(sqrt((2*muS/a1)-(2*muS/(a1+a3))) - sqrt(muS/a1)); % Calculating deltav for first Hohmann transfer manouver
deltav_escape = abs(sqrt((2*muT/rT)+(delta_v1^2)) - sqrt((muT/rT))); % Calculating deltav for launch body escape
delta_v2 = abs(sqrt(muS/a3) - sqrt((2*muS/a3) - (2*muS/(a1+a3)))); % Calculating deltav for second Hohmann transfer manouver
deltav_capture = abs(sqrt(((2*muL)/rL) + (delta_v2^2)) - sqrt(((2*muL)/rL)-(muL/a3))); % Calculating deltav for target body capture 
Totaldv = abs(deltav_escape + delta_v2 + deltav_capture); % Total deltav of the transfer

% Display results
fprintf('Delta v1 (km/s): %.4f\n', delta_v1);
fprintf('Delta v2 (km/s): %.4f\n', delta_v2);
fprintf('Delta v Total (km/s): %.4f\n', Totaldv);

% Propellant mass
m3 = 1500; % Spacecraft mass - no cargo in Kg
mi2 = mi1/(exp(-Totaldv/ve)); % Earth to Apollo Total Mass
mfuel2 = mi2 - m3; % Earth to Apollo Fuel Mass
fprintf('Propellant Mass Outgoing (kg): %.4f\n\n', mfuel2);

%Total Mass
mi_launch = mi1 + mi2;
mfuel_launch = mfuel1 + mfuel2;
fprintf('Total Launch Mass (kg): %.4f\n', mi_launch);
fprintf('Total Fuel Mass (kg): %.4f\n', mfuel_launch);

% Adjust periapsis distance for Apollo to a realistic value
close_approach_distance = 0.35 * AU; % Adjusted periapsis distance of Apollo during close approach

%% Hohmann transfer orbit for close approach
% Apollo's periapsis during close approach
rp = close_approach_distance; % Periapsis distance of Apollo at close approach

% Earth's position will remain at its circular orbit (a1)
% Semi-major axis of the Hohmann transfer orbit
a2 = (a1 + rp) / 2; % Semi-major axis for the Hohmann transfer

% Eccentricity of the transfer orbit
e2 = (rp - a1) / (rp + a1);

% Transfer time (half the period of the transfer orbit)
T_transfer = pi * sqrt(a2^3 / mu); % Time for half of the Hohmann transfer

%% Earth's position at launch
% Earth's mean motion
n_Earth = sqrt(mu / a1^3); % rad/s

% True anomaly of Earth at launch
theta_Earth_launch = n_Earth * T_transfer; % Angle by which Earth has moved during the transfer time

% X and Y coordinates of Earth at launch
x_Earth_launch = a1 * cos(theta_Earth_launch);
y_Earth_launch = a1 * sin(theta_Earth_launch);

%% Solving for True Anomaly at Close Approach for Apollo
n_Apollo = sqrt(mu / a3^3); % Mean motion of Apollo
M_Apollo_capture = M3 + n_Apollo * T_transfer; % Mean anomaly after transfer time

% Solve Kepler's equation for eccentric anomaly (E)
E_Apollo_capture = M_Apollo_capture; % Initial guess for E
for iter = 1:1000 % Iterate to solve for E
    E_Apollo_capture = M_Apollo_capture + e3 * sin(E_Apollo_capture); 
end

% True anomaly (nu) of Apollo at capture
nu_Apollo_capture = 2 * atan2(sqrt(1 + e3) * sin(E_Apollo_capture / 2), sqrt(1 - e3) * cos(E_Apollo_capture / 2));

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

% Earth's orbit (or the initial orbit)
x_Earth = a1 * cos(theta); % X-coordinates for Earth's orbit
y_Earth = a1 * sin(theta); % Y-coordinates for Earth's orbit
plot(x_Earth, y_Earth, 'b', 'LineWidth', 1.5); % Plot Earth's orbit
legend_entries{end+1} = 'Earth Orbit';

% Apollo orbit (or the target orbit)
x_Apollo = a3 * (cos(theta) - e3); % X-coordinates for Apollo orbit (correcting for center offset)
y_Apollo = a3 * sqrt(1 - e3^2) * sin(theta); % Y-coordinates for Apollo orbit
plot(x_Apollo, y_Apollo, 'r', 'LineWidth', 1.5); % Plot Apollo orbit
legend_entries{end+1} = 'Apollo Orbit';

% Plot the rotated transfer orbit
plot(x_transfer_rotated, y_transfer_rotated, 'g--', 'LineWidth', 1.5); % Plot transfer orbit
legend_entries{end+1} = 'Transfer Orbit';

% Plot the positions of the bodies at launch and capture
plot(x_Earth_launch, y_Earth_launch, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Earth at launch
plot(a3 * (cos(nu_Apollo_capture) - e3), a3 * sqrt(1 - e3^2) * sin(nu_Apollo_capture), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Apollo at capture
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
