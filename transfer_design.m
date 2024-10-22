clc; clear; close all;

% Constants
mu = 132712440018; % km^3/s^2
AU = 149597885.092797; % Astronomical Unit in km

% Orbital parameters
% Earth
a1 = 1 * AU; % Semi-major axis of Earth's orbit in km
e1 = 0; % Eccentricity of Earth's orbit
i1 = 0; % Inclination
r2 = 42157; % Orbit around Earth in Km (6371 [radius] + 35786 [GTO])
muT = 3.99e+5 ;% gravatational parameter of target body (earth)

% Asteroid: Apollo
a3 = 2.2e8; % Semi-major axis of Apollos orbit in km
e3 = 0.559871597; % Eccentricity of Apollo orbit
i3 = 6.352454347; % Inclination of Apollo orbit in degrees
Omega3 = 35.55648326; % Right Ascension of Apollo orbit in degrees
omega3 = 286.0329799; % Argument of periapsis of Apollo orbit in degrees
M3 = 41.059574; % Mean anomaly of Apollo orbit in degrees
rL = 5.725; % raduis of intial orbit around launch body (including raduis of body) (km)
rT = 3685; % raduis of capture orbit around target body (including raduis of body) (km)
muL = 1.33e-7;% gravatational parameter of launch body (apollo)

% Sun
muS = 1.33e+11; % gravatational parameter of sun

% Convert orbital elements position and velocity
% Cartesian --> elements
orb = [a3; e3; i3/180*pi; Omega3/180*pi; omega3/180*pi; M3/180*pi];
rv = E2C(orb, mu); % E2C is a function that converts elements to cartesian

% elements --> cartesian
orb_elements = C2E(rv, mu); % C2E is a function that converts cartesian to elements

% Hohmann transfer
% Calculating periapsis and apoapsis for Apollo
rp = a3 * (1 - e3); % Periapsis distance of Apollo orbit
ra = a3 * (1 + e3); % Apoapsis distance of Apollo orbit

% Semi-major axis of the Hohmann transfer orbit
a2 = (rp + ra) / 2;

% Velocity increments for Hohmann transfer 
delta_v1 = abs(sqrt((2*muS/a1)-(2*muS/(a1+a3))) - sqrt(muS/a1)); % Calculating deltav for first Hohmann transfer manouver
deltav_escape = abs(sqrt((2*muT/rT)+(delta_v1^2)) - sqrt((2*muT/rT)-(muT/a1))); % Calculating deltav for launch body escape
delta_v2 = abs(sqrt(muS/a3) - sqrt((2*muS/a3) - (2*muS/(a1+a3)))); % Calculating deltav for second Hohmann transfer manouver
deltav_capture = abs(sqrt(((2*muL)/rL) + (delta_v2^2)) - sqrt(((2*muL)/rL)-(muL/a3))); % Calculating deltav for target body capture 
Totaldv = abs(deltav_escape + delta_v2 + deltav_capture); % Total deltav of the transfer

% Display results
fprintf('Delta v1 (km/s): %.4f\n', delta_v1);
fprintf('Delta v2 (km/s): %.4f\n', delta_v2);
fprintf('Delta v Total (km/s): %.4f\n', Totaldv);

% Propellant mass
% delta_v = 2.43 ± 1.47 (assumed values)
Isp = 320; % Specific impulse of the rocket engine in seconds
g = 9.81; % Standard gravity in m/s^2
m2 = 3000; % kg, assumed spacecraft mass
m_Fuel = m2 * (exp((Totaldv)*10^3 / (g * Isp)) - 1);
fprintf('Propellant Mass Outgoing (kg): %.4f\n', m_Fuel);

% Velocity increments for Hohmann transfer Return Leg
delta_v1r = abs(sqrt((2*muS/a1)-(2*muS/(a1+a3))) - sqrt(muS/a1)); % Calculating deltav for first Hohmann transfer manouver
deltavr_escape = abs(sqrt((2*muL/rL)+(delta_v1r^2)) - sqrt((2*muL/rL)-(muL/a3))); % Calculating deltav for launch body escape
delta_v2r = abs(sqrt(muS/a1) - sqrt((2*muS/a1) - (2*muS/(a1+a3)))); % Calculating deltav for second Hohmann transfer manouver
deltavr_capture = abs(sqrt(((2*muT)/rT) + (delta_v2r^2)) - sqrt(((2*muT)/rT)-(muT/a1))); % Calculating deltav for target body capture 
Totaldvr = abs(deltavr_escape + delta_v2r + deltavr_capture); % Total deltav of the transfer

% Display results
fprintf('Delta v1r (km/s): %.4f\n', delta_v1r);
fprintf('Delta v2r (km/s): %.4f\n', delta_v2r);
fprintf('Delta vr Total (km/s): %.4f\n', Totaldvr);

% Propellant mass
% delta_v = 2.43 ± 1.47 (assumed values)
Isp = 320; % Specific impulse of the rocket engine in seconds
g = 9.81; % Standard gravity in m/s^2
m2 = 3000; % kg, assumed spacecraft mass
m_Fuel = m2 * (exp((Totaldvr)*10^3 / (g * Isp)) - 1);
fprintf('Propellant Mass Return (kg): %.4f\n', m_Fuel);

% Plot
% Create a theta array for plotting orbits
theta = linspace(0, 2*pi, 1000);         % For full orbits (Earth and Apollo)
theta_transfer = linspace(0, pi, 1000);  % For half of the transfer orbit (periapsis to apoapsis)

% Plot orbits with corrections
figure;
hold on;
grid on;
axis equal;

% Earth's orbit (or the initial orbit)
x_Earth = a1 * cos(theta); % X-coordinates for Earth's orbit
y_Earth = a1 * sin(theta); % Y-coordinates for Earth's orbit
plot(x_Earth, y_Earth, 'b', 'LineWidth', 1.5); % Plot Earth's orbit
legend_entries = {'Earth Orbit'};

% Apollo orbit (or the target orbit)
x_Apollo = a3 * (cos(theta) - e3); % X-coordinates for Apollo orbit (correcting for center offset)
y_Apollo = a3 * sqrt(1 - e3^2) * sin(theta); % Y-coordinates for Apollo orbit
plot(x_Apollo, y_Apollo, 'r', 'LineWidth', 1.5); % Plot Apollo orbit
legend_entries{end+1} = 'Apollo Orbit';

% Hohmann transfer orbit
x_transfer = a2 * cos(theta_transfer); % X-coordinates for transfer orbit
y_transfer = a2 * sin(theta_transfer); % Y-coordinates for transfer orbit
plot(x_transfer, y_transfer, 'g--', 'LineWidth', 1.5); % Plot transfer orbit
legend_entries{end+1} = 'Transfer Orbit';

% Plot the positions of the bodies at the starting point
plot(a1, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Earth starting point
plot(a3 * (1 - e3), 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Apophis periapsis point
legend_entries{end+1} = 'Earth Start';
legend_entries{end+1} = 'Apollo Periapsis';

% Set axis limits to improve visibility and keep the plot centered
xlim([-1.75 * 10^8, 1.75 * 10^8]); % Set X-axis limits based on the orbit size
ylim([-1.25 * 10^8, 1.25 * 10^8]); % Set Y-axis limits for better proportion

% Labels and title
xlabel('X position (km)');
ylabel('Y position (km)');
title('Hohmann Transfer between Earth and Apollo Orbits');
legend(legend_entries, 'Location', 'best');
hold off;
