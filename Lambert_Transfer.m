%% Define Constants

% General Constants
conv_rads = pi/180; % Converts degrees to radians
mu_sun = 1.33e+11; % km^3/s^3
au_to_km = 1.496e+8; % Conversion factor from au to km
AU = 149597885.092797; % Astronomical Unit in km
G = 6.67430e-11; % Gravitational constant [m^3 kg^-1 s^-2]

% Itokawa Asteroid 
e = 0.2802554393054337; % Eccentricity of Itokawa's orbit
a = 1.324136563113898 * au_to_km; % Semi-major axis of Itokawa's orbit in km (converted from AU)
Om = 69.07689978577788; % Longitude of ascending node of Itokawa's orbit in degrees
om = 162.8201670956601; % Argument of periapsis of Itokawa's orbit in degrees
M = 142.5740657740646; % Mean anomaly of Itokawa's orbit in degrees
r = 1.165; % Radius of initial orbit around Itokawa (Km)
mu = 2.33e-09; % Gravitational parameter of Itokawa
close_approach_distance = 0.04060 * AU;

% Earth
aE = 1 * AU; % Semi-major axis of Earth's orbit in km
eE = 0; % Eccentricity of Earth's orbit
iE = 0; % Inclination
muE = 3.99e+5 ; % gravitational parameter of Earth
rE = 24367.5; % Orbit around Earth in Km (6371 [radius] + 35786 [GTO])
rE1 = 6571; % 200Km LEO Orbit for end of mission

%% Lambert Transfer
% Velocity increments for Lambert transfer 
delta_v1 = 4.67217;
deltav_escape = abs(sqrt((2*muE/rE)+(delta_v1^2)) - sqrt((2*muE/rE)-(muE/aE))); % Calculating deltav for Earth escape
delta_v2 = 4.67217;
deltav_capture = abs(sqrt(((2*mu)/r) + (delta_v2^2)) - sqrt(((2*mu)/r)-(mu/a))); % Calculating deltav for Itokawa capture 
Totaldv = abs(deltav_escape + deltav_capture); % Total deltav of the transfer
fprintf('Delta v esacpe (km/s): %.4f\n', deltav_escape);
fprintf('Delta v capture (km/s): %.4f\n', deltav_capture);
fprintf('Delta v Total (km/s): %.4f\n', Totaldv);

% Velocity increments for Lambert transfer Return Leg
delta_v1r = 4.67512;
deltavr_escape = abs((sqrt((2*mu/r)+(delta_v1r^2)) - sqrt((2*mu/r)-(mu/a))) - delta_v1r); % Calculating deltav for Itokawa escape
Totaldvr = abs(deltavr_escape + delta_v1r); % Total deltav of the transfer
fprintf('Delta vr escape (km/s): %.4f\n', deltavr_escape);
fprintf('Delta vr Total (km/s): %.4f\n', Totaldvr)

Totaldv_Complete_Journey = Totaldv + Totaldvr;
fprintf('Full Journey Delta v (km/s): %.4f\n', Totaldv_Complete_Journey)

%% Fuel Budget
Isp = 461; % Specific impulse of the rocket engine in seconds
g = 9.81; % Standard gravity in m/s^2
m2 = 2232; % Spacecraft mass + RL10 Engine - with cargo in Kg
ve = (Isp*g)/1000; % Exhaust Velocity in Km/s
mi1 = m2/(exp(-Totaldvr/ve)); % Total Mass of Spacecraft (spacecraft + fuel)
mfuel1 = mi1 - m2; % Mass of Fuel for Return Leg
fprintf('Propellant Mass Return (kg): %.4f\n\n', mfuel1);

m3 = 1732; % Spacecraft mass + RL10 Engine - no cargo in Kg
mi2 = mi1/(exp(-Totaldv/ve)); % Earth to Itokawa Total Mass
mfuel2 = mi2 - m3; % Earth to Itokawa Fuel Mass
fprintf('Propellant Mass Outgoing (kg): %.4f\n\n', mfuel2);

% Total Mass
mfuel_launch = mfuel1 + mfuel2;
mi_launch = mfuel_launch + m3;
fprintf('Total Launch Mass (kg): %.4f\n', mi_launch);
fprintf('Total Fuel Mass (kg): %.4f\n\n', mfuel_launch);