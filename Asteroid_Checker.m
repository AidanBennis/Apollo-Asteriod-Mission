%% Orbital Parameters

% Constants
conv_rads = pi/180; % Converts degrees to radians
mu_sun = 1.33e+11; % km^3/s^3
au_to_km = 1.496e+8; % Conversion factor from au to km
AU = 149597885.092797; % Astronomical Unit in km
G = 6.67430e-11; % Gravitational constant [m^3 kg^-1 s^-2]

% Asteroid: 
e =  0.2802554393054337; % Eccentricity of Arwen' orbit
a = 1.324136563113898 * au_to_km; % Semi-major axis of Arwen' orbit in km (converted from AU)
Om = 69.07689978577788; % Longitude of ascending node of Arwen' orbit in degrees
om = 162.8201670956601; % Argument of periapsis of Arwen' orbit in degrees
M = 142.5740657740646; % Mean anomaly of Arwen' orbit in degrees
r = 2.2; % Radius of initial orbit around Arwen (Km)
mu = 1.88E-10; % Gravitational parameter of Arwen (assumed similar to Arwen for now)
close_approach_distance = 0.04060		 * AU;

% Earth
aE = 1 * AU; % Semi-major axis of Earth's orbit in km
eE = 0; % Eccentricity of Earth's orbit
iE = 0; % Inclination
muE = 3.99e+5 ; % gravatational parameter of Earth
rE = 42157; % Orbit around Earth in Km (6371 [radius] + 35786 [GTO])
rE1 = 6571; % 200Km LEO Orbit for end of mission

%% Hohmann Transfer
% Velocity increments for Hohmann transfer 
delta_v1 = abs(sqrt((2*mu_sun/aE)-(2*mu_sun/(aE+a))) - sqrt(mu_sun/aE)); % Calculating deltav for first Hohmann transfer manouver
deltav_escape = abs(sqrt((2*muE/rE)+(delta_v1^2)) - sqrt((muE/rE))); % Calculating deltav for Earth escape
delta_v2 = abs(sqrt(mu_sun/a) - sqrt((2*mu_sun/a) - (2*mu_sun/(aE+a)))); % Calculating deltav for second Hohmann transfer manouver
deltav_capture = abs(sqrt(((2*mu)/r) + (delta_v2^2)) - sqrt(((2*mu)/r)-(mu/a))); % Calculating deltav for Arwen capture 
Totaldv = abs(deltav_escape + delta_v2 + deltav_capture); % Total deltav of the transfer
fprintf('Delta v Total (km/s): %.4f\n', Totaldv);

% Velocity increments for Hohmann transfer Return Leg
delta_v1r = abs(sqrt((2*mu_sun/a)-(2*mu_sun/(aE+a))) - sqrt(mu_sun/a)); % Calculating deltav for first Hohmann transfer manouver
deltavr_escape = abs(sqrt((2*mu/r)+(delta_v1r^2)) - sqrt((2*mu/r)-(mu/a))); % Calculating deltav for Arwen escape
delta_v2r = abs(sqrt(mu_sun/aE) - sqrt((2*mu_sun/aE) - (2*mu_sun/(aE+a)))); % Calculating deltav for second Hohmann transfer manouver
deltavr_capture = abs(sqrt(((2*muE)/rE1) + (delta_v2r^2)) - sqrt(((2*muE)/rE1)-(muE/aE))); % Calculating deltav for Earth capture 
Totaldvr = abs(deltavr_escape + delta_v2r + deltavr_capture); % Total deltav of the transfer
fprintf('Delta v1r (km/s): %.4f\n', delta_v1r);
fprintf('Delta vr Total (km/s): %.4f\n', Totaldvr)

Totaldv_Complete_Journey = Totaldv + Totaldvr;
fprintf('Full Journey Delta v (km/s): %.4f\n', Totaldv_Complete_Journey)

% Hohmann Transfer Orbit for Close Approach
rp = close_approach_distance; % Periapsis distance of Arwen at close approach
a2 = (aE + rp) / 2; % Semi-major axis for the Hohmann transfer orbit
e2 = (rp - aE) / (rp + aE); % Eccentricity of the transfer orbit
T_transfer = pi * sqrt(a2^3 / mu_sun); % Time for half of the Hohmann transfer
fprintf('Hohmann Transfer Time: %.4f seconds', T_transfer);
fprintf(' (%.2f days)\n\n', T_transfer / 60 / 60 / 24);

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

%Total Mass
mfuel_launch = mfuel1 + mfuel2;
mi_launch = mfuel_launch + m3;
fprintf('Total Launch Mass (kg): %.4f\n', mi_launch);
fprintf('Total Fuel Mass (kg): %.4f\n\n', mfuel_launch);


%% Pork Chop Outward
% Constants
muS = 1.32712440018e11; % Sun's gravitational parameter in km^3/s^2
mu = muS;

% Date Setup
% Start and end of possible departure dates in Julian Date (JD)
jd_start = 2464328.9443287;  % January 1, 2035
jd_end = 2466519.9443171;  % December 31, 2040

% Convert Julian Date to Modified Julian Date (MJD)
td_start = jd_start - 2400000.5;  % Convert to Modified Julian Date (MJD)
td_end = jd_end - 2400000.5;  % Convert to Modified Julian Date (MJD)

% Time of Flight setup
ToF0 = T_transfer / 60 / 60 / 24; % Minimum Time of Flight (days) - from Hohmann
dToF = 730; % Maximum additional Time of Flight in days (2 years max ToF)
nsteps_i = 200; % Steps for departure date
nsteps_j = 200; % Steps for Time of Flight

% Initialise Variables
ti = zeros(1, nsteps_i); % Array for departure dates
ToF = zeros(1, nsteps_j); % Array for Time of Flight values
dv = zeros(nsteps_i, nsteps_j); % Delta-v array

% Iterate over possible departure dates and times of flight
for i = 1:nsteps_i
    % Calculate departure date in MJD
    ti(i) = td_start + (td_end - td_start) * (i - 1) / (nsteps_i - 1);
    for j = 1:nsteps_j    
        % Calculate Time of Flight in days
        ToF(j) = ToF0 + dToF * (j - 1) / (nsteps_j - 1);
        % Calculate arrival date in MJD
        tj = ti(i) + ToF(j);

        % Keplerian parameters of Earth and Apollo at the current times
        [RE_kep] = Earth_Ephemeris(ti(i));
        [RA_kep] = Itokawa_Ephemeris(tj);

        % Convert Keplerian parameters to Cartesian coordinates
        [RE, vE] = kep2cart(RE_kep, mu);
        [RA, vA] = kep2cart(RA_kep, mu);

        % Solve Lambert's problem to find initial and final velocities
        try
            [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(RE, RA, ToF(j) * 86400, muS, 0, 0, 0);
            % Calculate the total delta-v needed for the mission
            dv(i, j) = (norm(VI - vE)) + (norm(VF - vA));
        catch
            dv(i, j) = NaN; % Assign NaN if the Lambert solver fails
        end
    end
end

% Filter Delta-v values
% Set an upper threshold for delta-v to highlight smaller values
max_dv_threshold = 20; % For example, you want to emphasize values less than 20 km/s
dv(dv > max_dv_threshold) = NaN; % Assign NaN to higher values to ignore them in the plot

% Generate Pork-chop Plot with adjusted limits
figure;
contourf(ti, ToF, dv', 40, 'ShowText', 'on'); % Increase the number of contour levels to improve granularity
xlabel('Departure Date (MJD)');
ylabel('Time of Flight (days)');
title('Earth to Itokawa Transfer Pork-Chop Plot');
colorbar;
clim([0 max_dv_threshold]); % Set the color scale to emphasize smaller delta-v values
grid on;

%% Pork-Chop Plot for Return Transfer

% Constants
muS = 1.32712440018e11; % Sun's gravitational parameter in km^3/s^2
mu = muS;

% Date Setup
% Start and end of possible departure dates in Julian Date (JD)
jd_start = 2466521.216701;  % January 1, 2041 (departure from asteroid)
jd_end = 2468346.216701;  % December 31, 2045

% Convert Julian Date to Modified Julian Date (MJD)
td_start = jd_start - 2400000.5;  % Convert to Modified Julian Date (MJD)
td_end = jd_end - 2400000.5;  % Convert to Modified Julian Date (MJD)

% Time of Flight setup
ToF0 = T_transfer / 60 / 60 / 24; % Minimum Time of Flight (days)
dToF = 730; % Maximum additional Time of Flight in days (2 years max ToF)
nsteps_i = 100; % Steps for departure date
nsteps_j = 100; % Steps for Time of Flight

% Initialise Variables
ti = zeros(1, nsteps_i); % Array for departure dates
ToF = zeros(1, nsteps_j); % Array for Time of Flight values
dv = zeros(nsteps_i, nsteps_j); % Delta-v array

% Iterate over possible departure dates and times of flight
for i = 1:nsteps_i
    % Calculate departure date in MJD from the asteroid
    ti(i) = td_start + (td_end - td_start) * (i - 1) / (nsteps_i - 1);
    for j = 1:nsteps_j    
        % Calculate Time of Flight in days
        ToF(j) = ToF0 + dToF * (j - 1) / (nsteps_j - 1);
        % Calculate arrival date in MJD
        tj = ti(i) + ToF(j);

        % Keplerian parameters of Arwen (asteroid) at the departure time
        [RA_kep] = Asteroid_Ephemeris(ti(i)); % Asteroid at departure time
        [RE_kep] = Earth_Ephemeris(tj); % Earth at arrival time

        % Convert Keplerian parameters to Cartesian coordinates
        [RA, vA] = kep2cart(RA_kep, mu);
        [RE, vE] = kep2cart(RE_kep, mu);

        % Solve Lambert's problem to find initial and final velocities
        try
            [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(RA, RE, ToF(j) * 86400, muS, 0, 0, 0);
            % Calculate the total delta-v needed for the mission
            dv(i, j) = (norm(VI - vA)) + (norm(VF - vE));
        catch
            dv(i, j) = NaN; % Assign NaN if the Lambert solver fails
        end
    end
end

% Filter Delta-v values
% Set an upper threshold for delta-v to highlight smaller values
max_dv_threshold = 20; % For example, you want to emphasize values less than 20 km/s
dv(dv > max_dv_threshold) = NaN; % Assign NaN to higher values to ignore them in the plot

% Generate Pork-chop Plot with adjusted limits
figure;
contourf(ti, ToF, dv', 40, 'ShowText', 'on'); % Increase the number of contour levels to improve granularity
xlabel('Departure Date (MJD)');
ylabel('Time of Flight (days)');
title('Itokawa to Earth Transfer Pork-Chop Plot');
colorbar;
clim([0 max_dv_threshold]); % Set the color scale to emphasize smaller delta-v values
grid on;
