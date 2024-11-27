function [kep] = Itokawa_Ephemeris(time)
% Ephemeris of the Arwen Asteroid
% OUTPUT
%   kep     Keplerian parameters. It is a 6 entry vector:
%               [a e i Om om wom]
%           where:
%               a is the semimajor axis [km];
%               e is the eccentricity;
%               i is the inclination [rad];
%               Om is the anomaly of the ascending node [rad];
%               om is the anomaly of the pericentre [rad];
%               wom is the true anomaly (from the pericentre) [rad].
%               time_Mo is the time at which the true anomaly occurs [MJD]
%
%   mass    Mass of the NEO [kg]. It can be read from the database, or, if
%           not available, estimated by an approximate equation.
%   M       Mean anomaly at time [rad].

%% MAIN %%

conv_rads = pi/180; % Converts degrees to radians
mu_sun = 1.33e+11; % km^3/s^3
au_to_km = 1.496e+8; % Conversion factor from au to km

% Base Keplerian Elements of the Asteroid (updated values from image)
a = 1.324136563113898 * au_to_km; % Semi-major axis [km]
e = 0.2802554393054337; % Eccentricity
i = 1.621180377898326 * conv_rads; % Inclination [rad]
Om = 69.07689978577788 * conv_rads; % Longitude of ascending node [rad]
om = 162.8201670956601 * conv_rads; % Argument of periapsis [rad]
T = 556.5415504754865 * 86400; %Orbital Period days to Seconds

t = time;
M = 2 * pi * (t / T); % Mean anomaly (rad)
    
% Solve Kepler's equation: M = E - e*sin(E) iteratively
E = M; % Initial guess for eccentric anomaly
tol = 1e-8; % Convergence tolerance
while true
    dE = (M - (E - e * sin(E))) / (1 - e * cos(E));
    E = E + dE;
    if abs(dE) < tol
        break;
    end
end

% Step 4: Calculate True Anomaly (wom)
wom = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2)); % True anomaly [rad]

% Output Keplerian Elements
kep = [a, e, i, Om, om, wom];
end