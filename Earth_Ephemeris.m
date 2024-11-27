function [kep] = Earth_Ephemeris(time)
%Ephemeris of the Apollo Asteroid
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

%% MAIN %%

conv_rads = pi/180; % Converts degrees to radians
mu_sun = 1.33e+11; % Gravitational parameter of the Sun in km^3/s^2

% Base Keplerian Elements of Earth 
a = 1.496e+8; % Semi-major axis in km (1 AU)
e = 0.0167; % Eccentricity
i = 0; % Inclination in radians
Om = 0; % Longitude of ascending node in radians
om = 102.93768193 * conv_rads; % Argument of perihelion in radians
T = 365.256 * 86400; %orbital period days to seconds


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

% Return Keplerian elements
kep = [a, e, i, Om, om, wom];
end
