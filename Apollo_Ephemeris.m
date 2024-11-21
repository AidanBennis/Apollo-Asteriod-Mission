function [kep] = Apollo_Ephemeris(time)
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

% Base Keplerian Elements of Apollo
a = 2.200004e+8; % Semi-major axis [km]
e = 0.5598715965109902; % Eccentricity
i = 6.352454347151728 * conv_rads; % Inclination [rad]
Om = 35.55648326281653 * conv_rads; % Longitude of ascending node [rad]
om = 286.0329799312636 * conv_rads; % Argument of periapsis [rad]
Mo = 252.868692344222 * conv_rads; % Mean anomaly at reference time [rad]
time_Mo = 60633; % Reference time in MJD (2024-Nov-19)

% Mean motion
n = sqrt(mu_sun / a^3); % [rad/s]

% Calculate time difference in seconds from the reference time
timediff = (time - time_Mo) * 86400; % Convert time difference to seconds

% Calculate mean anomaly at the given time
meanAnomaly = Mo + n * timediff; % Mean anomaly at time [rad]
meanAnomaly = mod(meanAnomaly, 2 * pi); % Keep mean anomaly within 0 to 2*pi

% Calculate eccentric anomaly using iterative method (Kepler's Equation)
eccAnomaly = meanAnomaly; % Initial guess for eccentric anomaly
tolerance = 1e-6; % Set tolerance for iterative solution
maxIter = 100; % Maximum number of iterations

for k = 1:maxIter
    f = eccAnomaly - e * sin(eccAnomaly) - meanAnomaly;
    f_prime = 1 - e * cos(eccAnomaly);
    delta = -f / f_prime;
    eccAnomaly = eccAnomaly + delta;
    if abs(delta) < tolerance
        break;
    end
end

% Calculate true anomaly from eccentric anomaly
trueAnomaly = 2 * atan2(sqrt(1 + e) * sin(eccAnomaly / 2), sqrt(1 - e) * cos(eccAnomaly / 2));

% Update the true anomaly (wom)
wom = trueAnomaly;

% Return Keplerian elements
kep = [a, e, i, Om, om, wom];
end
