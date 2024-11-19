function [kep] = Earth_Ephemeris(time)
%Ephemeris of the Apollo Asteriod
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
%               time_Mo is the time at which the true anamoly ocurs [MJD]
%
%   mass    Mass of the NEO [kg]. It can be read from the database, or, if
%           not available, estimated by an approximate equation.
%   M       Mean anomaly at time [rad].

%% MAIN %%

conv_rads = pi/180; % Converets degrees to rads
mu_sun = 1.33e+11; %km3/s3

    % Base Keplerian Elements of Earth 
    a = 1.496e+8; % Semi-major axis in km (1 AU)
    e = 0.0167; % Eccentricity
    i = 0; % Inclination in radians
    Om = 0; % Longitude of ascending node in radians
    om = 102.93768193 * conv_rads; % Argument of perihelion in radians
    Mo = 358.617 * conv_rads; % Mean anomaly at reference time in radians
    time_Mo = 60633; % Reference time (MJD for 2024-Nov-19)


% True Anamoly From NEO Ephemeris does all this to get true Anamoly
% at the time being simulated in Lambert
n  = sqrt(mu_sun/a^3);
t0 = Mo*pi/180/n;
timediff = time_Mo - 51544.5; % Convert to MJD 2000
t = (time - timediff)*86400+t0;

p = 2*pi*sqrt(a^3/mu_sun);
np = floor(t/p);
t = t-p*np;
phi = n*t;
wom=2*atan(sqrt((1+e)/(1-e))*tan(phi*0.5));

kep = [a,e,i,Om,om,wom];
