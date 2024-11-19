function [kep] = Apollo_Ephemeris(time)
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

% Base Keplerian Elements of Apollo
a = 2.200004e+8;
e = 0.5598715965109902;
i = 6.352454347151728*conv_rads;
Om = 35.55648326281653*conv_rads;
om = 286.0329799312636*conv_rads;
Mo = 252.868692344222*conv_rads;
time_Mo = 60600; % 2024-Oct-17 in MJD

% True Anamoly From NEO Ephemeris does all this to get true Anamoly
% at the time being simulated in Lambert
n  = sqrt(mu_sun/a^3);
t0 = Mo*pi/180/n;
timediff = time_Mo - 51544.5; % Convert to MJD2000
t = (time - timediff)*86400+t0;

p = 2*pi*sqrt(a^3/mu_sun);
np = floor(t/p);
t = t-p*np;
phi = n*t;
ddf=(e*cos(phi)-1);
phi=phi-t*n/ddf+(phi-e*sin(phi))/ddf;
wom=2*atan(sqrt((1+e)/(1-e))*tan(phi*0.5));

kep = [a,e,i,Om,om,wom];
