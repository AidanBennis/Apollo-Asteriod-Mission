function [kep] = Toutatis_Ephemeris(time)
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
a = 3.80472e+8;
e = 0.6247568946161821;
i = 0.4480321679365896*conv_rads;
Om = 125.3589084385662*conv_rads;
om = 277.877068984381*conv_rads;
Mo = 339.675583583357*conv_rads;


% True Anamoly From NEO Ephemeris does all this to get true Anamoly
% at the time being simulated in Lambert
n  = sqrt(mu_sun/a^3);
t0 = Mo*pi/180/n;

M = n * (time - t0); % Mean Anomaly [rad]

E = M; % Initial guess for E
tolerance = 1e-8; % Convergence tolerance
while true
    E_next = M + e * sin(E); % Iterative solution
    if abs(E_next - E) < tolerance
        break;
    end
    E = E_next;
end

% Step 4: Calculate True Anomaly (nu)
wom = 2 * atan2(sqrt(1 + e) * sin(E_next / 2), sqrt(1 - e) * cos(E_next / 2)); % True anomaly [rad]

kep = [a,e,i,Om,om,wom];
