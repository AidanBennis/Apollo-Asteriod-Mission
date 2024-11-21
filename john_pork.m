%% Pork-Chop plot for Earth to Asteroid Transfer
clear all;

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
ToF0 = 101; % Minimum Time of Flight (days)
dToF = 730; % Maximum additional Time of Flight in days (2 years max ToF)
nsteps_i = 2191; % Steps for departure date
nsteps_j = 730; % Steps for Time of Flight

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
        [RA_kep] = Apollo_Ephemeris(tj);

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

% Generate Pork-chop Plot
figure;
contourf(ti, ToF, dv', 'ShowText', 'on');
xlabel('Departure Date (MJD)');
ylabel('Time of Flight (days)');
title('Earth to Asteroid Transfer Pork-Chop Plot');
colorbar;
grid on;


%% Date Finder
% Display corresponding date and time of flight for a specific cell
i = 1665; % 
j = 179; % 
fprintf('For dv(%d, %d): Departure Date = %.2f MJD, Time of Flight = %.2f days\n', i, j, ti(i), ToF(j));
