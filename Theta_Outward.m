%% Pork-Chop Outward

% Date Setup
jd_start = 2464769.334537;   % 16 March 2036
jd_end = 2464791.334537;     % 7 April 2036

% Convert Julian Date to Modified Julian Date (MJD)
td_start = jd_start - 2400000.5;  
td_end = jd_end - 2400000.5;

% Time of Flight setup
ToF0 = T_transfer / 60 / 60 / 24; % Minimum Time of Flight (days) - from Hohmann
dToF = 730;                      % Maximum additional Time of Flight in days (2 years max ToF)
nsteps_i = 22;                   % Steps for departure date
nsteps_j = 22;                   % Steps for Time of Flight

% Initialise Variables
ti = zeros(1, nsteps_i);          % Array for departure dates
ToF = zeros(1, nsteps_j);         % Array for Time of Flight values
THETA_values = zeros(nsteps_i, nsteps_j); % Array to store THETA values

% Iterate over possible departure dates and times of flight
for i = 1:nsteps_i
    % Calculate departure date in MJD
    ti(i) = td_start + (td_end - td_start) * (i - 1) / (nsteps_i - 1);
    for j = 1:nsteps_j    
        % Calculate Time of Flight in days
        ToF(j) = ToF0 + dToF * (j - 1) / (nsteps_j - 1);
        % Calculate arrival date in MJD
        tj = ti(i) + ToF(j);

        % Keplerian parameters of Earth and Itokawa at the current times
        [RE_kep] = Earth_Ephemeris(ti(i));
        [RA_kep] = Itokawa_Ephemeris(tj);

        % Convert Keplerian parameters to Cartesian coordinates
        [RE, vE] = kep2cart(RE_kep, mu_sun);
        [RA, vA] = kep2cart(RA_kep, mu_sun);

        % Solve Lambert's problem to find initial and final velocities
        try
            [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(RE, RA, ToF(j) * 86400, mu_sun, 0, 0, 0);
            % Store THETA value in degrees
            THETA_values(i, j) = THETA*(180/pi);
        catch
            % In case Lambert's solver fails, assign NaN to indicate no solution
            THETA_values(i, j) = NaN;
        end
    end
end

% Plot Heatmap of THETA
figure;
imagesc(ToF, ti, THETA_values);
colorbar;
xlabel('Time of Flight (days)');
ylabel('Departure Time (MJD)');
title('Theta plot in (Degrees) for Pork-Chop Plot outward');
set(gca, 'YDir', 'normal');
