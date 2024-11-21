% Constants
MuS = 1.32712440018e11; % Gravitational parameter of the Sun (Km^3/s^2)
MuA = 1.33e-7; % Gravitational parameter of Apollo (Km^3/s^2)
au = 1.496e11; % Astronomical Unit in meters

% Orbit parameters
orbit_start_date = datetime(2040, 4, 28);
orbit_end_date = datetime(2042, 7, 20);
time_step = 0.02083; % in days

% Communication parameters
eclipse_angle = 10; % degrees cone angle for eclipse analysis

% Initialise time array
time = orbit_start_date:days(time_step):orbit_end_date;
num_steps = length(time);

% Obtain ephemeris data for Apollo
r_apollo_sun = zeros(3, num_steps);
for i = 1:num_steps
    current_time = time(i);
    % Get Apollo positions from ephemeris (in heliocentric frame)
    ephemeris_data_apollo = Apollo_Ephemeris('Apollo');
    r_apollo_sun(:, i) = ephemeris_data_apollo(1:3); % Use first three elements as position vector
end

% Assuming initial conditions for spacecraft orbiting Apollo
% Set initial position and velocity around Apollo
r0_apollo = [0; 0; 5]; % Initial position 5 km away from Apollo
v0_apollo = [0; 0.21 / 1000; 0]; % Initial velocity for a stable orbit in km/s
y0 = double([r0_apollo; v0_apollo]); % Initial state vector (position and velocity)

% Integrate the spacecraft trajectory using initial conditions
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, state] = ode45(@(t, y) two_body_equation(t, y, MuA), double([0, days(orbit_end_date - orbit_start_date) * 86400]), y0, options); % Adjust time span for orbital period

% Extract spacecraft positions over time
r_spacecraft = state(:, 1:3)';

% Eclipse analysis around Apollo asteroid
R_apollo = 0.775; % Radius of Apollo in km
eclipse_status = ones(1, length(t)); % Initially assume no eclipse

for i = 1:length(t)
    % Find the closest time index in the ephemeris data for current spacecraft time
    ephem_index = min(i, num_steps);
    % Vector from Apollo to spacecraft
    r_apollo = r_apollo_sun(:, ephem_index); % Get Apollo position at the closest available ephemeris time step
    r_spacecraft_current = r_spacecraft(:, i);
    r_apollo_spacecraft = r_spacecraft_current - r_apollo;

    % Vector from Apollo to Sun (negative of Apollo's position in heliocentric frame)
    r_apollo_sun_vec = -r_apollo;

    % Calculate whether the spacecraft falls in Apollo's shadow
    % Step 1: Check if the spacecraft is behind Apollo relative to the Sun
    if dot(r_apollo_sun_vec, r_apollo_spacecraft) > 0
        % Step 2: Check if the spacecraft is within the shadow cone
        % Calculate the angular radius of Apollo as seen from the spacecraft
        D = norm(r_apollo_spacecraft); % Distance between Apollo and the spacecraft
        angular_radius_apollo = atan(R_apollo / D) * (180 / pi); % Angular radius in degrees

        % Calculate the angular separation between Apollo-Sun vector and Apollo-Spacecraft vector
        angular_separation = acosd(dot(r_apollo_spacecraft, r_apollo_sun_vec) / (norm(r_apollo_spacecraft) * norm(r_apollo_sun_vec)));

        % If the angular separation is less than the angular radius, the spacecraft is in shadow
        if angular_separation < angular_radius_apollo
            eclipse_status(i) = 0; % Eclipse occurs
        end
    end
end

% Plotting eclipse status around Apollo
time_days = t / 86400; % Convert seconds to days
figure;
plot(time_days, eclipse_status, 'LineWidth', 2);
ylabel('Eclipse Status');
yticklabels({'Eclipse', 'No Eclipse'});
yticks([0 1]);
xlabel('Time (days)');
title('Eclipse Window during Spacecraft Orbit around Apollo');
grid on;

% Function to define the two-body orbital motion equation
function dydt = two_body_equation(~, y, mu)
    if length(y) ~= 6
        error('State vector must have exactly 6 elements (position and velocity).');
    end
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);

    % Derivative of the state vector
    dydt = zeros(6,1);
    dydt(1:3) = v;
    dydt(4:6) = -mu * r / r_norm^3;
end
