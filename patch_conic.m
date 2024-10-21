% Aidan Bennis 
% Patched Conic Transfer Function

%% Inital Conditions

% Constants
muL = 1.33e-7;% gravatational parameter of launch body (apollo)
muT = 3.99e+5 ;% gravatational parameter of target body (earth)
muS = 1.33e+11; % gravatational parameter of sun
aL = 2.2e+8; % semi major axis of launch body (apollo) km
aT = 1.496e+8; % semi major axis of target body (earth) km


% Intial conditions
r1 = 2.2e+8; % raduis of launch body from sun (km)
r2 = 1.496e+8; % raduis of target body from sun (km)
rL = 5.725; % raduis of intial orbit around launch body (including raduis of body) (km)
rT = 3685; % raduis of capture orbit around target body (including raduis of body) (km)

%% Main

%Calculating deltav for first Hohmann transfer manouver
deltav1 = sqrt((2*muS/r2)-(2*muS/(r1+r2))) - sqrt(muS/r2);

%Calculating deltav for launch body escape
deltav_esacpe = sqrt((2*muL/rL)+(deltav1^2)) - sqrt((2*muL/rL)-(muL/aL));

%Calculating deltav for second Hohmann transfer manouver
deltav2 = sqrt(muS/r1) - sqrt((2*muS/r1) - (2*muS/(r1+r2)));

%Calculating deltav for target body capture 
deltav_capture = sqrt(((2*muT)/rT) + (deltav2^2)) - sqrt(((2*muT)/rT)-(muT/aT));

%Total deltav of the transfer
Totaldv = deltav_esacpe + deltav2 + deltav_capture;
