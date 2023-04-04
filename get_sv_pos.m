%% Calculate satellite position and clock bias based on emphemeris and time
function [x] = get_sv_pos(eph, pr, c, wedot, mu, F)

    % Input parameters:
    rcvr_tow    = eph(1);	% Receiver time of week (s)
    svid        = eph(2);	% Satellite PRN number (1-32)
    toc         = eph(3);	% Reference time of clock parameters (s)
    toe         = eph(4);	% Reference time of ephemeris parameters (s)
    af0         = eph(5);   % Clock correction coefficient - group delay (s)
    af1         = eph(6);	% Clock correction coefficient (s/s)
    af2         = eph(7);	% Clock correction coefficient (s/s/s)
    ura         = eph(8);	% User range accuracy (m)
    e           = eph(9);	% Eccentricity (-)
    sqrta       = eph(10);	% Square root of semi-major axis a (m^1/2)
    dn          = eph(11);	% Mean motion correction (r/s)
    m0          = eph(12);	% Mean anomaly at reference time (r)
    w           = eph(13);	% Argument of perigee (r)
    omg0        = eph(14);	% Right ascension (r)
    i0          = eph(15);  % Inclination angle at reference time (r)
    odot        = eph(16);  % Rate of right ascension (r/s)
    idot        = eph(17);  % Rate of inclination angle (r/s)
    cus         = eph(18);  % Argument of latitude correction, sine (r)
    cuc         = eph(19);  % Argument of latitude correction, cosine (r)
    cis         = eph(20);  % Inclination correction, sine (r)
    cic         = eph(21);  % Inclination correction, cosine (r)
    crs         = eph(22);  % Radius correction, sine (r)
    crc         = eph(23);  % Radius correction, cosine (r)
    iod         = eph(24);  % Issue of data number   
    
    % Calculation follows Table 20-IV in ICD file
    a = sqrta^2;            % Semi-major axis
    n0 = sqrt(mu/a^3);      % Computed mean motion (r/s)
    t = rcvr_tow - pr/c;    % Satellite signal transmition time
    tk = t - toe;           % Time from ephemeris reference epoch
    
    % Account for beginning or end of week crossovers
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -30240
        tk = tk + 604800;
    end
    
    n = n0+dn;              % Corrected mean motion
    mk = m0+n*tk;           % Mean anomaly
    
    % Solve eccentric anomaly
    syms ex
    eqn = ex - e*sin(ex) == mk;
    ek = vpasolve(eqn);
    ek = double(ek);
    clear ex eqn
    
    % True anomaly
    vk = atan2((sqrt(1-e^2)*sin(ek)/(1-e*cos(ek))),((cos(ek)-e)/(1-e*cos(ek))));
    
    % Eccentric anomaly
    phik = vk + w;          % Argument of latitude
    
    % Second harmonic perturbations
    duk = cus*sin(2*phik) + cuc*cos(2*phik);    % Argument of Latitude Correction
    drk = crs*sin(2*phik) + crc*cos(2*phik);    % Radius Correction
    dik = cis*sin(2*phik) + cic*cos(2*phik);    % Inclination Correction
    
    uk = phik + duk;                            % Corrected argument of latitude
    rk = a*(1 - e*cos(ek)) + drk;               % Corrected Radius
    ik = i0 + dik + idot*tk;                    % Corrected Inclination
    xkp = rk*cos(uk);                           % X position in orbital plane
    ykp = rk*sin(uk);                           % Y position in orbital plane
    omgk = omg0 + (odot - wedot)*tk -wedot*toe; % Corrected longitude of ascending node
    
    xk = xkp*cos(omgk) - ykp*cos(ik)*sin(omgk); % Earth-fixed X coordinates
    yk = xkp*sin(omgk) + ykp*cos(ik)*cos(omgk); % Earth-fixed Y coordinates
    zk = ykp*sin(ik);                           % Earth-fixed Z coordinates
    
    % Calculate satellite clock bias (s)
    dtsv = af0 + af1*(t - toc) + af2*(t - toc)^2 + F*e*sqrta*sin(ek);
    
    x = [xk yk zk dtsv];
end