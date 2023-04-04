% AAE6102 LAB SESSION LO, Li-yu 22039044R
clear all
clc;

% INITIAL POSITION & CLOCK BIAS
approx = [-2694685.473;  
          -4293642.366;
          3857878.924;
          0];
      
c = 299792458.0;            % SPEED OF LIGHT (m/s)
wedot = 7.2921151467e-5;    % WGS 84 VALUE OF EARTH’S ROTATION RATE (r/s)
mu = 3.986005e+14;          % WGS 84 VALUE OF EARTH'S UNIVERSAL GRAVITATION CONSTANT (m^3/s^2)
F = (-4.442807633e-10);     % RELATIVISTIC CORRECTION TERM CONSTANT
dx_threshold = 1e-4;        % THRESHOLD TO STOP THE ITERATION
         
eph = readmatrix('eph.dat');
rcvr = readmatrix('rcvr.dat');
eph_ = sortrows(eph,2);
rcvr_ = sortrows(rcvr,2);

% CALCULATE SATELLITE POSITION AND CLOCK BIAS
sv_num = size(eph_,1);       
sv_data = [];
for i = 1:sv_num
    sv_data(i,:) = get_sv_pos(eph_(i,:),rcvr_(i,3), c, wedot, mu, F);
end

sv_data(:,4) = rcvr_(:,3) + c .* sv_data(:,4);


% USE LEAST SQUARE TO CALCULATE USER POSITION AND CLOCK BIAS
dx = [-5710, 1080, -2610, 519450]; 

i = 0;                      

% BELOW START THE MAIN TASK
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('CALCULATING USER POSITION AND CLOCK BIAS BASED ON EPHEMERIS DATA AND PSEUDORANGE')
fprintf('\n')
disp('INITIAL STATE IN WGS 84 XYZ COORDINATE (m) and USER CLOCK BIAS (s): ')
disp(num2str(approx(1:3)))
disp(num2str(approx(4)/c))
fprintf('\n')
fprintf('\n')

while norm(dx(1:3)) > dx_threshold
    dx = estimate_user_pos(sv_data,approx);
    approx = approx + dx;   % Update approximation
    i = i+1;
    D = ['Iteration ',num2str(i),' result:'];
    disp(D)
    disp(num2str(approx(1:3)))
    disp(num2str(approx(4)/c))
    fprintf('\n')
    fprintf('\n')
end

% Print the results
disp('ERROR < PRESET THRESHOLD, ITERATION ENDS')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')








