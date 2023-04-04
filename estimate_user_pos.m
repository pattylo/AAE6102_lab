%% Estimate user position and clock bias based on pseudorange and satellite position
function [dx] = estimate_user_pos(sv_data,approx)
    
    % Input parameters
    sv_x = sv_data(:,1);    % Satellite X position
    sv_y = sv_data(:,2);    % Satellite Y position
    sv_z = sv_data(:,3);    % Satellite Z position
    pr = sv_data(:,4);      % Measured pseudorange
    x0 = approx(1);         % Initial X user position
    y0 = approx(2);         % Initial Y user position
    z0 = approx(3);         % Initial Z user position
    b0 = approx(4);         % Initial user clock bias (m)
    
    % Calculate delta rho
    rho = pr;               % Current (measured) pseudorange
    r = sqrt((sv_x - x0).^2 + (sv_y - y0).^2 + (sv_z - z0).^2);
                            % Geometric range between user and satellite
    rho_hat = r - b0;       % Last (approximated) pseudorange
    delta_rho = rho_hat - rho;
    
    % Calculate H matrix
    ax = (sv_x - x0)./r;
    ay = (sv_y - y0)./r;
    az = (sv_z - z0)./r;
    H = [ax ay az ones(size(sv_data,1),1)];
    
    % Calculate dx
    dx = inv(H'*H)*H'*delta_rho;
    
    lala =0;

end