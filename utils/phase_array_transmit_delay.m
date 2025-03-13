function emit_delay = phase_array_transmit_delay(probe, angle_rad, F, c)
    
    xx = F * sin(angle_rad);
    zz = F * cos(angle_rad);
    yy = 0 * xx;

    focus_point = [xx; yy; zz];
    focus_point = focus_point';

    aperture_positions = [probe.x_ele; probe.y_ele; probe.z_ele];
    aperture_positions = aperture_positions';

    % Calculate distances from each element to the focus point
    L = sqrt(sum((aperture_positions - focus_point).^2, 2));

    emit_delay = 5e-6 - (L - F) / c;
    emit_delay = emit_delay';
end



% element_num = 64;
% width = 0.12e-3;         % 阵元宽度 (m)
% kerf = 0.18e-3;          % 阵元间距 (m)
% pitch = width + kerf;

% % 阵元参数
% x_ele = ([0:element_num-1]-(element_num-1)/2).*pitch;
% z_ele  = zeros(1,length(x_ele));
% probe.element_num = element_num;
% probe.pitch = pitch;
% probe.x_ele = x_ele;
% probe.y_ele = z_ele;
% probe.z_ele = z_ele;
% probe.ele_pos = [x_ele; z_ele; z_ele]';

% F = 120e-3;
% scan_angle_deg = -30:2:30;
% num_lines = length(scan_angle_deg);

% for line = 1:num_lines
%     angle_deg = scan_angle_deg(line);
%     angle_rad = angle_deg * pi / 180;

%     emit_delay = transmit_delay(probe, angle_rad, F, 1540);

%     plot(emit_delay);
%     hold on;
% end