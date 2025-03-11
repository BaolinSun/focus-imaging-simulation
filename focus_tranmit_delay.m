function emit_delay = focus_tranmit_delay(probe, focus_point, sound_speed)
   
    aperture_positions = [probe.x_ele; probe.y_ele; probe.z_ele];
    aperture_positions = aperture_positions';

    % 计算每个阵元到焦点的距离
    distances = sqrt(sum((aperture_positions - focus_point).^2, 2));

    % 找到最大距离
    max_distance = max(distances);
    min_distance = min(distances);
    mid_distance = norm(focus_point);
    
        % 计算延时矩阵，使所有信号同时到达焦点
        delay_matrix = (max_distance - distances) / sound_speed;
    %    delay_matrix = -(distances - min_distance) / sound_speed;
    %    delay_matrix = -(distances - mid_distance) / sound_speed;
    
    emit_delay = delay_matrix';
end