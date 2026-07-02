function [point_position, point_amplitudes] = point_phantom()

    point_position(1,:) = [0 0 10e-3];
    point_position(2,:) = [0 0 20e-3];
    point_position(3,:) = [0 0 30e-3];
    point_position(4,:) = [0 0 40e-3];
    point_position(5,:) = [0 0 50e-3];
    point_position(6,:) = [0 0 60e-3];
    point_position(7,:) = [0 0 70e-3];
    point_position(8,:) = [0 0 80e-3];
    point_position(9,:) = [0 0 90e-3];
    point_position(10,:) = [0 0 100e-3];
    point_position(11,:) = [0 0 110e-3];
    
    point_position(12,:) = [-50e-3 0 40e-3];
    point_position(13,:) = [-40e-3 0 40e-3];
    point_position(14,:) = [-30e-3 0 40e-3];
    point_position(15,:) = [-20e-3 0 40e-3];
    point_position(16,:) = [-10e-3 0 40e-3];
    point_position(17,:) = [10e-3 0 40e-3];
    point_position(18,:) = [20e-3 0 40e-3];
    point_position(19,:) = [30e-3 0 40e-3];
    point_position(20,:) = [40e-3 0 40e-3];
    point_position(21,:) = [50e-3 0 40e-3];
    
    point_amplitudes = ones(size(point_position,1),1);

end

