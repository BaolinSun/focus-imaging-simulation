classdef Probe
    properties
        type
        N_elements
        pitch
        x_ele
        y_ele
        z_ele
        ele_pos
        tx_ori
    end

    methods
        function obj = Probe(probe_type)
            obj.type = probe_type;
        end
    end
end