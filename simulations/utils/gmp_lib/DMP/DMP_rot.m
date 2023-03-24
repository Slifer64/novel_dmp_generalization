% Implements the DMP from:
% "A novel DMP formulation for global and frame independent spatial scaling in the task space"
% DOI: 10.1109/RO-MAN47096.2020.9223500

classdef DMP_rot < DMP_classic
    
    methods (Access = public)
        
        function this = DMP_rot(n_dofs, N_kernels)
            
            this@DMP_classic(n_dofs, N_kernels);
            this.setSpatialScaling('rot-min');
            
        end
        
    end
    
end