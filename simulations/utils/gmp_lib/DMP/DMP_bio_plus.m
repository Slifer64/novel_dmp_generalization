% Implements the DMP from:
% "Overcoming some drawbacks of Dynamic Movement Primitives"
% DOI: https://doi.org/10.1016/j.robot.2021.103844

classdef DMP_bio_plus < DMP_bio
    
    methods (Access = public)
        
        function this = DMP_bio_plus(n_dofs, N_kernels)
            
            this@DMP_bio(n_dofs, N_kernels);
            
            this.setKernelFun('mollifier');
            this.setSpatialScaling('roto-dial');
            this.setGatingFun('linear');
        end
        
    end
    
end