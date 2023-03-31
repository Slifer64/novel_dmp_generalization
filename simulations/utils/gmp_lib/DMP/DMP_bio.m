% Implements the DMP from:
% "Biologically-inspired dynamical systems for movement generation: Automatic real-time goal adaptation and obstacle avoidance"
% DOI: 10.1109/ROBOT.2009.5152423

classdef DMP_bio < DMP_
    
    methods (Access = public)
        
        function this = DMP_bio(n_dofs, N_kernels)
            
            this@DMP_(n_dofs, N_kernels, 'spatial_scaling','none', 'kernel_fun','gaussian', 'gating','exp');
            
        end

        function f_g = goal_attractor(this, y, y_dot, tau)
           
            f_g = (this.K*(this.g - y) - tau*this.D*y_dot) / tau^2;
            
        end
        
        function f_s = shape_attractor(this, s, s_dot, s_ddot, tau)
           
            f_s = (this.K*this.forcing_term(s) - this.K*(this.g - this.y0)*this.gating_fun(s)) / tau^2;

        end
        
    end
    
    methods (Access = protected)
        
        function Fd = calcFd(this, s, y, dy, ddy, g, y0, tau)
            
            f_gat = this.gating_fun(s);
            Fd = (ddy*tau^2 + this.D*tau*dy - this.K*(g - y) + this.K*(g - y0).*f_gat) ./ (this.K * f_gat);
            
        end 

    end
    
end