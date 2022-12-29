classdef Obstacle < matlab.mixin.Copyable
    
    methods (Access = public)
        
        function this = Obstacle(d0, color)
 
            if (nargin < 2), color=[]; end
            
            this.color = color;
            this.d0 = d0;
            
        end
        
        function f_rep = repulsive_force(this, y, gain)
    
            psi = this.surface_fun(y);
            psi_grad = this.surface_fun_grad(y);
            
            if (psi < 0)
                warning('The obstacle has been penetrated!');
            end

            if (psi > this.d0)
                f_rep = zeros(size(y));
            else
                e = (psi-this.d0)^2 / this.d0^2;
                f_rep = -gain * (2*(psi-this.d0) / (this.d0^2*(1-e))) * psi_grad;
            end

        end
        
    end
    
    methods (Abstract, Access = public)

        plot(this, varargin)
        
    end
    
    methods (Abstract, Access = protected)
        
        surface_fun(this, y)
        
        surface_fun_grad(this, y)
        
    end
    
    properties (Access = protected)
       
        d0 % tolerance to activate the repulsive force
        color
        
    end
    
end