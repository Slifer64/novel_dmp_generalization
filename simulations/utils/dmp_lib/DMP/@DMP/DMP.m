%% DMP class
%

classdef DMP < DMP_
       
    methods (Access = public)
        %% DMP constructor.
        %  @param[in] N_kernels: the number of kernels
        %  @param[in] a_z: Parameter 'a_z' relating to the spring-damper system.
        %  @param[in] b_z: Parameter 'b_z' relating to the spring-damper system.
        %  @param[in] can_clock_ptr: Pointer to a DMP canonical system object.
        function this = DMP(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr)
            
            if (nargin < 4), can_clock_ptr = CanonicalClock(); end
            if (nargin < 5), shape_attr_gating_ptr=SigmoidGatingFunction(1.0, 0.5); end
            this@DMP_(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr);

        end
        
        
        function shape_attr = shapeAttractor(this, x, g)

            s_attr_gating = this.shapeAttrGating(x);
            f_scale = this.forcingTermScaling(g);
            shape_attr = s_attr_gating * f_scale * this.forcingTerm(x);

        end
        
    end
    
    methods (Access = protected)
        
        function Fd = calcFd(this, x, y, dy, ddy, g)

            tau = this.getTau();
            Fd = (ddy*tau^2 - this.goalAttractor(x, y, tau*dy, g));

        end
        
        
        function Fd = calcLearnedFd(this, x, g)
            
            Fd = this.shapeAttractor(x, g);

        end
        
        
        function f_scale = forcingTermScaling(this, g)

            f_scale = (g - this.y0);

        end

    end
    

end

