%% DMP_nss class
%

classdef DMP_nss < DMP
       
    methods  (Access = public)
        %% DMP constructor.
        %  @param[in] N_kernels: the number of kernels
        %  @param[in] a_z: Parameter 'a_z' relating to the spring-damper system.
        %  @param[in] b_z: Parameter 'b_z' relating to the spring-damper system.
        %  @param[in] can_clock_ptr: Pointer to a DMP canonical system object.
        function this = DMP_nss(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr)
            
            if (nargin < 4), can_clock_ptr = CanonicalClock(); end
            if (nargin < 5), shape_attr_gating_ptr=SigmoidGatingFunction(1.0, 0.5); end
            this@DMP(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr);

        end
        
    end
    
    methods  (Access = protected)
        
        function f_scale = forcingTermScaling(this, g)

            f_scale = 1.0;

        end
        
    end
    

end

