classdef DMP_classic < matlab.mixin.Copyable
    
    methods (Access = public)
        
        function this = DMP_classic(gmp)
            
            this.gmp = gmp.deepCopy();
            this.n_dof = this.gmp.numOfDoFs();
            this.gmp.setScaleMethod(TrajScale_Prop(this.n_dof));
            
        end
        
        function setAdaptToRobot(this, set)

        end
        
        function setRecursiveUpdate(this, set)

            
        end
        
        function reset(this)
            
        end
        
        function setOptMetric(this, metric)
           

        end
        
        function init(this, s, y0, yg, Tf)

            this.gmp.setY0(y0);
            this.gmp.setGoal(yg);
            
        end
            
        function update(this, yg, s, s_dot, y, y_dot)
            
            this.gmp.setGoal(yg);

        end
        
        function s = updateViapoint(this, s_v, pos, search_closest, reverse)
        

        end
        
        function p_ref = getRefPos(this, s)
            
            p_ref = this.gmp.getYd(s);
            
        end
        
        function v_ref = getRefVel(this, s, s_dot)
            
            v_ref = this.gmp.getYdDot(s, s_dot);
            
        end
        
        function a_ref = getRefAccel(this, s, s_dot, s_ddot)
            
            a_ref = this.gmp.getYdDDot(s, s_dot, s_ddot);
            
        end

        %% for compatibility with GMP_MPC:
        function out = getYd(this, s)
            out = this.getRefPos(s);
        end
        
        function out = getYdDot(this, s, s_dot)
            out = this.getRefVel(s, s_dot);
        end
        
        function n_dof = numOfDoFs(this)
            n_dof = this.n_dof;
        end
        
    end
    
    properties (Access = protected)
        
        gmp
        n_dof
        
    end
    
end