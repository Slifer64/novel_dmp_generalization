classdef DMP_classic_wrapper < matlab.mixin.Copyable
    
    methods (Access = public)
        
        function this = DMP_classic_wrapper(gmp, traj_scale)
            
            if (nargin < 2), traj_scale=TrajScale_Prop(gmp.numOfDoFs()); end
            
            this.gmp = gmp.deepCopy();
            this.n_dof = this.gmp.numOfDoFs();
            this.gmp.setScaleMethod(traj_scale);

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
        
        function f_g = goal_attractor(this, y, y_dot, tau)
           
            f_g = this.K*(this.gmp.getGoal() - y) - this.D*y_dot;
            
        end
        
        function f_s = shape_attractor(this, s, s_dot, s_ddot, tau)
           
            f_s = this.K*(this.getRefPos(s) - this.gmp.getGoal()) + this.D*this.getRefVel(s, s_dot) + this.getRefAccel(s, s_dot, s_ddot);
            
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
    
    properties (Access = public)
        
        K % DMP stiffness
        D % DMP damping
        
        gmp
        
    end
    
    properties (Access = protected)
        
        
        n_dof
        
    end
    
end