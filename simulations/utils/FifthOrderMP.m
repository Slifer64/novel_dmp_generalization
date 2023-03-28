classdef FifthOrderMP < matlab.mixin.Copyable
    
    methods (Access = public)
        
        function this = FifthOrderMP(n_dof)
            
            this.n_dof = n_dof;
            this.a_fifth = zeros(this.n_dof, 6);
            
        end
        
        function update(this, s1, h1, h1_dot, h1_ddot, s2, h2, h2_dot, h2_ddot, s_dot, s_ddot)
            
            A1 = this.fifth_pol_mat(s1, s_dot, s_ddot, 1, ~any(isnan(h1_dot(:))), ~any(isnan(h1_ddot(:))));
            A2 = this.fifth_pol_mat(s2, s_dot, s_ddot, 1, ~any(isnan(h2_dot(:))), ~any(isnan(h2_ddot(:))));
            A = [A1 A2];
            b = h1;
            if ~any(isnan(h1_dot(:))), b = [b h1_dot]; end
            if ~any(isnan(h1_ddot(:))), b = [b h1_ddot]; end
            b = [b h2];
            if ~any(isnan(h2_dot(:))), b = [b h2_dot]; end
            if ~any(isnan(h2_ddot(:))), b = [b h2_ddot]; end
            
%             err_0 = vecnorm(this.a_fifth * A - b, 2, 1);
            
            err_0 = vecnorm(this.a_fifth * A - b, 2, 1);

            this.a_fifth = b / A;
            
            err_1 = vecnorm(this.a_fifth * A - b, 2, 1);
            
%             disp('--------------------')
%             s1
%             A
%             err_0
%             err_1
%             pause
        end
        
        function A = fifth_pol_mat(this, s, s_dot, s_ddot, pos, vel, accel)
            
            A = [];
            if (pos), A = [A this.fifth_phi(s)]; end
            if (vel), A = [A this.fifth_phi_dot(s, s_dot)]; end
            if (accel), A = [A this.fifth_phi_ddot(s, s_dot, s_ddot)]; end
            
        end
        
            
        function phi = fifth_phi(this, s)
            
            phi = [1 s s^2 s^3 s^4 s^5]';
            
        end
        
        function phi_dot = fifth_phi_dot(this, s, s_dot)
            
            phi_dot = [0 1 2*s 3*s^2 4*s^3 5*s^4]'*s_dot;
            
        end
        
        function phi_ddot = fifth_phi_ddot(this, s, s_dot, s_ddot)
            
            phi_ddot = [0 0 2 6*s 12*s^2 20*s^3]'*s_dot^2 + [0 1 2*s 3*s^2 4*s^3 5*s^4]'*s_ddot;
            
        end
        
        function pos = getRefPos(this, s)
            
            pos = this.a_fifth * this.fifth_phi(s);
            
        end
        
        function vel = getRefVel(this, s, s_dot)
            
            vel = this.a_fifth * this.fifth_phi_dot(s, s_dot);
            
        end
        
        function accel = getRefAccel(this, s, s_dot, s_ddot)
            
            accel = this.a_fifth * this.fifth_phi_ddot(s, s_dot, s_ddot);
            
        end


    end

    
    properties (Access = public)
        
        n_dof
        a_fifth
        
    end
    
end