classdef VMP < matlab.mixin.Copyable
    
    methods (Access = public)
        
        function this = VMP(gmp)
            
            this.vp_map = containers.Map();
            
            this.K = 300;
            this.D = 2*sqrt(this.K + 10);
            
            this.gmp = gmp.deepCopy();
            this.n_dof = this.gmp.numOfDoFs();
            this.gmp.setScaleMethod(TrajScale_None(gmp.numOfDoFs()));
            
            this.a_fifth = zeros(this.n_dof, 6);
            
        end

        
        function init(this, s_start, y0, yg, Tf)
            
            this.y0 = y0;
            this.g = yg;
            this.Tf = Tf;
            this.s_dot = 1/this.Tf;
            this.s_ddot = 0;

            O_nof = zeros(this.n_dof, 1);

            s1 = s_start;
            s2 = 1.0;

            h1 = this.y0 - this.gmp.getYd(s1);
            h1_dot = O_nof - this.gmp.getYdDot(s1, this.s_dot);
            h1_ddot = O_nof - this.gmp.getYdDDot(s1, this.s_dot, this.s_ddot);

            h2 = this.g - this.gmp.getYd(s2);
            h2_dot = O_nof - this.gmp.getYdDot(s2, this.s_dot);
            h2_ddot = O_nof - this.gmp.getYdDDot(s2, this.s_dot, this.s_ddot);
            
            this.solve_a_fifth(s1, h1, h1_dot, h1_ddot, s2, h2, h2_dot, h2_ddot);
            
        end
        
        function solve_a_fifth(this, s1, h1, h1_dot, h1_ddot, s2, h2, h2_dot, h2_ddot)
            
            A1 = this.fifth_pol_mat(s1, this.s_dot, this.s_ddot, 1, ~any(isnan(h1_dot(:))), ~any(isnan(h1_ddot(:))));
            A2 = this.fifth_pol_mat(s2, this.s_dot, this.s_ddot, 1, ~any(isnan(h2_dot(:))), ~any(isnan(h2_ddot(:))));
            A = [A1 A2];
            b = h1;
            if ~any(isnan(h1_dot(:))), b = [b h1_dot]; end
            if ~any(isnan(h1_ddot(:))), b = [b h1_ddot]; end
            b = [b h2];
            if ~any(isnan(h2_dot(:))), b = [b h2_dot]; end
            if ~any(isnan(h2_ddot(:))), b = [b h2_ddot]; end
            
            err_0 = vecnorm(this.a_fifth * A - b, 2, 1);
            
            this.a_fifth = b / A;
            
            
            this.a_fifth = b / A;
            
            err_1 = vecnorm(this.a_fifth * A - b, 2, 1);
            
            disp('--------------------')
            err_0
            err_1
        end
        
        function A = fifth_pol_mat(this, s, s_dot, s_ddot, pos, vel, accel)
            
            A = [];
            
            if (pos)
                A = [A; [1  s    s^2   s^3     s^4      s^5]];
            end
            
            if (vel)
                phi = [0  1  2*s   3*s^2   4*s^3    5*s^4]*s_dot;
                A = [A; phi];
            end
            
            if (accel)
                phi = [0  0  2     6*s    12*s^2   20*s^3]*s_dot + [0  1  2*s   3*s^2   4*s^3    5*s^4]*s_ddot;
                A = [A; phi];
            end
            
            A = A';
            
        end
        
            
        function update(this, yg, s, s_dot, y, y_dot)
            
            this.g = yg;
            this.s_dot = s_dot;
            this.s_ddot = 0;
            
            s_via = [0 1];
            y_via = [this.y0 this.g];
            ydot_via = zeros(this.n_dof, 2);
            yddot_via = zeros(this.n_dof, 2);
            
            vp_groups = this.vp_map.values();
            
            for k=1:numel(vp_groups)
                s_via = [s_via vp_groups{k}.s];
                y_via = [y_via vp_groups{k}.pos];
                
%                 for j=1:length(vp_groups{k}.s)
%                     s_j = vp_groups{k}.s(j);
%                     dy_v = this.getRefVel(s_j, s_dot);
%                     ddy_v = this.getRefAccel(s_j, s_dot, s_ddot);
%                     ydot_via = [ydot_via dy_v];
%                     yddot_via = [yddot_via ddy_v];
%                 end
                
                ydot_via = [ydot_via vp_groups{k}.vel];
                yddot_via = [yddot_via vp_groups{k}.accel];
            end
            
            [~, ind] = sort(s_via);
            s_via = s_via(ind);
            y_via = y_via(:, ind);
            ydot_via = ydot_via(:, ind);
            yddot_via = yddot_via(:, ind);
            
            i1 = 1;
            i2 = length(s_via);
            for j=2:length(s_via)
               if (s_via(j) > s)
                   i1 = j-1;
                   i2 = j;
                   break;
               end
            end

            s1 = s_via(i1);
            s2 = s_via(i2);
            
%             fprintf('s=%.3f , s1=%.3f, s2=%.3f \n', s, s1, s2);
%             disp('-----------')

            h1 = y_via(:,i1) - this.gmp.getYd(s1);
            h1_dot = ydot_via(:,i1) - this.gmp.getYdDot(s1, this.s_dot);
            h1_ddot = yddot_via(:,i1) - this.gmp.getYdDDot(s1, this.s_dot, this.s_ddot);
            
%             h1 = this.getRefPos(s1);
%             h1_dot = this.getRefVel(s1, this.s_dot);
%             h1_ddot = this.getRefAccel(s1, this.s_dot, this.s_ddot);
            
%             h1 = this.getRefPos(s);
%             h1_dot = this.getRefVel(s, s_dot);
%             h1_ddot = this.getRefAccel(s, s_dot, s_ddot);
            
            h2 = y_via(:,i2) - this.gmp.getYd(s2);
            h2_dot = ydot_via(:,i2) - this.gmp.getYdDot(s2, this.s_dot);
            h2_ddot = yddot_via(:,i2) - this.gmp.getYdDDot(s2, this.s_dot, this.s_ddot);
            
            this.solve_a_fifth(s1, h1, h1_dot, h1_ddot, s2, h2, h2_dot, h2_ddot);

            fprintf('s=%.3f , s1=%.3f, s2=%.3f \n', s, s1, s2);
%             pause
            
%             y1_hat = this.getRefPos(s1);
%             y2_hat = this.getRefPos(s2);
%             
%             y1_err = norm(y_via(:,i1) - y1_hat)
%             y2_err = norm(y_via(:,i2) - y2_hat)
            
        end
        
        
        function updateViapoints(this, s, vp_pos, vp_name)
            
            s0 = s;
            
            n_points = size(vp_pos,2);
            vp_s = zeros(1, n_points);
            vp_vel = nan(this.n_dof, n_points);
            vp_accel = nan(this.n_dof, n_points);
            for j=1:n_points
                s = this.findClosest(s, 1.0, 60, vp_pos(:,j));
                vp_s(j) = s;
                vp_vel(:,j) = nan(this.n_dof,1);
                vp_accel(:,j) = nan(this.n_dof,1);
%                 vp_vel(:,j) = this.getRefVel(s, s_dot);
%                 vp_accel(:,j) = this.getRefAccel(s, s_dot, s_ddot);
%                 vp_vel(:,j) = this.gmp.getYdDot(s, s_dot);
%                 vp_accel(:,j) = this.gmp.getYdDDot(s, s_dot, s_ddot);
            end
            this.vp_map(vp_name) = struct('s',vp_s, 'pos',vp_pos, 'vel',vp_vel, 'accel',vp_accel);
            
            this.update(this.g, s0, this.s_dot, this.getRefPos(s), this.getRefVel(s, this.s_dot));
        end
        
        function p_ref = getRefPos(this, s)
            
            t_fifth = [1 s s^2 s^3 s^4 s^5]';
            p_ref = this.gmp.getYd(s) + this.a_fifth*t_fifth;
            
        end
        
        function v_ref = getRefVel(this, s, s_dot)
            
            t_fifth = [0 1 2*s 3*s^2 4*s^3 5*s^4]'*s_dot;
            v_ref = this.gmp.getYdDot(s, s_dot)  + this.a_fifth*t_fifth;
            
        end
        
        function a_ref = getRefAccel(this, s, s_dot, s_ddot)
            
            t_fifth = [0 0 2 6*s 12*s^2 20*s^3]'*s_dot + [0 1 2*s 3*s^2 4*s^3 5*s^4]'*s_ddot;
            a_ref = this.gmp.getYdDDot(s, s_dot, s_ddot)  + this.a_fifth*t_fifth;
            
        end
        
        function f_g = goal_attractor(this, y, y_dot, tau)
           
            f_g = this.K*(this.gmp.getGoal() - y) - this.D*y_dot;
            
        end
        
        function f_s = shape_attractor(this, s, s_dot, s_ddot, tau)
            
            f_s = this.K*(this.getRefPos(s) - this.gmp.getGoal()) + this.D*this.getRefVel(s, s_dot) + this.getRefAccel(s, s_dot, s_ddot);
            
        end
 
        function n_dof = numOfDoFs(this)
            n_dof = this.n_dof;
        end
        
        function [Time, y_data, dy_data, ddy_data] = generate_trajectory(this, y0, g, Tf, dt)
           
            Time = 0:dt:Tf;
            if (Time(end) < Tf), Time = [Time Time(end)+dt]; end
            
            this.init(0, y0, g, Tf);
            
            y_data = zeros(this.numOfDoFs(), length(Time));
            dy_data = zeros(size(y_data));
            ddy_data = zeros(size(y_data));

            s_data = Time / Time(end);
            s_dot = 1/Tf;
            s_ddot = 0;
            for j=1:length(Time)
                s = s_data(j);
                y_data(:, j) = this.getRefPos(s);
                dy_data(:, j) = this.getRefVel(s, s_dot);
                ddy_data(:, j) = this.getRefAccel(s, s_dot, s_ddot);
            end
            
        end
        
    end
    
    methods (Access = protected)
        
        function s1 = findClosest(this, s0, sf, n_points, pos)

            % choose whether to search forward or backwards in time
            if (s0 < sf), s_data = linspace(s0, sf, n_points);
            else, s_data = -linspace(-s0, -sf, n_points);
            end

            min_dist = 1e8;
            i_min = 0;
            for i=1:length(s_data)
                s = s_data(i);
                dist = norm(pos - this.getRefPos(s));
                if (dist < min_dist)
                    i_min = i;
                    min_dist = dist;
                end
            end
            
            s1 = s_data(i_min);
            
        end

    end
    
    properties (Access = public)
        
        K % DMP stiffness
        D % DMP damping
        
        gmp
        
        rv
        rf
        r1
        
        vp_map
        
    end
    
    properties (Access = protected)
        
        n_dof
        
        y0
        g
        
        Tf
        s_dot
        s_ddot
        
        a_fifth
        
    end
    
end