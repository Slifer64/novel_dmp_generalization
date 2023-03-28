classdef VMP < matlab.mixin.Copyable
    
    methods (Access = public)
        
        function this = VMP(gmp)
            
            this.vp_map = containers.Map();
            
            this.K = 300;
            this.D = 2*sqrt(this.K + 10);
            
            this.gmp = gmp.deepCopy();
            this.n_dof = this.gmp.numOfDoFs();
            this.gmp.setScaleMethod(TrajScale_None(gmp.numOfDoFs()));
            
            this.pol5_mp = FifthOrderMP(this.n_dof);
            
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
            
            this.pol5_mp.update(s1, h1, h1_dot, h1_ddot, s2, h2, h2_dot, h2_ddot, this.s_dot, this.s_ddot);
            
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
                ydot_via = [ydot_via vp_groups{k}.vel];
                yddot_via = [yddot_via vp_groups{k}.accel];
            end
            
            [~, ind] = sort(s_via);
            s_via = s_via(ind);
            y_via = y_via(:, ind);
            ydot_via = ydot_via(:, ind);
            yddot_via = yddot_via(:, ind);
            
%             s_via
%             pause
            
            i2 = length(s_via);
            for j=2:length(s_via)
               if (s_via(j) > s)
                   i2 = j;
                   break;
               end
            end
            i1 = i2-1;

            s1 = s_via(i1);
            s2 = s_via(i2);
            
            disp('-----------')
            fprintf('s=%.3f , s1=%.3f, s2=%.3f \n', s, s1, s2);
%             pause
            
            yv1 = y_via(:,i1);
            yv1_dot = ydot_via(:,i1);
            yv1_ddot = yddot_via(:,i1);
            
            s1 = s;
            if (abs(s1 - s2) < 5e-3), s1 = s2 - 5e-3; end
            yv1 = this.getRefPos(s1);
            yv1_dot = 1*this.getRefVel(s1, this.s_dot);
            yv1_ddot = 1*this.getRefAccel(s1, this.s_dot, this.s_ddot);

            yv2 = y_via(:,i2);
            yv2_dot = 1*this.gmp.getYdDot(s2, this.s_dot);
            yv2_ddot = 1*this.gmp.getYdDDot(s2, this.s_dot, this.s_ddot);
            
%             if (s2 > 0.5 && s2 < 0.95)
%                 yv2_dot = norm(yv2_dot)*[0; -1];
%                 yv2_ddot = norm(yv2_ddot)*[0; 1];
%             end
            
%             yv2_dot = ydot_via(:,i2);
%             yv2_ddot = yddot_via(:,i2);
            
%             if any(isnan(ydot_via(:,i1)))
%                 ydot_via(:,i1) = this.getRefVel(s1, this.s_dot);
%             end
%             if any(isnan(yddot_via(:,i1)))
%                 yddot_via(:,i1) = this.getRefAccel(s1, this.s_dot, this.s_ddot);
%             end

            h1 = yv1 - this.gmp.getYd(s1);
            h1_dot = yv1_dot - this.gmp.getYdDot(s1, this.s_dot);
            h1_ddot = yv1_ddot - this.gmp.getYdDDot(s1, this.s_dot, this.s_ddot);
            
%             h1 = this.pol5_mp.getRefPos(s1);
%             h1_dot = this.pol5_mp.getRefVel(s1, this.s_dot);
%             h1_ddot = this.pol5_mp.getRefAccel(s1, this.s_dot, this.s_ddot);
            
            h2 = yv2 - this.gmp.getYd(s2);
            h2_dot = yv2_dot - this.gmp.getYdDot(s2, this.s_dot);
            h2_ddot = yv2_ddot - this.gmp.getYdDDot(s2, this.s_dot, this.s_ddot);
            
            this.pol5_mp.update(s1, h1, h1_dot, h1_ddot, s2, h2, h2_dot, h2_ddot, this.s_dot, this.s_ddot);

        end
        
        
        function vp_s = updateViapoints(this, s, vp_pos, vp_name)
            
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
            end
            this.vp_map(vp_name) = struct('s',vp_s, 'pos',vp_pos, 'vel',vp_vel, 'accel',vp_accel);
            
            this.update(this.g, s0, this.s_dot, this.getRefPos(s), this.getRefVel(s, this.s_dot));
        end
        
        function p_ref = getRefPos(this, s)

            p_ref = this.gmp.getYd(s) + this.pol5_mp.getRefPos(s);
            
        end
        
        function v_ref = getRefVel(this, s, s_dot)

            v_ref = this.gmp.getYdDot(s, s_dot) + this.pol5_mp.getRefVel(s, s_dot);
            
        end
        
        function a_ref = getRefAccel(this, s, s_dot, s_ddot)
            
            a_ref = this.gmp.getYdDDot(s, s_dot, s_ddot) + this.pol5_mp.getRefAccel(s, s_dot, s_ddot);
            
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
                yd0 = this.gmp.getYd(0);
                gd = this.gmp.getYd(1);
                Ks = diag((this.g - this.y0) ./ (gd - yd0));
                yd = this.gmp.getYd(s);
                yd = Ks*(yd - yd0) + this.y0;
                dist = norm(pos - yd);
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
        
        pol5_mp
        
    end
    
end