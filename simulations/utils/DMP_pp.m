classdef DMP_pp < matlab.mixin.Copyable
    
    methods (Access = public)
        
        function this = DMP_pp(gmp, traj_scale)
            
            if (nargin < 2), traj_scale=TrajScale_None(gmp.numOfDoFs()); end
            
            this.vp_map = containers.Map();
            
            this.K = 300;
            this.D = 2*sqrt(this.K + 10);
            
            this.gmp = gmp.deepCopy();
            this.n_dof = this.gmp.numOfDoFs();
            this.gmp.setScaleMethod(traj_scale);
            
            this.gmp_up = GMP_Update(this.gmp);
            this.W0 = gmp.W;
            
            this.prev_state = struct('initialized',false, 'target',[], 'sf',nan, 's',nan, 's_dot',nan);
            
            e1 = 1e-0;
            this.r0 = e1*[1e-9, 1e-7, 1e-7];
            this.rf = e1*[1e-9, 1e-7, 1e-7];
            this.r1 = e1*[1e-6 , 1e-6, 1e-4];
            this.rv = 1e-7;

            this.setAdaptToRobot(false);
            this.setRecursiveUpdate(true);
            
            this.reset();
            
        end
        
        function setAdaptToRobot(this, set)
           
            this.adapt_to_robot = set;
            
        end
        
        function setRecursiveUpdate(this, set)
            
            this.recursive_update = set;
            
        end
        
        function reset(this)
            
            this.gmp.W = this.W0;
            
            this.gmp_up.syncUpdate(false);
            this.gmp_up.recursiveUpdate(this.recursive_update);
            this.setOptMetric('accel');
            
            this.prev_state.initialized = false;
        end
        
        function setOptMetric(this, metric)
           
            if (strcmpi(metric,'pos'))
                this.gmp_up.initSigmaWfromMsr(0:0.01:1, 1e-1);
            elseif (strcmpi(metric,'vel'))
                this.gmp_up.initSigmaWfromVelMsr(0:0.01:1, 1, 1e-1);
            elseif (strcmpi(metric,'accel'))
                this.gmp_up.initSigmaWfromAccelMsr(0:0.01:1, 1, 1e-1);
            elseif (strcmpi(metric,'weights'))
                this.gmp_up.initSigmaw();  
            else
                error(['Unsupported metric "' metric '"']);
            end
            
        end
        
        function init(this, s_start, y0, yg, Tf)
            
            s0 = 0;
            sf = 1;
            s_dot = 1 / Tf;
            s_ddot = 0;
            
            this.y0 = y0;

            O_zeros = zeros(this.n_dof, 1);
            
            this.gmp.setY0(y0);
            this.gmp.setGoal(yg);

            % initial state constraints: pos, vel, accel
            this.gmp_up.updatePos(s0, this.y0, this.r0(1));
            this.gmp_up.updateVel(s0, s_dot, O_zeros, this.r0(2));
            this.gmp_up.updateAccel(s0, s_dot, s_ddot, O_zeros, this.r0(3));

            % final state constraints: pos, vel, accel
            this.gmp_up.updateVel(sf, s_dot, O_zeros, this.rf(2));
            this.gmp_up.updateAccel(sf, s_dot, s_ddot, O_zeros, this.rf(3));

%             this.gmp_up.updateNow();
%             this.Sigma_w1 = this.gmp_up.getSigmaW();
%             this.W1 = this.gmp.W;
            
            this.gmp_up.updatePos(sf, yg, this.rf(1));
            this.gmp_up.updateNow();

            this.prev_state.target = yg;
            this.prev_state.sf = sf;
            this.prev_state.s = s_start;
            this.prev_state.s_dot = s_dot;
            this.prev_state.initialized = true;
            
        end
            
        function update(this, yg, s, s_dot, y, y_dot)
            
            if (~this.prev_state.initialized), error('You need to call DMP_pp::init first!'); end

            s0 = 0;
            sf = 1; %max([s, 1.0]);
            s_ddot = 0;
            O_zeros = zeros(this.n_dof, 1);

            y1 = y;
            y1_dot = y_dot;
            y1_ddot = this.gmp.getYdDDot(this.prev_state.s, this.prev_state.s_dot, s_ddot);

            yd = this.gmp.getYd(s);
            yd_dot = this.gmp.getYdDot(s, s_dot);

            if (~this.adapt_to_robot)
                y1 = yd;
                y1_dot = yd_dot;
            else
                % Add a small tolerance to avoid numerical integration drift
                if (norm(y1 - yd) < 5e-4), y1 = yd; end
                if (norm(y1_dot - yd_dot) < 1e-3), y1_dot = yd_dot; end
            end

            % ---------- downdate --------------
            this.gmp_up.updatePos(this.prev_state.sf, this.prev_state.target, -this.rf(1));
            this.gmp_up.updateNow();
%             this.gmp_up.setSigmaW(this.Sigma_w1);
%             this.gmp.W = this.W1;

            % ---------- update --------------
            this.gmp.setGoal(yg);
            
            % current state constraints: pos, vel, accel
            this.gmp_up.updatePos(s, y1, this.r1(1));
            this.gmp_up.updateVel(s, s_dot, y1_dot, this.r1(2));
            this.gmp_up.updateAccel(this.prev_state.s, this.prev_state.s_dot, s_ddot, y1_ddot, this.r1(3));
            
            if (~this.recursive_update)
                this.gmp_up.updatePos(s0, this.y0, this.r0(1));
                this.gmp_up.updateVel(s0, s_dot, O_zeros, this.r0(2));
                this.gmp_up.updateAccel(s0, s_dot, s_ddot, O_zeros, this.r0(3));
                this.gmp_up.updateVel(sf, s_dot, O_zeros, this.rf(2));
                this.gmp_up.updateAccel(sf, s_dot, s_ddot, O_zeros, this.rf(3));
            end
            
%             this.gmp_up.updateNow();
%             this.Sigma_w1 = this.gmp_up.getSigmaW();
%             this.W1 = this.gmp.W;
            
            % final state constraints: pos, vel, accel
            this.gmp_up.updatePos(sf, yg, this.rf(1));
            this.gmp_up.updateNow();
            
            this.prev_state.target = yg;
            this.prev_state.sf = sf;
            this.prev_state.s = s;
            this.prev_state.s_dot = s_dot;
            this.prev_state.initialized = true;
        end
        
        function updateViapoints(this, s, via_points, vp_name)
            
            if this.vp_map.isKey(vp_name)
                % remove the previous via-points
                prev_vp = this.vp_map(vp_name);
                for j=1:length(prev_vp.s)
                    this.gmp_up.updatePos(prev_vp.s(j), prev_vp.pos(:,j), -this.rv);
                end
                this.gmp_up.updateNow();
            end
            
            n_points = size(via_points,2);
            vp_s = zeros(1, n_points);
            for j=1:n_points
                s = this.findClosest(s, 1.0, 60, via_points(:,j));
                vp_s(j) = s;
            end
            this.vp_map(vp_name) = struct('s',vp_s, 'pos',via_points);
            
%             dmp_pos = zeros(this.n_dof ,n_points);
%             for j=1:n_points
%                dmp_pos(:,j) = this.getRefPos(vp_s(j)); 
%             end
%             vp_s
%             via_points
%             dmp_pos
            
            for j=1:size(via_points,2)
                this.gmp_up.updatePos(vp_s(j), via_points(:,j), this.rv);
            end
            this.gmp_up.updateNow();

        end
        
        
        function s = updateViapoint(this, s_v, pos, search_closest, reverse)
        
            if (nargin < 4), search_closest=true; end
            if (nargin < 5), reverse=false; end
            
            if (search_closest)
                if (reverse)
                    s = this.findClosest(s_v, 0, 60, pos);
                    ds = s*0.1 + 1e-9;
                    s = this.findClosest(min([s+ds, 1.0]), max([s_v, s-ds]), 20, pos);
                else
                    s = this.findClosest(s_v, 1.0, 60, pos);
                    ds = (1-s)*0.1 + 1e-9;
                    s = this.findClosest(max([s_v, s-ds]), min([s+ds, 1.0]), 20, pos);
                end
            else
                s = s_v;
            end
            
            this.gmp_up.updatePos(s, pos, this.rv);
            this.gmp_up.updateNow();

        end
        
        function p_ref = getRefPos(this, s)
            
            p_ref = this.gmp.getYd(s);
            %if (s > 1), p_ref = this.prev_state.target; end
            
        end
        
        function v_ref = getRefVel(this, s, s_dot)
            
            v_ref = this.gmp.getYdDot(s, s_dot);
            %if (s > 1), v_ref = 0*v_ref; end
            
        end
        
        function a_ref = getRefAccel(this, s, s_dot, s_ddot)
            
            a_ref = this.gmp.getYdDDot(s, s_dot, s_ddot);
            %if (s > 1), a_ref = 0*a_ref; end
            
        end
        
        function f_g = goal_attractor(this, y, y_dot, tau)
           
            f_g = this.K*(this.gmp.getGoal() - y) - this.D*y_dot;
            
        end
        
        function f_s = shape_attractor(this, s, s_dot, s_ddot, tau)
            
            %this.K*(this.getRefPos(s) - this.gmp.getGoal())*(s<1)
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
                dist = norm(pos - this.gmp.getYd(s));
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

        r0
        rf
        r1
        rv
        
        gmp
        
    end
    
    properties (Access = protected)
        
        vp_map
        
        y
        y_dot
        
        gmp_up
        
        
        W0
        
        Sigma_w1
        W1
        
        prev_state
        
        n_dof
        
        y0
        
        adapt_to_robot
        recursive_update
        
    end
    
end