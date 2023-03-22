% Generic DMP implementation

classdef DMP_ < matlab.mixin.Copyable
    
    properties (Access = public)
        
        K % DMP stiffness
        D % DMP damping
        
    end
    
    methods (Access = public)
        
        function this = DMP_(n_dofs, N_kernels, spatial_scaling, kernel_fun, gating)
            
            this.K = 400;
            this.D = 2*sqrt(this.K + 10);
            this.W = zeros(n_dofs, N_kernels);
            this.c = linspace(0, 1, N_kernels);
            this.h_scale = 1.0;
            
            this.setKernelFun(kernel_fun);
            this.setGatingFun(gating);
            this.setSpatialScaling(spatial_scaling);

        end

        function train_mse = train(this, Time, y_data, dy_data, ddy_data)
           
            n_dofs = size(y_data, 1);
            if (n_dofs ~= this.numOfDoFs()), error('N_dofs mismatch!'); end
            
            n_data = length(Time);
            
            dTime = diff(Time);
            if (nargin < 4), dy_data = [diff(y_data, 1, 2) ./ dTime, zeros(n_dofs,1)]; end
            if (nargin < 5), ddy_data = [diff(dy_data, 1, 2) ./ dTime, zeros(n_dofs,1)]; end
            
            g = y_data(:, end);
            y0 = y_data(:, 1);
            Time = Time - Time(1);
            tau = Time(end);
            s_data = Time / Time(end);
            
            Fd = this.calcFd(s_data, y_data, dy_data, ddy_data, g, y0, tau);
            
            Phi = zeros(this.numOfKernels(), n_data);
            for j=1:n_data, Phi(:,j) = this.phi_fun(s_data(j)); end
            
            this.W = Fd / Phi;
            
%             Fd_hat = this.W * Phi;
%             figure;
%             n_dofs = this.numOfDoFs();
%             for i=1:n_dofs
%                subplot(n_dofs, 1, i); hold on;
%                plot(s_data, Fd(i,:), 'color','blue', 'LineWidth',2);
%                plot(s_data, Fd_hat(i,:), 'color','magenta');
%             end
%             Fd_mse = mean((this.W*Phi - Fd).^2, 2)
            
            this.sp_scaling.setNominalStartFinalPos(g, y0);
            this.sp_scaling.setNewStartFinalPos(g, y0);
            this.g = g;
            this.y0 = y0;
            
            [Time_hat, y_hat_data, ~, ~] = this.generate_trajectory(y0, g, Time(end), 0.005);
            y_hat_data = interp1(Time_hat, y_hat_data', Time, 'linear')';
            train_mse = mean((y_hat_data - y_data).^2, 2);
            
            %train_mse = mean((this.W*Phi - Fd).^2, 2);
            
        end
        
        function n_ker = numOfKernels(this)
           
            n_ker = length(this.c);
            
        end
        
        function n_dofs = numOfDoFs(this)
           
            n_dofs = size(this.W, 1);
            
        end
        
        function setY0(this, y0)

            this.sp_scaling.setNewStartFinalPos(this.g, y0);
            this.y0 = y0;
            
        end
        
        function setGoal(this, g)
            
            this.sp_scaling.setNewStartFinalPos(g, this.y0);
            this.g = g;
            
        end

        function [Time, y_data, dy_data, ddy_data] = generate_trajectory(this, y0, g, Tf, dt)
            
            Time = 0:dt:Tf;
            if (Time(end) < Tf), Time = [Time Time(end)+dt]; end
            y_data = zeros(this.numOfDoFs(), length(Time));
            dy_data = zeros(size(y_data));
            ddy_data = zeros(size(y_data));
            
            % store previous values
            g_prev = this.g;
            y0_prev = this.y0;
            
            % set new init/final pos
            this.setY0(y0);
            this.setGoal(g);
            
            tau = Tf;
            s_data = Time / Time(end);
            s_dot = 1/tau;
            s_ddot = 0;
            
            y = y0;
            dy = zeros(size(y));
            
            for j=1:length(Time)
                
                s = s_data(j);
                ddy = this.goal_attractor(y, dy, tau) + this.shape_attractor(s, s_dot, s_ddot, tau);
                
                y_data(:, j) = y;
                dy_data(:, j) = dy;
                ddy_data(:, j) = ddy;
                
                y = y + dy*dt;
                dy = dy + ddy*dt;
                
            end
            
            % restore previous values
            this.setY0(y0_prev);
            this.setGoal(g_prev);
        end
        
        function plotKernels(this, s_data)
            
            if (nargin < 2), s_data = linspace(0, 1, 400); end
            
            n_data = length(s_data);
            psi_data = zeros(this.numOfKernels(), n_data);
            for j=1:n_data, psi_data(:, j) = this.kernel_fun(s_data(j)); end

            phi_data = psi_data ./ (sum(psi_data, 1) + 1e-16);

            figure;
            subplot(2, 1, 1);
            plot(s_data, psi_data')
            ylabel('$\psi$', 'interpreter','latex', 'fontsize',14);
            subplot(2, 1, 2);
            plot(s_data, phi_data')
            ylabel('$\phi = \psi / \sum_{k=1}^K \psi_k$', 'interpreter','latex', 'fontsize',14);
            
        end
    end
    
    %% ============  Abstract  ==============
    methods (Abstract, Access = public)
        
        f_g = goal_attractor(this, y, y_dot, tau)
        
        f_s = shape_attractor(this, s, s_dot, s_ddot, tau)

    end
    
    methods (Abstract, Access = protected)

        Fd = calcFd(this, s, y, dy, ddy, g, y0, tau)

    end

    %% ============  Protected methods ==============
    
    methods (Access = protected)
        
        function setKernelFun(this, kernel_fun)
            
            if strcmpi(kernel_fun, 'gaussian')
                this.h = ones(size(this.c)) * this.h_scale / (this.c(2) - this.c(1))^2;
                this.kernel_fun = @(s)this.gaussian_kernel(s, this.c, this.h);
            elseif strcmpi(kernel_fun, 'mollifier')
                this.h = ones(size(this.c)) * this.h_scale / (this.c(2) - this.c(1));
                this.kernel_fun = @(s)this.mollifier_kernel(s, this.c, this.h);
            else
                error(['Unsupported kernel function: "' kernel_fun '"']);
            end
            
        end
        
        function setGatingFun(this, gating)
            
            if strcmpi(gating, 'exp')
                this.gating_fun = @this.exp_gating;
            elseif strcmpi(gating, 'linear')
                this.gating_fun = @this.linear_gating;
            else
                error(['Unsupported gating: "' gating '"']);
            end
            
        end
        
        function setSpatialScaling(this, spatial_scaling)
            
            n_dofs = this.numOfDoFs();
            
            if strcmpi(spatial_scaling, 'none')
                this.sp_scaling = TrajScale_None(n_dofs);
            elseif strcmpi(spatial_scaling, 'diag')
                this.sp_scaling = TrajScale_Prop(n_dofs);
            elseif strcmpi(spatial_scaling, 'roto-dial')
                this.sp_scaling = TrajScale_roto_dial(n_dofs);
            elseif strcmpi(spatial_scaling, 'rot-min')
                if (n_dofs ~= 3), error('"rot-min" scaling applies only for 3-DoFs'); end
                this.sp_scaling = TrajScale_Rot_min();
            elseif strcmpi(spatial_scaling, 'rot-wb')
                if (n_dofs ~= 3), error('"rot-wb" scaling applies only for 3-DoFs'); end
                this.sp_scaling = TrajScale_Rot_wb();
            else
                error(['Unsupported spatial scaling: "' spatial_scaling '"']);
            end
            
        end
        
        function phi = phi_fun(this, s)
            
            psi = this.kernel_fun(s);
            phi = psi / (sum(psi) + 1e-16);
            
        end
        
        %% S_scale * W * Phi * gat(s)
        function f_term = forcing_term(this, s)
            
            phi = this.phi_fun(s);
            f_term = this.sp_scaling.getScaling() * this.W*phi*this.gating_fun(s);
            
        end

    end
    
    %% ============  Protected properties  ==============
    
    properties (Access = protected)
        
        W % DMP weights
        sp_scaling % object that returns the spatial scaling matrix
        gating_fun % gating function pointer
        kernel_fun % kernel function pointer
        g % target position
        y0 % initial position
        
        c % kernel centers
        h % kernel inv widths
        h_scale % scaling of the kernel widths
        
    end
    
    %% ============  Static  ==============
    
    methods (Static, Access = protected)
        
        %% Gating functions
        
        function u = exp_gating(s)
        
            u0 = 1.0;
            u_end = 0.005;
            a = -log(u_end/u0);
            u = u0*exp(-a*s);
            
        end
        
        function u = linear_gating(s)
        
            u0 = 1.0;
            u_end = 0.005;
            a = u0 - u_end;
            u = u0 - a*s;
            u(u<0) = 0;
            
        end
        
        %% Kernel functions
        
        function psi = gaussian_kernel(s, c, h)

            psi = exp(-h.*((s-c).^2));
            psi = psi(:);

        end
        
        function psi = mollifier_kernel(s, c, h)
            
            x = abs(h .* (s - c));
            psi = zeros(size(c));
            i = find(x < 1.0);
            psi(i) = exp(-1.0 ./ (1.0 - x(i).^2));
            psi = psi(:);
            
        end

    end

    
end