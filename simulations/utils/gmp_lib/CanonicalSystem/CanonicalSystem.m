classdef CanonicalSystem < handle
    
    %% ============================
    %% ========  public  ==========
    %% ============================
    
    methods (Access = public)
        
        %% Constructor
        % param[in] T: time duration for the phase to go from 0 to 1.
        % param[in] Ds: rate of tracking the time duration. (default = 30).
        %               Set to inf to allow instant (discontinuous changes)
        function this = CanonicalSystem(T, Ds)
            
            if (nargin < 2), Ds = 30; end
            
            this.s0 = 0;
            this.sf = 1;
            this.Ds = Ds;

            this.sd_dot = 1/T;
            
            this.s = 0;
            this.s_dot = this.sd_dot;

        end
        
        %% Integrate the phase.
        % @param[in] t0: initial time instant.
        % @param[in] tf: final time instant.
        function integrate(this, t0, tf)
        
            if (isinf(this.Ds))
                
                this.s_dot = this.sd_dot;
                this.s = this.s + this.s_dot*(tf - t0);
                
                if (this.s > this.sf), this.s = this.sf; 
                elseif (this.s < this.s0) this.s = this.s0;
                end
                
                return;
                
            end
                
            % ode
            ode_fun = @(t, state) [state(2); this.getPhaseDDot(state(1), state(2))];
            [~,state_data] = ode45(@(t,state)ode_fun(t,state), [t0 tf], [this.s; this.s_dot]);
            this.s = state_data(end,1);
            this.s_dot = state_data(end,2);
            
        end
        
        %% Set the time duration for the phase to reach 1, given the current time instant.
        % param[in] T: time duration
        % param[in] t: current time instant (default = 0).
        function setDuration(this, T, t)
        
            if (nargin < 3), t = 0; end
            
            if (T < t), error('The current time has already exceeded the duration'); end
                
            this.sd_dot = (this.sf - this.s) / (T - t);
            
        end
        
        %% Set the remaing time duration for the phase to reach 1.
        % param[in] T: reaming time duration
        function setRemainingDuration(this, T)
                
            this.sd_dot = (this.sf - this.s) / T;
            
        end
           
        %% Set the duration
        % param[in] s: phase.
        % param[in] s_dot: phase 1st time derivative.
        % return: the 2nd time derivative of the phase
        function s_ddot = getPhaseDDot(this, s, s_dot)
            
            if (nargin < 3), s_dot = this.s_dot; end
            if (nargin < 2), s = this.s; end
            
            if (s>=this.s0 && s < this.sf), s_ddot = -this.Ds*(s_dot - this.sd_dot);
            elseif (s>=1),                  s_ddot = -400*s_dot - 1000*(s-this.sf);
            else,                           s_ddot = -400*s_dot - 1000*(s-this.s0);
            end
            
        end
        
        %% Set state values to their initial values.
        function reset(this)
                
            this.s = this.s0;
            this.s_dot = this.sd_dot;
            
        end
           
           
    end
    
    properties (Access = public)

        s % phase
        s_dot % phase 1st time derivative
        sd_dot % desired speed of the phase evolution

    end
    
    %% ===============================
    %% ========  protected  ==========
    %% ===============================
    
    properties (Access = protected)

        s0 % initial value of the phase
        sf % final value of the phase
        
        Ds % rate of tracking the desired speed of the phase evolution

    end
    
end