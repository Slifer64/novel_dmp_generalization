%% Linear Canonical Clock class
%  Implements a linear canonical clock, x = f(t), t:[0 tau] -> x:[x0 x_end].
%  The clock's evolution is defined as:
%     tau*dx = a_x
%  where x is the phase variable (clock's output), a_x a constant and
%  tau is a scaling factor defining the total time duration. The phase
%  variable can exceed x_end, which in turn means the time t exceeded tau.
%

classdef CanonicalClock < matlab.mixin.Copyable
   
   methods (Access = public)
      %% Linear Canonical Clock Constructor.
      %  @param[in] tau: Total time duration (used to set/scale the clocks time duration).
      function this = CanonicalClock(tau)

          if (nargin < 1)
            tau = 1.0;      
          end
          
          this.x0 = 0.0;
          this.x_end = 1.0;
          this.setTau(tau);
          this.a_x = this.x0 - this.x_end;
          
      end 

      
      %% Sets/scales the canonical clock's time duration.
      %  @param[in] tau: Total time duration (used to set/scale the clocks time duration).
      function setTau(this, tau)
          
          this.tau = tau;
          
      end
      
      %% Returns the canonical clock's time duration.
      %  @param[out] tau: Total time duration (used to set/scale the clocks time duration).
      function tau = getTau(this)
          
          tau = this.tau;
          
      end

      %% Returns the phase variable derivative of the canonical clock for 
      %% the specified phase variable value.
      %  @param[in] x: Vector of the phase variable values.
      %  @param[out] dx: Vector of phase variable derivatives.
      function dx = getPhaseDot(this, x)
          
           dx = -this.a_x*ones(size(x))/this.getTau();

      end
      
      %% Returns the output of the canonical clock for the specified
      %% timestamps.
      %  @param[in] t: Vector of timestamps.
      %  @param[out] x: Vector with the phase variable values at each timestamp in 't'.
      function x = getPhase(this, t)
          
          x = this.x0 - this.a_x*t/this.getTau();

      end

   end
   
   properties (Access = private)
       
       x0 % initial value of the phase variable
       x_end % final value of the phase variable
       a_x % the rate of evolution of the phase variable
       tau % total time duration
       
   end
   
end




