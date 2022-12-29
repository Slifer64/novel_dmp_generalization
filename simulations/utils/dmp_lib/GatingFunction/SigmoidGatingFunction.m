%% Sigmoid Gating Function class
%  Implements a sigmoidal gating function, u=f(x), x:[0 1]->u:[u0 u_end],
%  where u0 is the initial and u_end the final value.
%  The output of the gating function is:
%     u = u0 * ( 1 / (1 + exp(-a_u*(x-c)) ) );
%    du = -a_u*u0 * ( exp(-a_u*(x-c)) / (1 + exp(-a_u*(x-c)) )^2 );
%

classdef SigmoidGatingFunction < GatingFunction
   
   methods (Access = public)
       
      %% Sigmoid Gating Function Constructor.
      %  @param[in] u0: Initial value of the gating function (optional, default = 1.0).
      %  @param[in] u_end: Final value of the gating function (optional, default = 0.99).
      %  @param[out] this: Gating function object.
      function this = SigmoidGatingFunction(u0, u_end)

          if (nargin < 1), u0 = 1.0; end
          if (nargin < 2), u_end = 0.99; end
          
          this@GatingFunction();

          this.u0 = u0;
          this.u_end = u_end;
          this.setSteepness(700);
          this.calcCenter();

      end


      %% Returns the gating function's output for the specified timestamps.
      %  @param[in] x: Vector of timestamps.
      %  @param[out] u: Vector of values of the gating function's output.
      function u = getOutput(this, x)

          exp_t = exp((this.a_u)*(x-this.c));
          u = this.u0 * 1.0 ./ (1.0 + exp_t);

      end

      
      %% Returns the gating function's derivated output for the specified timestamps.
      %  @param[in] x: Vector of timestamps.
      %  @param[out] u: Vector of values of the gating function's derivated output.
      function du = getOutputDot(this, x)

          exp_t = exp(this.a_u*(x-this.c));
          du = -this.u0 * (this.a_u) * exp_t ./ (1.0 + exp_t).^2;

      end
      
      
      %% Returns the partial derivative of the gating output wrt 1/tau
      %  @param[in] t: timestamp
      %  @param[in] x: phase variable
      %  @param[out] u: partial derivative of the gating wrt 1/tau.
      function u = getPartDev_1oTau(this, t, x)

          h = this.getOutput(x);
          
          u = -this.a_u * t * h * (1-h);

      end
      
      
      %% Sets the sigmoid's steepness.
      %  @param[in] a_u: steepness.
      function setSteepness(this, a_u)
          
          this.a_u = a_u;
          this.calcCenter();
          
      end

   end
   
   methods (Access = protected)
        
       function calcCenter(this)

          this.c = 1.0 - (1.0/this.a_u)*log((this.u0-this.u_end)/this.u_end);
          
       end
          
       
   end
   
   properties (Access = private)
        
       u_end % final value of the gating function (at x=1)
       u0 % initial value of the gating function (at x=0)
       a_u % the rate of evolution of the gating function
       c % center of the exponential in the sigmoid
       
   end
   
end
