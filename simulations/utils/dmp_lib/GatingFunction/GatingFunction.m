%% Linear Gating Function class
%  Implements a linear gating function, u=f(x), x:[0 1]->u:[u0 u_end],
%  where u0 is the initial and u_end the final value.
%  The output of the gating function is:
%     u = u0 - a_u*x;
%    du = -a_u;
%

classdef GatingFunction < matlab.mixin.Copyable

   methods
      %% Linear Gating Function Constructor.
      %  @param[in] u0: Initial value of the gating function (optional, default = 1.0).
      %  @param[in] u_end: Final value of the gating function (optional, default = 0.005).
      %  @param[out] gating_fun: Gating function object.
      function this = GatingFunction()

      end
      
   end

   methods (Abstract)
       
      %% Returns the gating function's output for the specified timestamps.
      %  @param[in] x: Vector of timestamps.
      %  @param[out] u: Vector of values of the gating function's output.
      u = getOutput(this, x)

      
      %% Returns the gating function's derivated output for the specified timestamps.
      %  @param[in] x: Vector of timestamps.
      %  @param[out] u: Vector of values of the gating function's derivated output.
      du = getOutputDot(gating_fun, x)
      
      
      %% Returns the partial derivative of the gating output wrt 1/tau
      %  @param[in] t: timestamp
      %  @param[in] x: phase variable
      %  @param[out] u: partial derivative of the gating wrt 1/tau.
      u = getPartDev_1oTau(this, t, x)

   end
end
