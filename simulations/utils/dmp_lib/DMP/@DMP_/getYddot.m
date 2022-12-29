function ddy = getYddot(this, tau_dot, yc_dot)
            
    if (nargin < 2), tau_dot = 0; end
    if (nargin < 3), yc_dot = 0; end
    ddy = (this.getZdot() + yc_dot - tau_dot*this.getYdot()) / this.getTau();

end