function ddy = calcYddot(this, x, y, dy, g, tau_dot, yc, zc, yc_dot)

    if (nargin < 6), tau_dot = 0; end
    if (nargin < 7), yc = 0; end
    if (nargin < 8), zc = 0; end
    if (nargin < 9), yc_dot = 0; end

    tau = this.getTau();
    z = dy*tau - yc;
    
    shape_attr = this.shapeAttractor(x, g);
    goal_attr = this.goalAttractor(x, y, z, g);
    dz = ( goal_attr + shape_attr + zc) / tau;

    ddy = (dz + yc_dot - tau_dot*dy)/tau;

end