function update(this, x, y, z, g, y_c, z_c)

    if (nargin < 6), y_c=0; end
    if (nargin < 7), z_c=0; end

    tau = this.getTau();
    shape_attr = this.shapeAttractor(x, g);
    goal_attr = this.goalAttractor(x, y, z, g);

    this.dz = ( goal_attr + shape_attr + z_c) / tau;
    this.dy = ( z + y_c) / tau;
    this.dx = this.phaseDot(x);

end