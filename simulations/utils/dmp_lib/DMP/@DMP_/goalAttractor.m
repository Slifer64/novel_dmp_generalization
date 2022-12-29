function goal_attr = goalAttractor(this, x, y, z, g)

    g_attr_gating = this.goalAttrGating(x);
    goal_attr = g_attr_gating * this.a_z*(this.b_z*(g-y)-z);

end