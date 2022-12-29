function cp_obj = deepCopy(this)
            
    % Make a shallow copy of all properties
    cp_obj = this.copy();
    % Make a deep copy of the pointers
    cp_obj.can_clock_ptr = this.can_clock_ptr.copy();
    cp_obj.shape_attr_gating_ptr = this.shape_attr_gating_ptr.copy();

end