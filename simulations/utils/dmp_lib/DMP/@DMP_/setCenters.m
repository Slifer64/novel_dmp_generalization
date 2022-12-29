function setCenters(this)

    t = ((1:this.N_kernels)-1)/(this.N_kernels-1);
    x = this.phase(t*this.getTau());
    this.c = x(1,:)';

end