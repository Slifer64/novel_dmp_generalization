function setStds(this, kernelStdScaling)
            
    if (nargin < 2), kernelStdScaling=1.0; end

    this.h = 1./(kernelStdScaling*(this.c(2:end)-this.c(1:end-1))).^2;
    this.h = [this.h; this.h(end)];

end