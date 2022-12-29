function psi = kernelFunction(this,x)

    n = length(x);
    psi = zeros(this.N_kernels, n);

    for j=1:n
        psi(:,j) = exp(-this.h.*((x(j)-this.c).^2));
    end 

end