function f = forcingTerm(this, x)

    Psi = this.kernelFunction(x);
    f = dot(Psi,this.w) / (sum(Psi)+this.zero_tol); % add 'zero_tol' to avoid numerical issues

end