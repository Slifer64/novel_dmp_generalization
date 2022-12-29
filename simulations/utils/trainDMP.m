function gmp = trainDMP(Timed, Pd_data)

    % includes...
    import_gmp_lib();

    %% =============  Set GMP params  =============
    train_method = 'LS'; % {'LS', 'LWR'}
    N_kernels = 25;  % number or kernels (Gaussians)
    kernels_std_scaling = 1.5; % scaling of the width of each Gaussian. The greater 
                               % the scaling the more the overlapping). A vaule 
                               % in [1 1.5] is usually good.

    n_dof = size(Pd_data, 1);

    %% =============  Create/Train GMP  =============
    gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
    t_start = tic;
    sd_data = Timed/Timed(end);
    offline_train_mse = gmp.train(train_method, sd_data, Pd_data);
    offline_train_mse
    toc(t_start)
    
end