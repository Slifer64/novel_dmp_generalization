function [train_error, F, Fd] = train(this, train_method, Time, yd_data, dyd_data, ddyd_data)
            
    n_data = length(Time);

    tau = Time(end);
    y0 = yd_data(:,1);
    g = yd_data(:,end);

    this.setTau(tau);
    this.setY0(y0);
    
    x = zeros(1, n_data);
    s = zeros(1, n_data);
    Fd = zeros(1,n_data);
    Psi = zeros(this.numOfKernels(), n_data);

    for i=1:n_data
        x(i) = this.phase(Time(i));
        s(i) = this.forcingTermScaling(g) * this.shapeAttrGating(x(i));
        Fd(i) = this.calcFd(x(i), yd_data(i), dyd_data(i), ddyd_data(i), g);
        Psi(:,i) = this.kernelFunction(x(i));
    end

    if (strcmpi(train_method,'LWR') == 1), this.w = LWR(Psi, s, Fd, this.zero_tol);
    elseif (strcmpi(train_method,'LS') == 1), this.w = leastSquares(Psi, s, Fd, this.zero_tol);
    else, error('[DMP_::train]: Unsopported training method...');
    end

    if (nargout > 0)
        F = zeros(size(Fd));
        for i=1:size(F,2)
            F(i) = this.calcLearnedFd(x(i), g);
        end
        train_error = norm(F-Fd)/length(F);
    end

%             figure;
%             hold on
%             for i=1:size(Psi,1)
%                 plot(x, Psi(i,:));
%             end
%             bar(this.c, this.w/(max(abs(this.w))));
%             axis tight;
%             hold off

%             figure
%             plot(Time, F, Time, Fd)
%             legend('F','Fd');
%             pause

end