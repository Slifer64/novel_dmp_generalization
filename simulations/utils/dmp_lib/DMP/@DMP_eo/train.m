function train_err = train(this, train_method, Time, Quat_data, rotVel_data, rotAccel_data)

    tau = Time(end);
    this.setTau(tau);

    n_data = length(Time);
    eo_data = zeros(3, n_data);
    deo_data = zeros(3, n_data);
    ddeo_data = zeros(3, n_data);
    
    Qg = Quat_data(:,end);

    for j=1:n_data
       Qe = quatProd(Qg, quatInv(Quat_data(:,j)) );
       eo_data(:,j) = quatLog(Qe);
       deo_data(:,j) = DMP_eo.rotVel2deo(rotVel_data(:,j), Qe);
       ddeo_data(:,j) = DMP_eo.rotAccel2ddeo(rotAccel_data(:,j), rotVel_data(:,j), Qe);
    end

%             Ts = Time(end) - Time(end-1);
%             dTime = [diff(Time) Ts];
%             for i=1:3
%                 % deo_data(i,:) = [diff(eo_data(i,:)) 0] ./ dTime;
%                 % ddeo_data(i,:) = [diff(deo_data(i,:)) 0] ./ dTime;
%             end
%             
%             deo_data2 = zeros(3, n_data);
%             ddeo_data2 = zeros(3, n_data);
%             rotVel_data2 = zeros(3, n_data);
%             rotAccel_data2 = zeros(3, n_data);
%             
%             for j=floor(n_data/2):n_data
%                
%                Qe = quatProd(Qg, quatInv(Quat_data(:,j)) );
%                
%                rotVel = rotVel_data(:,j);
%                rotAccel = rotAccel_data(:,j);
%                
%                J_dq_deo = this.jacobDquatDeo(Qe);
%                dJ_dq_deo = this.jacobDotDquatDeo(Qe, rotVel);
%                J_deo_dq = this.jacobDeoDquat(Qe);
%                dJ_deo_dq = this.jacobDotDeoDquat(Qe, rotVel);
%                P = J_dq_deo * J_deo_dq;
%                P_dot = dJ_dq_deo*J_deo_dq + J_dq_deo*dJ_deo_dq;
% 
%                rotVelQ = [0; rotVel];
%                rotAccelQ = [0; rotAccel];
%                
%                % norm( P*quatProd(Qe,rotVelQ) - quatProd(Qe,rotVelQ) )
%                %norm( quatProd(P*quatProd(Qe,rotVelQ),rotVelQ) - quatProd(quatProd(Qe,rotVelQ),rotVelQ) )
%                
%                % norm( (P*quat2mat(Qe))*rotVelQ - quatProd(Qe,rotVelQ) )
%                % norm( quatProd(P*Qe,rotVelQ) - quatProd(Qe,rotVelQ) )
%                
% %                norm( P_dot*quatProd(Qe,rotVelQ) - 0.5*P*quatProd(quatProd(Qe,rotVelQ),rotVelQ) + P*quatProd(Qe,rotAccelQ) ...
% %                      + 0.5*quatProd(quatProd(Qe,rotVelQ),rotVelQ) - quatProd(Qe, rotAccelQ) )
%                
%                temp1 = P*quatProd(Qe,rotAccelQ);
%                temp2 = quatProd(Qe, rotAccelQ) + 0.5*(P-eye(4,4))*quatProd(quatProd(Qe,rotVelQ),rotVelQ) - P_dot*quatProd(Qe,rotVelQ);
%                norm(temp1 - temp2)
%                  
%                pause
% 
%                % deo_data2(:,j) = this.rotVel2deo(rotVel_data(:,j), Qe);
%                % ddeo_data2(:,j) = this.rotAccel2ddeo(rotAccel_data(:,j), rotVel_data(:,j), Qe);
%                
%                % rotVel_data2(:,j) = this.deo2rotVel(deo_data2(:,j), Qe);
%                % rotAccel_data2(:,j) = this.ddeo2rotAccel(ddeo_data(:,j), rotVel_data(:,j), Qe);
%                
%             end
% 
%             figure;
%             subplot(2,1,1);
%             hold on;
%             plot(Time, rotVel_data(1,:), 'LineWidth',1.5, 'Color','blue');
%             plot(Time, rotVel_data2(1,:), 'LineWidth',1.5, 'Color','magenta');
%             legend({'$\omega$','$\omega_2$'}, 'interpreter','latex', 'fontsize',14);
%             axis tight;
%             hold off;
%             subplot(2,1,2);
%             plot(Time, abs(rotVel_data(1,:)-rotVel_data2(1,:)), 'LineWidth',1.5, 'Color','red');
%             legend('error');
%             axis tight;
%             
%             figure;
%             subplot(2,1,1);
%             hold on;
%             plot(Time, deo_data(1,:), 'LineWidth',1.5, 'Color','blue');
%             plot(Time, deo_data2(1,:), 'LineWidth',1.5, 'Color','magenta');
%             legend({'$\dot{e}_o$','$\dot{e}_{o2}$'}, 'interpreter','latex', 'fontsize',14);
%             axis tight;
%             hold off;
%             subplot(2,1,2);
%             plot(Time, abs(deo_data(1,:)-deo_data2(1,:)), 'LineWidth',1.5, 'Color','red');
%             legend('error');
%             axis tight;
%             
%             figure;
%             subplot(2,1,1);
%             hold on;
%             plot(Time, rotAccel_data(1,:), 'LineWidth',1.5, 'Color','blue');
%             plot(Time, rotAccel_data2(1,:), 'LineWidth',1.5, 'Color','magenta');
%             legend({'$\dot{\omega}$','$\dot{\omega}_2$'}, 'interpreter','latex', 'fontsize',14);
%             axis tight;
%             hold off;
%             subplot(2,1,2);
%             plot(Time, abs(rotAccel_data(1,:)-rotAccel_data2(1,:)), 'LineWidth',1.5, 'Color','red');
%             legend('error');
%             axis tight;
%             
%             figure;
%             subplot(2,1,1);
%             hold on;
%             plot(Time, ddeo_data(1,:), 'LineWidth',1.5, 'Color','blue');
%             plot(Time, ddeo_data2(1,:), 'LineWidth',1.5, 'Color','magenta');
%             legend({'$\ddot{e}_o$','$\ddot{e}_{o2}$'}, 'interpreter','latex', 'fontsize',14);
%             axis tight;
%             hold off;
%             subplot(2,1,2);
%             plot(Time, abs(ddeo_data(1,:)-ddeo_data2(1,:)), 'LineWidth',1.5, 'Color','red');
%             legend('error');
%             axis tight;
%             
%             error('stop')

    if (nargout > 0)
        train_err = zeros(3,1);
        for i=1:3
            train_err(i) = this.dmp{i}.train(train_method, Time, eo_data(i,:), deo_data(i,:), ddeo_data(i,:));
        end
    else
        for i=1:3
            this.dmp{i}.train(train_method, Time, eo_data(i,:), deo_data(i,:), ddeo_data(i,:));
        end
    end

end
        