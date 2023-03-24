clc;
close all;
clear;

%% ==============================

dt = 0.002;

n_dof = 2;
N_kernels = 30;

y0 = [0; 0];
g1 = [1; 1];
Tf = 4.0;

Timed = 0:dt:Tf;

%% ==============================
P1_data = get5thOrderPol(y0, g1, Timed);

gmp = GMP(n_dof, N_kernels);
t_start = tic;
train_mse = gmp.train('LS', Timed/Timed(end), P1_data);

%% ==============================

g2 = g1 + [-0.2; 0.3];
P2_data = get5thOrderPol(y0, g2, Timed);


%% ==============================

Time = [];
P_data = [];
dP_data = [];
ddP_data = [];
F_data = [];

t = 0;
y = y0;
y_dot = zeros(size(y));
g = g1;

can_sys = CanonicalSystem(Tf);

traj_scale = TrajScale_None(n_dof);
% traj_scale = TrajScale_Prop(n_dof);
dmp_pp = DMP_pp(gmp, traj_scale);
% dmp_pp.setRecursiveUpdate(false);
dmp_pp.setOptMetric('accel');

dmp_pp.r1 = 1e3*[1e-9 , 1e-4, 1e-3];
dmp_pp.rf = 1e6*[1e-1, 1e-0, 1e-0];

dmp_pp.init(can_sys.s, y0, g1, Tf);
dmp_pp.setAdaptToRobot(true);

K = 300;
D = 2*sqrt(K + 10);
f_gain = 400;
dt = 0.005;

Pd_data = [];

while (true)

    dmp_pp.update(g1, can_sys.s, can_sys.s_dot, y, y_dot);
    
    y_ref = dmp_pp.getRefPos(can_sys.s);
    y_ref_dot = dmp_pp.getRefVel(can_sys.s, can_sys.s_dot);
    y_ref_ddot = dmp_pp.getRefAccel(can_sys.s, can_sys.s_dot, can_sys.getPhaseDDot());
    
    [f, yd] = get_force(t, y, f_gain, Timed, P2_data);
    
    fmax = 0.5;
    a = 0;
    if (norm(f) > fmax)
        a = (min([norm(f), fmax]) / fmax)^0.5;
    end
    a = 1;
%     a = 0;
    
    y_ddot = (1-a)*y_ref_ddot + D*(y_ref_dot - y_dot) + K*(y_ref - y) + f;

    %g = g + 0.02*f;
    
    % log data
    Time = [Time t];
    P_data = [P_data y];
    dP_data = [dP_data y_dot];
    ddP_data = [ddP_data y_ddot];
    F_data = [F_data f];
    
    Pd_data = [Pd_data yd];
    
    % check termination
    if (t > Tf), break; end
    
    % numerical integration
    can_sys.integrate(t, t+dt);
    t = t + dt;
    y = y + y_dot*dt;
    y_dot = y_dot + y_ddot*dt;
    
end


figure; 
ax = axes();
hold on; grid on; box on;
plot(P1_data(1, :), P1_data(1, :), 'linewidth',2, 'color',0.4*[1 1 1], 'DisplayName','demo');
plot(P2_data(1, :), P2_data(2, :), 'linewidth',2, 'color','green', 'DisplayName','target');
plot(P_data(1, :), P_data(2, :), 'linewidth',2, 'color','blue', 'LineStyle',':', 'DisplayName','DMP++');
plot(y0(1), y0(2), 'Linestyle','none', 'Marker','o', 'linewidth',4, 'MarkerSize',10, 'color',[0 0.6 0], 'HandleVisibility','off');
plot(g1(1), g1(2), 'Linestyle','none', 'Marker','x', 'linewidth',2, 'MarkerSize',13, 'color',[1 0.7 0.7], 'HandleVisibility','off');
plot(g2(1), g2(2), 'Linestyle','none', 'Marker','x', 'linewidth',2, 'MarkerSize',13, 'color',[0.9 0 0], 'HandleVisibility','off');
% plot(Pd_data(1, :), Pd_data(2, :), 'linewidth',2, 'color',[0.85, 0.33, 0.1], 'LineStyle','--', 'DisplayName','yd');
legend({}, 'fontsize',15, 'orientation','horizontal', 'Position',[0.2172 0.9355 0.5875 0.0624], 'box','off');
xlabel('X [$m$]', 'interpreter','latex', 'fontsize',14);
ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',14);
axis tight;
ax.XLim = ax.XLim + 0.05*(ax.XLim(2) - ax.XLim(1)) * [-1 1];
ax.YLim = ax.YLim + 0.05*(ax.YLim(2) - ax.YLim(1)) * [-1 1];

F_norm = vecnorm(F_data, 2, 1);
figure;
plot(Time, F_norm, 'linewidth',2, 'Color','magenta');
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',14);
ylabel('[$N$]', 'interpreter','latex', 'fontsize',14);
title('$||F_{ext}||$', 'interpreter','latex', 'fontsize',14);
axis tight;

%% ====================================================

function [f, yd] = get_force(t, y, gain, Time, Pd_data)
    
    yd = interp1(Time, Pd_data', t)';
    
%     [~, k] = min(vecnorm(Pd_data - y, 2, 1));
%     yd = Pd_data(:,k);
    
    f = gain*(yd - y);

end


function [y, y_dot, y_ddot] = get5thOrderPol(y0, yf, Time)

    y0 = y0(:);
    yf = yf(:);
    Time = Time(:)';
    
    T = Time(end);
    t = Time/T;
    
    y = y0 + (yf - y0) * (10*t.^3 - 15*t.^4 + 6*t.^5 );

    if (nargout > 1)
        y_dot = (yf - y0) * (30*t.^2 - 60*t.^3 + 30*t.^4 ) / T;
    end

    if (nargout > 2)
        y_ddot = (yf - y0) * (60*t - 180*t.^2 + 120*t.^3 ) / T^2;
    end

end