%% Plot results

ax_font = 13;
x_font = 16;
y_font = 16;
legend_font = 17;


for i=1:n_dof
    fig = figure;
    ax_vec = [];
    fig.Position(3:4) = [581 656];
    % Plot Position
    ax = subplot(3,1,1); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    if exist('Pg_data','var') == 1
        plot(Time, Pg_data(i,:), 'LineWidth',1.5, 'LineStyle','-', 'Color',[1, 0.0, 0.0], 'HandleVisibility','off');
    end
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue', 'DisplayName', 'DMP$^{++}$');
    dmp_pl = plot(Time2, P2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta', 'DisplayName','DMP');
    plot(Timed, Pd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green', 'DisplayName', 'demo');
    plot(Time(end), gd(i), 'LineWidth',3, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color',[1 0.6 0.6], 'DisplayName','$g_d$');
    plot(Time(end), g(i), 'LineWidth',3, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color','red', 'DisplayName','$g$');
    ax.FontSize = ax_font;
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',y_font);
    % title(['DoF $' num2str(i) '$'], 'interpreter','latex', 'fontsize',18);
    legend({}, 'interpreter','latex', 'fontsize',legend_font, 'Position',[0.2251 0.9396 0.5708 0.0425], 'Orientation','horizontal', 'Box','off');
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    % Plot Velocity
    ax = subplot(3,1,2); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time2, dP2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    plot(Timed, dPd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green');
    ax.FontSize = ax_font;
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',y_font);
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    % Plot Acceleration
    ax = subplot(3,1,3); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time2, ddP2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    plot(Timed, ddPd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green');
    ax.FontSize = ax_font;
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',y_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',x_font);
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    linkaxes(ax_vec, 'x');
    for i=1:length(ax_vec)
       ax = ax_vec(i);
       ax.Box = 'on';
       ax.YLim = ax.YLim + 0.07*(ax.YLim(2)-ax.YLim(1))*[-1 1];
    end
    ax = ax_vec(i);
    ax.XLim(2) = ax.XLim(2) + 0.05;
end