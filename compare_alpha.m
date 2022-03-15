clear
close all

run load_data.m
run set_figdims.m
run compare_alpha_FIGSPECS.m

fig_num = 9;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  10  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(2, 2);

    axa = nexttile;
    box on;
    set(gca, 'XScale', 'log')
    hold(axa, 'on');
      for i = 1:length(FB1.exp)
        plot(axa, FB1.exp(i).omega, FB1.exp(i).alpha, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      end
    ylabel('$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',12)
    % legend('Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)
    textbox_a = annotation('textbox', textbox_pos2_a_NW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'XScale', 'log')
    hold(axb, 'on');
    for i = 1:length(FB2.exp)
        plot(axb, FB2.exp(i).omega, FB2.exp(i).alpha, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
    end
    xlabel('$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',12)
    % legend('Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)
    textbox_b = annotation('textbox', textbox_pos2_b_NW,   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

    axc = nexttile;
    box on;
    hold(axc, 'on');
      for i = 1:length(FB1.exp)
          plot(axc, FB1.exp(i).q, FB1.exp(i).powerfit.m, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label)
      end
      ylabel('$$\alpha$$ fitted [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_c = annotation('textbox', textbox_pos22_c_NE,  'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
    textbox_c.FontSize = 16;

    axd = nexttile;
    box on;
    hold(axd, 'on');
      for i = 1:length(FB2.exp)
          plot(axd, FB2.exp(i).q, FB2.exp(i).powerfit.m, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label)
      end
    xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_d = annotation('textbox', textbox_pos22_d_NE, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none');
    textbox_d.FontSize = 16;

    axis([axa axb], [1e-2, 2e2, -1, 3])
    axis(axc,[0 2 -1 3])
    axis(axd,[0 16 -1 3])
    tile_object.TileSpacing = 'compact';
    tile_object.Padding = 'compact';
%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  11  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(2, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axa, 'on');
      for i = 1:length(FB1.exp)
          fplot(axa, @(Re) (FB1.exp(i).powerfit.b).*(Re).^(FB1.exp(i).powerfit.m), [71 10000],'-', 'Color', FB1.exp(i).color,'Linewidth', 2, 'DisplayName', FB1.exp(i).label)
      end
    ylabel('$$G$$', 'Interpreter', 'LaTeX','FontSize',12)
    % xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_a = annotation('textbox', textbox_pos22_a_NE,  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axb, 'on');
      for i = 1:length(FB2.exp)
        fplot(axb, @(Re) (FB2.exp(i).powerfit.b).*(Re).^(FB2.exp(i).powerfit.m), [71 10000],'-', 'Color', FB2.exp(i).color,'Linewidth', 2, 'DisplayName', FB2.exp(i).label)
      end
    % xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_b = annotation('textbox', textbox_pos22_b_NE,  'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

    axc = nexttile;
    box on;
    hold(axc, 'on');
      for i = 1:length(FB1.exp)
          plot(axc, FB1.exp(i).q, FB1.exp(i).powerfit.m, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label)
      end
    ylabel('$$\alpha$$ ', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_c = annotation('textbox', textbox_pos22_c_NE,  'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
    textbox_c.FontSize = 16;

    axd = nexttile;
    box on;
    hold(axd, 'on');
      for i = 1:length(FB2.exp)
          plot(axd, FB2.exp(i).q, FB2.exp(i).powerfit.m, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label)
      end
    xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_d = annotation('textbox', textbox_pos22_d_NE, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none');
    textbox_d.FontSize = 16;

  axis([axa, axb],Res_G_range)
  axis(axc,[0 2 -1 3])
  axis(axd,[0 16 -1 3])
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  12  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 1);
  alpha_high = 3.0;
  alpha_low = -0.5;
    axa = nexttile;
    view(view_mat(1, :));
    box on;

    set(gca, 'XScale', 'log')

    hold(axa, 'on');

      for i=1:FB1.len
        legend_set_1a(i) = scatter3(axa, FB1.exp(i).omega, FB1.exp(i).q*ones(size(FB1.exp(i).omega)), FB1.exp(i).alpha, 's','CData', FB1.exp(i).alpha, 'LineWidth', FB1.LW_L, 'DisplayName', FB1.exp(i).label);
      end
      colormap('Cool')
    ylabel('$$q$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    zlabel('$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)
    zlim([-0.5 3])

  % axis(axa,omega_tau_range)
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  13  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 1);

  axa = nexttile;
  view(view_mat(1, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(axa, 'on');

    for i=1:FB2.len
      legend_set_1a(i) = scatter3(axa, FB2.exp(i).omega, FB2.exp(i).q*ones(size(FB2.exp(i).omega)), FB2.exp(i).alpha, 'o','CData', FB2.exp(i).q*ones(size(FB2.exp(i).omega)), 'LineWidth', FB1.LW_L, 'DisplayName', FB2.exp(i).label);
    end
    colormap('Cool')
  ylabel('$$q$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
  zlabel('$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)
  zlim([-1 3])

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  14  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 1);

  axa = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axa, 'on');

  yyaxis(axa,'left')
  set(gca, 'YScale', 'log')

    plot(axa, PF1.omega, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
    ylabel('$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

  yyaxis(axa,'right')

    plot(axa, PF1.omega, PF1.alpha, PF1.specs,'Color', PFR.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
    % ylim([-0.5, 2]);
    ylabel('$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)

  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  15  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 1);

  axa = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axa, 'on');

  yyaxis(axa,'left')
  set(gca, 'YScale', 'log')

    plot(axa, PFR.omega, PFR.G_rat, PFR.specs,'Color', PF1.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
    ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)

  yyaxis(axa,'right')

    plot(axa, PFR.omega, PFR.alpha, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
    % ylim([-0.5, 2]);
    ylabel('$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)

  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  16  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 1);

  axa = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axa, 'on');

  yyaxis(axa,'left')
  set(gca, 'YScale', 'log')
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  yyaxis(axa,'right')
  ylabel('$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)
  ylim([-0.5, 2]);
    for i=1:FB1.len
      yyaxis(axa,'left')
      plot(axa, FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.exp(i).specs,'Color', blue14(i, :), 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
      yyaxis(axa,'right')
      plot(axa, FB1.exp(i).omega, FB1.exp(i).alpha, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
    end

  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  17  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 1);

  axa = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axa, 'on');

  yyaxis(axa,'left')
  set(gca, 'YScale', 'log')
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  yyaxis(axa,'right')
  ylabel('$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)
  ylim([-0.5, 2]);
    for i=1:FB2.len
      yyaxis(axa,'left')
      plot(axa, FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.exp(i).specs,'Color', blue14(i, :), 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
      yyaxis(axa,'right')
      plot(axa, FB2.exp(i).omega, FB2.exp(i).alpha, FB2.exp(i).specs,'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
    end



  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)


%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------
