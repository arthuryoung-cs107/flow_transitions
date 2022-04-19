clear
close all

run load_data.m
run set_figdims.m
run main_FIGSPECS.m

fig_num = 0;

% output information
write_figs = true;
write_all_figs = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';
%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------
fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axa, 'on');
      legend_set_1a(1) = errorbar(axa, PF1.omega, PF1.mu_torque, PF1.sigma_torque, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
      legend_set_1a(2) = errorbar(axa, PFR.omega, PFR.mu_torque, PFR.sigma_torque, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
    ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    legend(legend_set_1a,'Location', 'SouthEast', 'Interpreter', 'Latex');
    textbox_a = annotation('textbox', textbox_pos2_a_NW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axb, 'on');
      legend_set_1b(1) = errorbar(axb, NB1.omega_full, NB1.mu_torque_full, NB1.sigma_torque_full, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
      legend_set_1b(2) = errorbar(axb, NB2.omega, NB2.mu_torque, NB2.sigma_torque, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
      legend_set_1b(3) = errorbar(axb, NB3.omega, NB3.mu_torque, NB3.sigma_torque, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
    xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    legend(legend_set_1b,'Location', 'SouthEast', 'Interpreter', 'Latex');
    textbox_b = annotation('textbox', textbox_pos2_b_NW,   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

  axis([axa axb],omega_tau_range)
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axa, 'on');

       plot(axa, PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
       plot(axa, PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
       plot(axa, NB1.Re_s_noKD, NB1.G_rat_noKD, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
       plot(axa, NB2.Re_s_noKD, NB2.G_rat_noKD, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
       plot(axa, NB3.Re_s_noKD, NB3.G_rat_noKD, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)

    ylabel('$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$Re_s \textrm{ (PRE-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_a = annotation('textbox', textbox_pos2_a_SW,  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axb, 'on');

      legend_set_b(1) = plot(axb, PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
      legend_set_b(2) = plot(axb, PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
      legend_set_b(3) = plot(axb, NB1.Re_s, NB1.G_rat, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
      legend_set_b(4) = plot(axb, NB2.Re_s, NB2.G_rat, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
      legend_set_b(5) = plot(axb, NB3.Re_s, NB3.G_rat, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
      legend_set_b(6) = plot(axb, RM10.Re_s, RM10.G_rat, RM10.specs, 'Color', RM10.color, 'LineWidth', RM10.LW, 'DisplayName', RM10.label);
      legend_set_b(7) = plot(axb, RM20.Re_s, RM20.G_rat, RM20.specs, 'Color', RM20.color, 'LineWidth', RM20.LW, 'DisplayName', RM20.label);

    xlabel('$$Re_s \textrm{ (POST-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)
    legend(legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);
    textbox_b = annotation('textbox', textbox_pos2_b_SW,  'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

  axis([axa axb],Res_Grat_range)
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});

    axc = axes('Position', ax21(1, :));
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axc, 'on');
      fplot(@(Re) 1./(eta.*(Re)), [1e-1 5000], '-. k', 'Linewidth', 2)
      plot(axc, PF1.Re_s, PF1.cf, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
      plot(axc, PFR.Re_s, PFR.cf, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
      plot(axc, NB1.Re_s, NB1.cf, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
      plot(axc, NB2.Re_s, NB2.cf, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
      plot(axc, NB3.Re_s, NB3.cf, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
      plot(axc, RK.Re_s, RK.cf, RK.specs,'Color', RK.color, 'LineWidth', RK.LW)
    ylabel('$$c_f$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
    labels = {'$$ \frac{1}{\eta Re_s} $$', PF1.label, PFR.label, NB1.label, NB2.label, NB3.label, RK.label};
    legend(labels,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 3);
    textbox_c = annotation('textbox', textbox_pos21_c_SW,   'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
    textbox_c.FontSize = 16;

    axa = axes('Position', ax21(2, :));
    box on;
    set(gca, 'XScale', 'log')
    hold(axa, 'on');
      plot(axa, RV.Re_s, RV.alpha, RV.specs,'Color', RV.color, 'LineWidth', RV.LW, 'MarkerSize', RV.MS)
      plot(axa, LSa.Re_s, LSa.alpha, LSa.specs,'Color', LSa.color, 'LineWidth', LSa.LW)
      plot(axa, LSb.Re_s, LSb.alpha, LSb.specs,'Color', LSb.color, 'LineWidth', LSb.LW)
      plot(axa, PF1.Re_s_alpha, PF1.alpha, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
      plot(axa, PFR.Re_s_alpha, PFR.alpha, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
    ylabel('$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
    labels = {RV.label, LSa.label, LSb.label, PF1.label, PFR.label};
    legend(labels,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
    textbox_a = annotation('textbox', textbox_pos21_a_SW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;
    %
    axb = axes('Position', ax21(3, :));
    box on;
    set(gca, 'XScale', 'log')
    hold(axb, 'on');
      plot(axb, NB1.Re_s_alpha, NB1.alpha, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
      plot(axb, NB2.Re_s_alpha, NB2.alpha, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
      plot(axb, NB3.Re_s_alpha, NB3.alpha, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
    xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
    labels = {NB1.label, NB2.label, NB3.label};
    legend(labels,'Location', 'SouthEast', 'Interpreter', 'Latex');
    textbox_b = annotation('textbox', textbox_pos21_b_SW,   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

    axis(axc,Res_cf_range)
    axis([axa axb], Res_alpha_range)


%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axa, 'on');
      errorbar(axa, UB1.omega, UB1.mu_torque, UB1.sigma_torque, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS)
      errorbar(axa, UB2.omega, UB2.mu_torque, UB2.sigma_torque, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS)
    ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    labels = {UB1.label, UB2.label};
    legend(labels,'Location', 'SouthEast', 'Interpreter', 'Latex');
    textbox_a = annotation('textbox', textbox_pos2_a_NW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axb, 'on');
      errorbar(axb, XB1.omega, XB1.mu_torque, XB1.sigma_torque, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS)
      errorbar(axb, XB2.omega, XB2.mu_torque, XB2.sigma_torque, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS)
    xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    labels = {XB1.label, XB2.label};
    legend(labels,'Location', 'SouthEast', 'Interpreter', 'Latex');
    textbox_b = annotation('textbox', textbox_pos2_b_NW,   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

  axis([axa axb],omega_tau_range)
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  axc = axes('Position', ax21(1, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(axc, 'on');

    fplot(@(Re) 1./(eta.*(Re)), [1e-1 5000], '-. k', 'Linewidth', 2)
    plot(axc, UB1.Re_s, UB1.cf, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS)
    plot(axc, UB2.Re_s, UB2.cf, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS)
    plot(axc, XB1.Re_s, XB1.cf, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS)
    plot(axc, XB2.Re_s, XB2.cf, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS)

  ylabel('$$c_f$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels = {'$$ \frac{1}{\eta Re_s} $$', UB1.label, UB2.label, XB1.label, XB2.label};
  legend(labels,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
  textbox_c = annotation('textbox', textbox_pos21_c_SW, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
  textbox_c.FontSize = 16;

  axa = axes('Position', ax21(2, :));
  box on;
  set(gca, 'XScale', 'log')
  hold(axa, 'on');

    plot(axa, UB1.Re_s_alpha, UB1.alpha, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS)
    plot(axa, UB2.Re_s_alpha, UB2.alpha, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS)

  ylabel('$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels = {UB1.label, UB2.label};
  legend(labels,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox_a = annotation('textbox', textbox_pos21_a_NW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox_a.FontSize = 16;
  %
  axb = axes('Position', ax21(3, :));
  box on;
  set(gca, 'XScale', 'log')
  hold(axb, 'on');

    plot(axb, XB1.Re_s_alpha, XB1.alpha, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS)
    plot(axb, XB2.Re_s_alpha, XB2.alpha, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS)

  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels = {XB1.label, XB2.label};
  legend(labels,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox_b = annotation('textbox', textbox_pos21_b_NW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox_b.FontSize = 16;

  axis(axc,Res_cf_range)
  axis([axa axb], Res_alpha_range)

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   6  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(3, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axa, 'on');

      fplot(axa, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
      legend_set_6a(1) = plot(axa, PF1.Re_s, PF1.G, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
      legend_set_6a(2) = plot(axa, PFR.Re_s, PFR.G, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
      fplot(axa, @(Re) (PF1.powerfit.b).*(Re).^(PF1.powerfit.m), [71 10000],'-', 'Color', blue5,'Linewidth', 2, 'DisplayName', 'PF1 $$\beta Re_s^{\alpha}$$')
      fplot(axa, @(Re) (PFR.powerfit.b).*(Re).^(PFR.powerfit.m), [71 10000],'-', 'Color', PFR.color,'Linewidth', 2, 'DisplayName', 'PF2 $$\beta Re_s^{\alpha}$$')

    ylabel('$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    legend(legend_set_6a,'Location', 'NorthWest', 'Interpreter', 'Latex');
    textbox_a = annotation('textbox', textbox_pos222_a, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axb, 'on');

      fplot(axb, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
      legend_set_6b(1) = plot(axb, NB1.Re_s, NB1.G, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
      legend_set_6b(2) = plot(axb, NB2.Re_s, NB2.G, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
      legend_set_6b(3) = plot(axb, NB3.Re_s, NB3.G, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
      fplot(axb, @(Re) (NBall.powerfit.b).*(Re).^(NBall.powerfit.m), [71 10000],'-', 'Color', NBall.color,'Linewidth', 2, 'DisplayName', 'NB $$\beta Re_s^{\alpha}$$')

    legend(legend_set_6b,'Location', 'NorthWest', 'Interpreter', 'Latex');
    textbox_b = annotation('textbox', textbox_pos222_b, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

    axc = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axc, 'on');

      fplot(axc, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
      legend_set_6c(1) = plot(axc, UB1.Re_s, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
      legend_set_6c(2) = plot(axc, UB2.Re_s, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
      fplot(axc, @(Re) (UBall.powerfit.b).*(Re).^(UBall.powerfit.m), [71 10000],'-', 'Color', UBall.color,'Linewidth', 2, 'DisplayName', 'UB $$\beta Re_s^{\alpha}$$')

    ylabel('$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    legend(legend_set_6c,'Location', 'NorthWest', 'Interpreter', 'Latex');
    textbox_c = annotation('textbox', textbox_pos222_c, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
    textbox_c.FontSize = 16;

    axd = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axd, 'on');

      fplot(axd, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
      legend_set_6d(1) = plot(axd, XB1.Re_s, XB1.G, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
      legend_set_6d(2) = plot(axd, XB2.Re_s, XB2.G, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);
      fplot(axd, @(Re) (XBall.powerfit.b).*(Re).^(XBall.powerfit.m), [71 10000],'-', 'Color', XBall.color,'Linewidth', 2, 'DisplayName', 'XB $$\beta Re_s^{\alpha}$$')

    legend(legend_set_6d,'Location', 'NorthWest', 'Interpreter', 'Latex');
    textbox_d = annotation('textbox', textbox_pos222_d, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none');
    textbox_d.FontSize = 16;

    axe = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axe, 'on');

      fplot(axe, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
      for i=1:length(FB1.exp)
        legend_set_6e(i) = plot(axe, FB1.exp(i).Re_s, FB1.exp(i).G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
      end
      fplot(axe, @(Re) (FB1.powerfit.b).*(Re).^(FB1.powerfit.m), [71 10000],'-', 'Color', FB1.color,'Linewidth', 2, 'DisplayName', 'FB1 $$\beta Re_s^{\alpha}$$');

    ylabel('$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
    legend(legend_set_6e,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);
    textbox_e = annotation('textbox', textbox_pos222_e, 'Interpreter', 'LaTeX', 'String', 'e)', 'LineStyle', 'none');
    textbox_e.FontSize = 16;

    axf = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axf, 'on');

      fplot(axf, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
      for i=1:length(FB2.exp)
        legend_set_6f(i) = plot(axf, FB2.exp(i).Re_s, FB2.exp(i).G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
      end
      fplot(axf, @(Re) (FB2.powerfit.b).*(Re).^(FB2.powerfit.m), [71 10000],'-', 'Color', FB2.color,'Linewidth', 2, 'DisplayName', 'FB2 $$\beta Re_s^{\alpha}$$');
    legend(legend_set_6f,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);
    xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_f = annotation('textbox', textbox_pos222_f, 'Interpreter', 'LaTeX', 'String', 'f)', 'LineStyle', 'none');
    textbox_f.FontSize = 16;

  axis([axa axb axc axd axe axf],Res_G_range)
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   7  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axa, 'on');
      errorbar(axa, EC000.omega, EC000.mu_torque, EC000.sigma_torque, EC000.specs, 'Color', EC000.color, 'LineWidth', EC000.LW, 'MarkerSize', EC000.MS, 'DisplayName', EC000.label)
      errorbar(axa, EC050.omega, EC050.mu_torque, EC050.sigma_torque, EC050.specs, 'Color', EC050.color, 'LineWidth', EC050.LW, 'MarkerSize', EC050.MS, 'DisplayName', EC050.label)
      errorbar(axa, EC075.omega, EC075.mu_torque, EC075.sigma_torque, EC075.specs, 'Color', EC075.color, 'LineWidth', EC075.LW, 'MarkerSize', EC075.MS, 'DisplayName', EC075.label)
      errorbar(axa, EC100.omega, EC100.mu_torque, EC100.sigma_torque, EC100.specs, 'Color', EC100.color, 'LineWidth', EC100.LW, 'MarkerSize', EC100.MS, 'DisplayName', EC100.label)
      for i = 1:length(FB1.exp)
        errorbar(axa, FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.exp(i).sigma_torque, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      end
    ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    legend('Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)
    textbox_a = annotation('textbox', textbox_pos2_a_SW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axb, 'on');
    for i = 1:length(FB2.exp)
        errorbar(axb, FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.exp(i).sigma_torque, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
    end
    xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    legend('Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)
    textbox_b = annotation('textbox', textbox_pos2_b_SW,   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

    axis([axa axb],omega_tau_range)
    tile_object.TileSpacing = 'compact';
    tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   8  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(1, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axa, 'on');
      for i = 1:length(FB1.exp)
        plot(axa, FB1.exp(i).omega, FB1.exp(i).appmu, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      end
    ylabel('$$\mu_{app}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$\omega_{i}$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    legend('Show', 'Location', 'NorthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)
    textbox_a = annotation('textbox', textbox_pos2_a_NW,  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold(axb, 'on');
      for i = 1:length(FB2.exp)
          plot(axb, FB2.exp(i).omega, FB2.exp(i).appmu, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
      end
    xlabel('$$\omega_{i}$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
    legend('Show', 'Location', 'NorthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)
    textbox_b = annotation('textbox', textbox_pos2_b_NW,  'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

  axis([axa axb],omega_appmu_range)
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   9  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(2, 2);

    axa = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    hold(axa, 'on');
      for i = 1:length(FB1.exp)
          plot(axa, FB1.exp(i).q, (FB1.exp(i).tau_y)/(FB1.exp(i).tau_static), FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label)
      end
    ylabel('$$\tau_y/\tau_{q = 0}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    % xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_a = annotation('textbox', textbox_pos22_a_NE,  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
    textbox_a.FontSize = 16;

    axb = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    hold(axb, 'on');
      for i = 1:length(FB2.exp)
          plot(axb, FB2.exp(i).q, (FB2.exp(i).tau_y)/(FB2.exp(i).tau_static), FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label)
      end
    % xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_b = annotation('textbox', textbox_pos22_b_NE,  'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
    textbox_b.FontSize = 16;

    axc = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    hold(axc, 'on');
      for i = 1:length(FB1.exp)
          plot(axc, FB1.exp(i).q, FB1.exp(i).mu_p, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label)
      end
    ylabel('$$\tilde{\mu}_{p}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
    xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_c = annotation('textbox', textbox_pos22_c_NE,  'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
    textbox_c.FontSize = 16;

    axd = nexttile;
    box on;
    set(gca, 'YScale', 'log')
    hold(axd, 'on');
      for i = 1:length(FB2.exp)
          plot(axd, FB2.exp(i).q, FB2.exp(i).mu_p, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label)
      end
    xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    textbox_d = annotation('textbox', textbox_pos22_d_NE, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none');
    textbox_d.FontSize = 16;

  axis(axa,[0 2 1e-5 1.0])
  axis(axb,[0 16 1e-5 1.0])
  axis(axc,[0 2 5e-2 1.2])
  axis(axd,[0 16 5e-2 1.2])
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

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
%%%%%%%% -----------------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

if (write_figs)
  if (write_all_figs)
    figs_to_write = 1:length(figs);
  end
  AYfig.save_figs(figs, figs_to_write, save_type, save_dir);
end
