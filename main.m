clear
close all

run main_RUNSCRIPT.m

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

figs(1) = figure('Name', 'PF_NB_T_vs_omega' ,'Renderer', 'painters', 'Position', [pos11, dim1]);
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  ylim([1e-8, 1e-3]);
  xlim([0.0628  100]);
  hold on
  box on;
    errorbar(PF1.omega, PF1.mu_torque, PF1.sigma_torque, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
    errorbar(PFR.omega, PFR.mu_torque, PFR.sigma_torque, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
    errorbar(NB1.omega_full, NB1.mu_torque_full, NB1.sigma_torque_full, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
    errorbar(NB2.omega, NB2.mu_torque, NB2.sigma_torque, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
    errorbar(NB3.omega, NB3.mu_torque, NB3.sigma_torque, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_{i}$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
  labels1 = {PF1.label, PFR.label, NB1.label, NB2.label, NB3.label};
  legend(labels1,'Location', 'NorthWest', 'Interpreter', 'Latex');
  hold off

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

figs(2) = figure('Name', 'PF_NB_Grat_vs_Res' ,'Renderer', 'painters', 'Position', [pos21 dim2]);
  tile_object = tiledlayout(1, 2);
  ax1 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax1, 'on');
    plot(ax1, PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
    plot(ax1, PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
    plot(ax1, NB1.Re_s_noKD, NB1.G_rat_noKD, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
    plot(ax1, NB2.Re_s_noKD, NB2.G_rat_noKD, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
    plot(ax1, NB3.Re_s_noKD, NB3.G_rat_noKD, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
  ylabel('$$G_{rat} = \frac{G}{G_{cc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s \textrm{ (PRE-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)
  textbox2_a = annotation('textbox', [0.15, 0.825, 0.1, 0.1],  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox2_a.FontSize = 16;

  ax2 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax2, 'on');
    plot(ax2, PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
    plot(ax2, PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
    plot(ax2, NB1.Re_s, NB1.G_rat, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
    plot(ax2, NB2.Re_s, NB2.G_rat, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
    plot(ax2, NB3.Re_s, NB3.G_rat, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
    plot(ax2, RM10.Re_s, RM10.G_rat, RM10.specs, 'Color', RM10.color, 'LineWidth', 4)
    plot(ax2, RM20.Re_s, RM20.G_rat, RM20.specs, 'Color', RM20.color, 'LineWidth', 4)
  xlabel('$$Re_s \textrm{ (POST-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels2 = {PF1.label, PFR.label, NB1.label, NB2.label, NB3.label, RM10.label, RM20.label};
  legend(labels2,'Location', 'NorthWest', 'Interpreter', 'Latex');
  textbox2_b = annotation('textbox', [0.85, 0.825, 0.1, 0.1],  'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox2_b.FontSize = 16;
  axis([ax1 ax2],[1e-1 2e4 0.8 1e2])


%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

figs(3) = figure('Name', 'PF_NB_cf_vs_Res' ,'Renderer', 'painters', 'Position', [pos13 dim21]);
  ax1 = axes('Position', ax21(1, :));
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax1, 'on');
    fplot(@(Re) 1./(eta.*(Re)), [1e-1 5000], '-. k', 'Linewidth', 2)
    plot(ax1, PF1.Re_s, PF1.cf, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
    plot(ax1, PFR.Re_s, PFR.cf, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
    plot(ax1, NB1.Re_s, NB1.cf, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
    plot(ax1, NB2.Re_s, NB2.cf, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
    plot(ax1, NB3.Re_s, NB3.cf, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
    plot(ax1, RK.Re_s, RK.cf, RK.specs,'Color', RK.color, 'LineWidth', 2)
  ylabel('$$c_f$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels3_c = {'$$ \frac{1}{\eta Re_s} $$', PF1.label, PFR.label, NB1.label, NB2.label, NB3.label, RK.label};
  legend(labels3_c,'Location', 'NorthEast', 'Interpreter', 'Latex');
  textbox3_c = annotation('textbox', [0.12, 0.07, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
  textbox3_c.FontSize = 16;

  ax2 = axes('Position', ax21(2, :));
  box on;
  set(gca, 'XScale', 'log')
  hold(ax2, 'on');
    plot(ax2, RV.Re_s, RV.alpha, RV.specs,'Color', RV.color, 'LineWidth', 2, 'MarkerSize', 10)
    plot(ax2, LSa.Re_s, LSa.alpha, LSa.specs,'Color', LSa.color, 'LineWidth', 2)
    plot(ax2, LSb.Re_s, LSb.alpha, LSb.specs,'Color', LSb.color, 'LineWidth', 2)
    plot(ax2, PF1.Re_s_alpha, PF1.alpha, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
    plot(ax2, PFR.Re_s_alpha, PFR.alpha, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
  ylabel('$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels3_a = {RV.label, LSa.label, LSb.label};
  legend(labels3_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox3_a = annotation('textbox', [0.12, 0.875, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox3_a.FontSize = 16;
  hold(ax2, 'off');
  %
  ax3 = axes('Position', ax21(3, :));
  box on;
  set(gca, 'XScale', 'log')
  hold(ax3, 'on');
    plot(ax3, NB1.Re_s_alpha, NB1.alpha, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
    plot(ax3, NB2.Re_s_alpha, NB2.alpha, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
    plot(ax3, NB3.Re_s_alpha, NB3.alpha, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels3_b = {NB1.label, NB2.label, NB3.label};
  legend(labels3_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox3_b = annotation('textbox', [0.55, 0.875, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox3_b.FontSize = 16;
  axis(ax1,[1e-1 1e4 1e-4 1e4])
  axis([ax2 ax3],[1e0 5e5 0.9 2])
  hold(ax3, 'off');


%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


figs(4) = figure('Name', 'UB_XB_T_vs_omega' ,'Renderer', 'painters', 'Position', [pos23 dim2]);
  tile_object = tiledlayout(1, 2);
  ax1 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax1, 'on');
    errorbar(ax1, UB1.omega, UB1.mu_torque, UB1.sigma_torque, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS)
    errorbar(ax1, UB2.omega, UB2.mu_torque, UB2.sigma_torque, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS)
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
  labels4_a = {UB1.label, UB2.label};
  legend(labels4_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox4_a = annotation('textbox', [0.10, 0.85, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox4_a.FontSize = 16;
  hold(ax1, 'off');

  ax2 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax2, 'on');
    errorbar(ax2, XB1.omega, XB1.mu_torque, XB1.sigma_torque, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS)
    errorbar(ax2, XB2.omega, XB2.mu_torque, XB2.sigma_torque, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
  labels4_b = {XB1.label, XB2.label};
  legend(labels4_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox4_b = annotation('textbox', [0.575, 0.85, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox4_b.FontSize = 16;
  axis([ax1 ax2],[1e-2 2e2 1e-6 1e-2])
  hold(ax2, 'off');

  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

figs(5) = figure('Name', 'UB_XB_cf_vs_Res' ,'Renderer', 'painters', 'Position', [pos12 dim21]);
  ax1 = axes('Position', ax21(1, :));
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax1, 'on');
    fplot(@(Re) 1./(eta.*(Re)), [1e-1 5000], '-. k', 'Linewidth', 2)
    plot(ax1, UB1.Re_s, UB1.cf, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS)
    plot(ax1, UB2.Re_s, UB2.cf, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS)
    plot(ax1, XB1.Re_s, XB1.cf, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS)
    plot(ax1, XB2.Re_s, XB2.cf, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS)
  ylabel('$$c_f$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels5_c = {'$$ \frac{1}{\eta Re_s} $$', UB1.label, UB2.label, XB1.label, XB2.label};
  legend(labels5_c,'Location', 'NorthEast', 'Interpreter', 'Latex');
  textbox5_c = annotation('textbox', [0.12, 0.07, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
  textbox5_c.FontSize = 16;

  ax2 = axes('Position', ax21(2, :));
  box on;
  set(gca, 'XScale', 'log')
  hold(ax2, 'on');
    plot(ax2, UB1.Re_s_alpha, UB1.alpha, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS)
    plot(ax2, UB2.Re_s_alpha, UB2.alpha, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS)
  ylabel('$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels5_a = {UB1.label, UB2.label};
  legend(labels5_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox5_a = annotation('textbox', [0.12, 0.875, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox5_a.FontSize = 16;
  hold(ax2, 'off');
  %
  ax3 = axes('Position', ax21(3, :));
  box on;
  set(gca, 'XScale', 'log')
  hold(ax3, 'on');
    plot(ax3, XB1.Re_s_alpha, XB1.alpha, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS)
    plot(ax3, XB2.Re_s_alpha, XB2.alpha, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  labels5_b = {XB1.label, XB2.label};
  legend(labels5_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
  textbox5_b = annotation('textbox', [0.55, 0.875, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox5_b.FontSize = 16;

  axis(ax1,[1e-1 1e4 1e-4 1e4])
  axis([ax2 ax3],[1e0 5e5 -1 2])
  hold(ax3, 'off');

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   6  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

figs(6) = figure('Name', 'PF_NB_UB_XB_FB_G_vs_Res' ,'Renderer', 'painters', 'Position', [pos23 dim222]);
  ax1 = axes('Position', ax222(1, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax1, 'on');
    fplot(ax1, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2)
    plot(ax1, PF1.Re_s, PF1.G, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
    plot(ax1, PFR.Re_s, PFR.G, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
    fplot(ax1, @(Re) (PF1.powerfit.b).*(Re).^(PF1.powerfit.m), [71 10000],'-', 'Color', blue5,'Linewidth', 2)
    fplot(ax1, @(Re) (PFR.powerfit.b).*(Re).^(PFR.powerfit.m), [71 10000],'-', 'Color', PFR.color,'Linewidth', 2)
  ylabel('$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  labels6_a = {'$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$', PF1.label, PFR.label, 'PF1 $$\beta Re_s^{\alpha}$$', 'PF2 $$\beta Re_s^{\alpha}$$'};
  legend(labels6_a,'Location', 'NorthWest', 'Interpreter', 'Latex');
  textbox18_a = annotation('textbox', [0.4, 0.675, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox18_a.FontSize = 16;

  ax2 = axes('Position', ax222(2, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax2, 'on');
    fplot(ax2, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2)
    plot(ax2, NB1.Re_s, NB1.G, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
    plot(ax2, NB2.Re_s, NB2.G, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
    plot(ax2, NB3.Re_s, NB3.G, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
    fplot(ax2, @(Re) (NBall.powerfit.b).*(Re).^(NBall.powerfit.m), [71 10000],'-', 'Color', NBall.color,'Linewidth', 2)
  labels6_b = {'$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$', NB1.label, NB2.label, NB3.label, 'NB $$\beta Re_s^{\alpha}$$'};
  legend(labels6_b,'Location', 'NorthWest', 'Interpreter', 'Latex');
  textbox18_b = annotation('textbox', [0.9, 0.675, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox18_b.FontSize = 16;

  ax3 = axes('Position', ax222(3, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax3, 'on');
    fplot(ax3, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2)
    plot(ax3, UB1.Re_s, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS)
    plot(ax3, UB2.Re_s, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS)
    fplot(ax3, @(Re) (UBall.powerfit.b).*(Re).^(UBall.powerfit.m), [71 10000],'-', 'Color', UBall.color,'Linewidth', 2)
  ylabel('$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  labels6_c = {'$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$', UB1.label, UB2.label, 'UB $$\beta Re_s^{\alpha}$$'};
  legend(labels6_c,'Location', 'NorthWest', 'Interpreter', 'Latex');
  textbox18_c = annotation('textbox', [0.4, 0.35, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none');
  textbox18_c.FontSize = 16;

  ax4 = axes('Position', ax222(4, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax4, 'on');
    fplot(ax4, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2)
    plot(ax4, XB1.Re_s, XB1.G, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS)
    plot(ax4, XB2.Re_s, XB2.G, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS)
    fplot(ax4, @(Re) (XBall.powerfit.b).*(Re).^(XBall.powerfit.m), [71 10000],'-', 'Color', XBall.color,'Linewidth', 2)
  labels6_d = {'$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$', XB1.label, XB2.label, 'XB $$\beta Re_s^{\alpha}$$'};
  legend(labels6_d,'Location', 'NorthWest', 'Interpreter', 'Latex');
  textbox18_d = annotation('textbox', [0.9, 0.35, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none');
  textbox18_d.FontSize = 16;

  ax5 = axes('Position', ax222(5, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax5, 'on');
    fplot(ax5, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2)
    for i=1:length(FB1.exp)
      plot(ax5, FB1.exp(i).Re_s, FB1.exp(i).G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
    end
    fplot(ax5, @(Re) (FB1.powerfit.b).*(Re).^(FB1.powerfit.m), [71 10000],'-', 'Color', FB1.color,'Linewidth', 2, 'DisplayName', 'FB1 $$\beta Re_s^{\alpha}$$');
  ylabel('$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  textbox18_e = annotation('textbox', [0.4, 0.025, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'e)', 'LineStyle', 'none');
  textbox18_e.FontSize = 16;

  ax6 = axes('Position', ax222(6, :));
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax6, 'on');
    fplot(ax6, @(Re) G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2)
    for i=1:length(FB2.exp)
      plot(ax6, FB2.exp(i).Re_s, FB2.exp(i).G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', 1, 'MarkerSize', 5.0, 'MarkerSize', 5)
    end
    fplot(ax6, @(Re) (FB2.powerfit.b).*(Re).^(FB2.powerfit.m), [71 10000],'-', 'Color', FB2.color,'Linewidth', 2, 'DisplayName', 'FB2 $$\beta Re_s^{\alpha}$$');
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  textbox18_f = annotation('textbox', [0.9, 0.025, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'f)', 'LineStyle', 'none');
  textbox18_f.FontSize = 16;

  axis([ax1 ax2 ax3 ax4 ax5 ax6],[1e-2 2e4 1e-1 1e8])

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   7  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

figs(7) = figure('Name', 'FB_T_vs_omega' ,'Renderer', 'painters', 'Position',  [pos22 dim2]);
  tile_object = tiledlayout(1, 2);
  ax1 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax1, 'on');
  legend_set = [];
    errorbar(ax1, EC000.omega, EC000.mu_torque, EC000.sigma_torque, EC000.specs, 'Color', EC000.color, 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', EC000.label)
    errorbar(ax1, EC050.omega, EC050.mu_torque, EC050.sigma_torque, EC050.specs, 'Color', EC050.color, 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', EC050.label)
    errorbar(ax1, EC075.omega, EC075.mu_torque, EC075.sigma_torque, EC075.specs, 'Color', EC075.color, 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', EC075.label)
    errorbar(ax1, EC100.omega, EC100.mu_torque, EC100.sigma_torque, EC100.specs, 'Color', EC100.color, 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', EC100.label)
    for i = 1:length(FB1.exp)
      errorbar(ax1, FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.exp(i).sigma_torque, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
    end
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
  legend('Show', 'Location', 'SouthEast')
  textbox7_a = annotation('textbox', [0.1, 0.2, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox7_a.FontSize = 16;

  ax2 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax2, 'on');
  for i = 1:length(FB2.exp)
      errorbar(ax2, FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.exp(i).sigma_torque, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
  end
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
  legend('Show', 'Location', 'SouthEast')
  textbox7_b = annotation('textbox', [0.6, 0.2, 0.1, 0.1],   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox7_b.FontSize = 16;
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';
  axis([ax1 ax2],[1e-2 2e2 1e-8 1e-2])

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% ----------------------------   SUB FIG 1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

sub_figs(1) = figure('Name', 'FB_muapp_vs_omega' ,'Renderer', 'painters', 'Position', [900 100 525 250]);
  tile_object = tiledlayout(1, 2);
  ax1 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax1, 'on');
  legend_set = [];
    for i = 1:length(FB1.exp)
      plot(ax1, FB1.exp(i).omega, FB1.exp(i).appmu, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
    end
  xlim([1e-2 2e2])
  ylim([1e-2 2e4])
  ylabel('$$\mu_{app} [Pa.s]$$', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_{i} [rad/s]$$', 'Interpreter', 'LaTeX','FontSize',12)
  textbox13_a = annotation('textbox', [0.1, 0.2, 0.1, 0.1],  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox13_a.FontSize = 16;

  ax2 = nexttile;
  box on;
  set(gca, 'YScale', 'log')
  set(gca, 'XScale', 'log')
  hold(ax2, 'on');
    for i = 1:length(FB2.exp)
        plot(ax2, FB2.exp(i).omega, FB2.exp(i).appmu, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
    end
  xlim([1e-2 2e2])
  ylim([1e-2 2e4])
  ylabel('$$\mu_{app} [Pa.s]$$', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_{i} [rad/s]$$', 'Interpreter', 'LaTeX','FontSize',12)
  textbox13_b = annotation('textbox', [0.6, 0.2, 0.1, 0.1],  'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox13_b.FontSize = 16;
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% ----------------------------   SUB FIG 2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------
sub_figs(2) = figure('Name', 'FB_mup_vs_q' ,'Renderer', 'painters', 'Position', [900 100 525 250]);
  tile_object = tiledlayout(1, 2);
  ax1 = nexttile;
  box on;
  hold(ax1, 'on');
  legend_set = [];
    for i = 1:length(FB1.exp)
        plot(ax1, FB1.exp(i).q, FB1.exp(i).mu_p, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
    end
  % xlim([0 16])
  xlim([0 2])
  ylim([0 1])
  ylabel('$$\tilde{\mu}_{p} [Pa.s]$$', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  textbox13_a = annotation('textbox', [0.1, 0.85, 0.1, 0.1],  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox13_a.FontSize = 16;
  legend('Show', 'Location', 'NorthEast', 'Interpreter', 'LaTeX')

  ax2 = nexttile;
  box on;
  hold(ax2, 'on');
    for i = 1:length(FB2.exp)
        plot(ax2, FB2.exp(i).q, FB2.exp(i).mu_p, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
    end
  xlim([0 16])
  ylim([0 1])
  ylabel('$$\tilde{\mu}_{p} [Pa.s]$$', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  legend('Show', 'Location', 'NorthEast', 'Interpreter', 'LaTeX')
  textbox13_b = annotation('textbox', [0.6, 0.85, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox13_b.FontSize = 16;
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% ----------------------------   SUB FIG 3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

sub_figs(3) = figure('Name', 'FB_tauy_vs_q' ,'Renderer', 'painters', 'Position', [900 100 525 250]);
  tile_object = tiledlayout(1, 2);
  ax1 = nexttile;
  box on;
  hold(ax1, 'on');
  legend_set = [];
    for i = 1:length(FB1.exp)
        plot(ax1, FB1.exp(i).q, (FB1.exp(i).tau_y)/(FB1.exp(i).tau_static), FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
    end
  xlim([0 2])
  ylim([0 1.2])
  ylabel('$$\frac{\tau_y}{\tau_{q = 0}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  textbox13_a = annotation('textbox', [0.1, 0.85, 0.1, 0.1],  'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none');
  textbox13_a.FontSize = 16;
  legend('Show', 'Location', 'NorthEast', 'Interpreter', 'LaTeX')

  ax2 = nexttile;
  box on;
  hold(ax2, 'on');
    for i = 1:length(FB2.exp)
        plot(ax2, FB2.exp(i).q, (FB2.exp(i).tau_y)/(FB2.exp(i).tau_static), FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
    end
  ylim([0 1.2])
  ylabel('$$\frac{\tau_y}{\tau_{q = 0}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  legend('Show', 'Location', 'NorthEast', 'Interpreter', 'LaTeX')
  textbox13_b = annotation('textbox', [0.6, 0.85, 0.1, 0.1],  'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
  textbox13_b.FontSize = 16;
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';


save_dir = '~/Desktop/JFM1_15Jan_2021_plots/';

for i=1:length(figs)
  fig_it = figs(i);
  saveas(fig_it, [save_dir fig_it.Name], 'pdf');
end
