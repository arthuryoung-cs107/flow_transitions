clear
close all

run load_data.m
run set_figdims.m
run compare_alpha_FIGSPECS.m

fig_num = 0;

% output information
write_figs = true;
write_all_figs = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';
%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(2, 3);
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

  axa = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axa, 'on');
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',12)
  axis(axa, [omega_range torque_range])

  axb = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axb, 'on');
  ylabel('$$G$$ [dimensions]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  axis(axb, [Re_s_range, G_range])

  axc = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axc, 'on');
  ylabel('$$c_f$$ [dimensions]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  axis(axc, [Re_s_range cf_range])

  axd = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axd, 'on');
  ylabel('$$\alpha_{T}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',12)

  axe = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axe, 'on');
  ylabel('$$\alpha_{G}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

  axf = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axf, 'on');
  ylabel('$$\alpha_{cf}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

  axis([axe axf], [Re_s_range alpha_range])
  axis(axd, [omega_range alpha_range])
    for i = 1:length(FB1.exp)
      plot(axa, FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      plot(axb, FB1.exp(i).Re_s, FB1.exp(i).G, [FB1.exp(i).specs, '-'],'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      plot(axc, FB1.exp(i).Re_s, FB1.exp(i).cf, [FB1.exp(i).specs, '-'],'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      plot(axd, FB1.exp(i).omega, FB1.exp(i).alpha_T, [FB1.exp(i).specs, '-'],'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      plot(axe, FB1.exp(i).Re_s, FB1.exp(i).alpha_G, [FB1.exp(i).specs, '-'],'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      plot(axf, FB1.exp(i).Re_s, FB1.exp(i).alpha_cf, [FB1.exp(i).specs, '-'],'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
    end

  legend(axa, 'Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(2, 3);
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

  axa = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axa, 'on');
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',12)
  axis(axa, [omega_range torque_range])

  axb = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axb, 'on');
  ylabel('$$G$$ [dimensions]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  axis(axb, [Re_s_range, G_range])

  axc = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axc, 'on');
  ylabel('$$c_f$$ [dimensions]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
  axis(axc, [Re_s_range cf_range])

  axd = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axd, 'on');
  ylabel('$$\alpha_{T}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',12)

  axe = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axe, 'on');
  ylabel('$$\alpha_{G}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

  axf = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axf, 'on');
  ylabel('$$\alpha_{cf}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

  axis([axe axf], [Re_s_range alpha_range])
  axis(axd, [omega_range alpha_range])
    for i = 1:length(FB2.exp)
      plot(axa, FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.exp(i).specs,'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
      plot(axb, FB2.exp(i).Re_s, FB2.exp(i).G, [FB2.exp(i).specs, '-'],'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
      plot(axc, FB2.exp(i).Re_s, FB2.exp(i).cf, [FB2.exp(i).specs, '-'],'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
      plot(axd, FB2.exp(i).omega, FB2.exp(i).alpha_T, [FB2.exp(i).specs, '-'],'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
      plot(axe, FB2.exp(i).Re_s, FB2.exp(i).alpha_G, [FB2.exp(i).specs, '-'],'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
      plot(axf, FB2.exp(i).Re_s, FB2.exp(i).alpha_cf, [FB2.exp(i).specs, '-'],'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
    end

  legend(axa, 'Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(2, 5);
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

    for i=1:9
      ax(i) = nexttile;
      box on;
      set(gca, 'XScale', 'log')
      hold(ax(i), 'on');
      % xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
      xlim([1e0, 2e4]);

      yyaxis(ax(i),'left')
      set(gca, 'YScale', 'log')
      % ylabel('$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
      ylim([1e0, 3e1]);

      yyaxis(ax(i),'right')
      % ylabel('$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)
      ylim([1, 2]);

    end


    xlabel(ax(6:9), '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
    yyaxis(ax(6),'left')
    ylabel(ax(6), '$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
    yyaxis(ax(9),'right')
    ylabel(ax(9), '$$\alpha$$', 'Interpreter', 'LaTeX','FontSize',12)

    yyaxis(ax(1),'left')
      plot(ax(1), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
    yyaxis(ax(1),'right')
      plot(ax(1), PF1.Re_s, PF1.alpha, PF1.specs,'Color', PFR.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
    yyaxis(ax(2),'left')
      plot(ax(2), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PF1.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
    yyaxis(ax(2),'right')
      plot(ax(2), PFR.Re_s, PFR.alpha, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
    yyaxis(ax(3),'left')
      plot(ax(3), NB1.Re_s, NB1.G_rat, NB1.specs,'Color', PF1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
    yyaxis(ax(3),'right')
      plot(ax(3), NB1.Re_s, NB1.alpha, NB1.specs,'Color', PFR.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
    yyaxis(ax(4),'left')
      plot(ax(4), NB2.Re_s, NB2.G_rat, NB2.specs,'Color', PF1.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
    yyaxis(ax(4),'right')
      plot(ax(4), NB2.Re_s, NB2.alpha, NB2.specs,'Color', PFR.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
    yyaxis(ax(5),'left')
      plot(ax(5), NB3.Re_s, NB3.G_rat, NB3.specs,'Color', PF1.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
    yyaxis(ax(5),'right')
      plot(ax(5), NB3.Re_s, NB3.alpha, NB3.specs,'Color', PFR.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
    yyaxis(ax(6),'left')
      plot(ax(6), UB1.Re_s, UB1.G_rat, UB1.specs,'Color', PF1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
    yyaxis(ax(6),'right')
      plot(ax(6), UB1.Re_s, UB1.alpha, UB1.specs,'Color', PFR.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
    yyaxis(ax(7),'left')
      plot(ax(7), UB2.Re_s, UB2.G_rat, UB2.specs,'Color', PF1.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
    yyaxis(ax(7),'right')
      plot(ax(7), UB2.Re_s, UB2.alpha, UB2.specs,'Color', PFR.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
    yyaxis(ax(8),'left')
      plot(ax(8), XB1.Re_s, XB1.G_rat, XB1.specs,'Color', PF1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
    yyaxis(ax(8),'right')
      plot(ax(8), XB1.Re_s, XB1.alpha, XB1.specs,'Color', PFR.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
    yyaxis(ax(9),'left')
      plot(ax(9), XB2.Re_s, XB2.G_rat, XB2.specs,'Color', PF1.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);
    yyaxis(ax(9),'right')
      plot(ax(9), XB2.Re_s, XB2.alpha, XB2.specs,'Color', PFR.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
figs(fig_num) = AYfig.figure(fig_specs{fig_num});
  tile_object = tiledlayout(2, 2);
  tile_object.TileSpacing = 'compact';
  tile_object.Padding = 'compact';

  axa = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axa, 'on');
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

  axb = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  set(gca, 'YScale', 'log')
  hold(axb, 'on');
  ylabel('$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

  axc = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axc, 'on');
  ylabel('$$\alpha_{T}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

  axd = nexttile;
  box on;
  set(gca, 'XScale', 'log')
  hold(axd, 'on');
  ylabel('$$\alpha_{T}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
  xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

  axis([axa axb], [omega_range torque_range])
  axis([axc axd], [omega_range alpha_range])
    for i = 1:length(FB1.exp)
      plot(axa, FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
      plot(axc, FB1.exp(i).omega, FB1.exp(i).alpha_T, [FB1.exp(i).specs, '-'],'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label)
    end

    for i = 1:length(FB2.exp)
      plot(axb, FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.exp(i).specs,'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
      plot(axd, FB2.exp(i).omega, FB2.exp(i).alpha_T, [FB2.exp(i).specs, '-'],'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label)
    end

  legend(axa, 'Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)
  legend(axb, 'Show', 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'NumColumns', 2)

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

if (write_figs)
  if (write_all_figs)
    figs_to_write = 1:length(figs);
  end
  AYfig.save_figs(figs, figs_to_write, save_type, save_dir);
end
