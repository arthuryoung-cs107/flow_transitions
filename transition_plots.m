classdef transition_plots < main_plots
    properties

    end
    methods
        function obj = transition_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@main_plots(write_figs_, write_all_figs_, figs_to_write_);

        end
        function fig_out = FB1_FB2_T_alphaT_vs_omega(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i=1:length(FB1.exp)
              legend_set_a(i) = plot(axs(1), FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
              legend_set_b(i) = plot(axs(2), FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end
            for i=1:length(FB1.exp)
              legend_set_c(i) = plot(axs(3), FB1.exp(i).omega, FB1.exp(i).alpha_T, [FB1.specs '-'], 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
              legend_set_d(i) = plot(axs(4), FB2.exp(i).omega, FB2.exp(i).alpha_T, [FB2.specs '-'], 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            set(axs,'XScale', 'log');
            set(axs(1:2),'YScale', 'log');

            ylabel(axs(1), '$$T$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3), '$$\alpha_T$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:4), '$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns',2);
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns',2);
            legend(axs(3), legend_set_c,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',2);
            legend(axs(4), legend_set_d,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',2);

            axis(axs(1:2), [obj.omega_range obj.torque_range])
            axis(axs(3:4), [obj.omega_range -0.5 obj.alpha_range(2)])

            fig_out = AYfig_;
        end
        function fig_out = ALL_G_vs_Res_fits(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall)
            AYfig_.init_tiles([2,3]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i = 1:6
                fplot(axs(i), @(Re) obj.G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
            end

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.Re_s, PF1.G, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = plot(axs(1), PFR.Re_s, PFR.G, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            fplot(axs(1), @(Re) (PF1.powerfit.b).*(Re).^(PF1.powerfit.m), [71 10000],'-', 'Color', PF1.color,'Linewidth', 2)
            fplot(axs(1), @(Re) (PFR.powerfit.b).*(Re).^(PFR.powerfit.m), [71 10000],'-', 'Color', PFR.color,'Linewidth', 2)

            legend_set_b(1) = plot(axs(2), NB1.Re_s, NB1.G, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = plot(axs(2), NB2.Re_s, NB2.G, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = plot(axs(2), NB3.Re_s, NB3.G, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            fplot(axs(2), @(Re) (NBall.powerfit.b).*(Re).^(NBall.powerfit.m), [71 10000],'-', 'Color', NBall.color,'Linewidth', 2, 'DisplayName', 'NB $$\beta Re_s^{\alpha}$$')

            legend_set_c(1) = plot(axs(3), UB1.Re_s, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_c(2) = plot(axs(3), UB2.Re_s, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
            fplot(axs(3), @(Re) (UBall.powerfit.b).*(Re).^(UBall.powerfit.m), [71 10000],'-', 'Color', UBall.color,'Linewidth', 2, 'DisplayName', 'UB $$\beta Re_s^{\alpha}$$')

            legend_set_d(1) = plot(axs(4), XB1.Re_s, XB1.G, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_d(2) = plot(axs(4), XB2.Re_s, XB2.G, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);
            fplot(axs(4), @(Re) (XBall.powerfit.b).*(Re).^(XBall.powerfit.m), [71 10000],'-', 'Color', XBall.color,'Linewidth', 2, 'DisplayName', 'XB $$\beta Re_s^{\alpha}$$')

            for i=1:length(FB1.exp)
              legend_set_e(i) = plot(axs(5), FB1.exp(i).Re_s, FB1.exp(i).G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            fplot(axs(5), @(Re) (FB1.powerfit.b).*(Re).^(FB1.powerfit.m), [71 10000],'-', 'Color', FB1.color,'Linewidth', 2,'DisplayName', 'FB1 $$\beta Re_s^{\alpha}$$');

            for i=1:length(FB2.exp)
              legend_set_f(i) = plot(axs(6), FB2.exp(i).Re_s, FB2.exp(i).G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end
            fplot(axs(6), @(Re) (FB2.powerfit.b).*(Re).^(FB2.powerfit.m), [71 10000],'-', 'Color', FB2.color,'Linewidth', 2, 'DisplayName', 'FB2 $$\beta Re_s^{\alpha}$$');

            % textbox_a = annotation('textbox', obj.textbox_pos222_a, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_b = annotation('textbox', obj.textbox_pos222_b, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_c = annotation('textbox', obj.textbox_pos222_c, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_d = annotation('textbox', obj.textbox_pos222_d, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_e = annotation('textbox', obj.textbox_pos222_e, 'Interpreter', 'LaTeX', 'String', 'e)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_f = annotation('textbox', obj.textbox_pos222_f, 'Interpreter', 'LaTeX', 'String', 'f)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:3:4), '$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(4:6), '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(4), legend_set_d,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(5), legend_set_e,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(6), legend_set_f,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs,obj.Res_G_range)

            fig_out = AYfig_;
        end
    end
end
