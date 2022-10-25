classdef main_plots
    properties
        write_figs;
        write_all_figs;
        figs_to_write;

        eta;
        G_obs_Res_slope;

        dim2;
        dim21;
        dim222;
        dim22;
        ax21;

        r_i = 0.01208;
        r_o = 0.025;
        L = 0.036;
        tau2T = 2*pi*(0.01208)*(0.01208)*(0.036);

        omega_tau_range = [1e-2 1e2 1e-9 1e-2];
        Res_G_range = [1e-2 1e4 1e0 1e8];
        Res_alpha_range = [1e0 5e5 -1 2.1];
        Res_cf_range = [1e-1 1e4 1e-4 1e4];
        Res_Grat_range = [1e-1 2e4 0.8 1e2];
        omega_appmu_range = [1e-2 2e2 6e-2 1e7];

        omega_range = [1e-2, 2e2];
        Re_s_range = [1e-4, 2e4];

        torque_range = [1e-7, 1e-2];
        tau_range = [1e-2 1e3];
        G_range = [1e-1, 1e5];
        Grat_range = [5e-1, 3e1];
        cf_range = [1e1, 1e10];
        alpha_range = [0, 3];

        pos11 = [0 500];
        pos12 = [450 500];
        pos13 = [900 500];
        pos21 = [0 100];
        pos22 = [450 100];
        pos23 = [900 100];

        pos_spread = [0 500; 180 500; 360 500; 540 500; 720 500; 900 500; 0 100; 180 100; 360 100; 540 100; 720 100; 900 100];

        textbox_pos2_a_NW = [0.09, 0.85, 0.1, 0.1];
        textbox_pos2_b_NW = [0.575, 0.85, 0.1, 0.1];
        textbox_pos2_a_SW = [0.09, 0.175, 0.1, 0.1];
        textbox_pos2_b_SW = [0.575, 0.175, 0.1, 0.1];

        textbox_pos21_a_NW = [0.12, 0.85, 0.1, 0.1];
        textbox_pos21_b_NW = [0.55, 0.85, 0.1, 0.1];
        textbox_pos21_c_NW = [0.12, 0.07, 0.1, 0.1];
        textbox_pos21_a_SW = [0.12, 0.675, 0.1, 0.1];
        textbox_pos21_b_SW = [0.55, 0.675, 0.1, 0.1];
        textbox_pos21_c_SW = [0.12, 0.1, 0.1, 0.1];

        textbox_pos222_a = [0.4, 0.625, 0.1, 0.1];
        textbox_pos222_b = [0.9, 0.625, 0.1, 0.1];
        textbox_pos222_c = [0.4, 0.32, 0.1, 0.1];
        textbox_pos222_d = [0.9, 0.32, 0.1, 0.1];
        textbox_pos222_e = [0.4, 0.008, 0.1, 0.1];
        textbox_pos222_f = [0.9, 0.008, 0.1, 0.1];

        textbox_pos22_a_NE = [0.425, 0.85, 0.1, 0.1];
        textbox_pos22_b_NE = [0.9, 0.85, 0.1, 0.1];
        textbox_pos22_c_NE = [0.425, 0.375, 0.1, 0.1];
        textbox_pos22_d_NE = [0.9, 0.375, 0.1, 0.1];

        dim2_short = [580 250];
        dim2_tall = [580 330];

        dim21_short = [580 330];
        dim21_tall = [580 450];

        dim222_short = [580 600];
        dim222_tall = [580 800];

        dim22_short = [580 300];
        dim22_tall = [580 500];

        dim1 = [580 325];

        % posdimfull = [1 1 1438 796];
        % posdimfull = [1 1 1728 1000];
        posdimfull = [0 0 1728 1000];

        pos_top_row = [1 551 1728 460];
        pos_bottom_row = [0 1 1728 460];

        ax21_short = [0.1 0.1 0.85 0.475; 0.1 0.675 0.4 0.3; 0.545 0.675 0.4 0.3];
        ax21_tall = [0.1 0.1 0.85 0.475; 0.1 0.675 0.4 0.3; 0.545 0.675 0.4 0.3];
    end
    methods
        function obj = main_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj.write_figs = write_figs_;
            obj.write_all_figs = write_all_figs_;
            obj.figs_to_write = figs_to_write_;

            obj.eta = obj.r_i/obj.r_o; % gap ratio
            obj.G_obs_Res_slope = (2*pi*obj.r_i*obj.r_o)/((obj.r_o-obj.r_i)^2);

            obj.dim2 = obj.dim2_short;
            obj.dim21 = obj.dim21_short;
            obj.dim222 = obj.dim222_tall;
            obj.dim22 = obj.dim22_short;

            obj.ax21 = obj.ax21_short;

        end
        function fig_out = PF_NB_T_vs_omega(obj, AYfig_, PF1, PFR, NB1, NB2, NB3)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            %% handling pure fluid plot component
            legend_set_a(1) = errorbar(axs(1), PF1.omega, PF1.mu_torque, PF1.sigma_torque, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = errorbar(axs(1), PFR.omega, PFR.mu_torque, PFR.sigma_torque, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            %% handling neutrally buoyant stuff
            legend_set_b(1) = errorbar(axs(2), NB1.omega_full, NB1.mu_torque_full, NB1.sigma_torque_full, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = errorbar(axs(2), NB2.omega, NB2.mu_torque, NB2.sigma_torque, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = errorbar(axs(2), NB3.omega, NB3.mu_torque, NB3.sigma_torque, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            textbox_a = annotation('textbox', obj.textbox_pos2_a_NW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_NW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            ylabel(axs(1), '$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');

            set(axs, 'XScale', 'log', 'YScale', 'log');

            axis(axs,obj.omega_tau_range)

            fig_out = AYfig_;
        end
        function fig_out = PF_NB_Grat_vs_Res(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, RM10, RM20)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            %% pre KD-effective viscosity
            plot(axs(1), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
            plot(axs(1), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
            plot(axs(1), NB1.Re_s_noKD, NB1.G_rat_noKD, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
            plot(axs(1), NB2.Re_s_noKD, NB2.G_rat_noKD, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
            plot(axs(1), NB3.Re_s_noKD, NB3.G_rat_noKD, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)

            %% post KD-effective viscosity
            legend_set_b(1) = plot(axs(2), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_b(2) = plot(axs(2), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            legend_set_b(3) = plot(axs(2), NB1.Re_s, NB1.G_rat, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(4) = plot(axs(2), NB2.Re_s, NB2.G_rat, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(5) = plot(axs(2), NB3.Re_s, NB3.G_rat, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            legend_set_b(6) = plot(axs(2), RM10.Re_s, RM10.G_rat, RM10.specs, 'Color', RM10.color, 'LineWidth', RM10.LW, 'DisplayName', RM10.label);
            legend_set_b(7) = plot(axs(2), RM20.Re_s, RM20.G_rat, RM20.specs, 'Color', RM20.color, 'LineWidth', RM20.LW, 'DisplayName', RM20.label);

            textbox_a = annotation('textbox', obj.textbox_pos2_a_SW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_SW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            ylabel(axs(1), '$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(1), '$$Re_s \textrm{ (PRE-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel('$$Re_s \textrm{ (POST-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex');

            set(axs, 'XScale', 'log', 'YScale', 'log');

            axis(axs,obj.Res_Grat_range)

            fig_out = AYfig_;
        end
        function fig_out = PF_NB_cf_alpha_vs_Res(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, RV, LSa, LSb, RK)
            AYfig_.ax_sub = set_alpha_cf_vs_Res_axes(AYfig_,obj.ax21);
            axs = AYfig_.ax_sub;
            hold(axs, 'on');
            box(axs,'on');

            %% PF alpha vs Res
            legend_set_a(1) = plot(axs(1), RV.Re_s, RV.alpha, RV.specs,'Color', RV.color, 'LineWidth', RV.LW, 'MarkerSize', RV.MS, 'DisplayName', RV.label);
            legend_set_a(2) = plot(axs(1), LSa.Re_s, LSa.alpha, LSa.specs,'Color', LSa.color, 'LineWidth', LSa.LW, 'DisplayName', LSa.label);
            legend_set_a(3) = plot(axs(1), LSb.Re_s, LSb.alpha, LSb.specs,'Color', LSb.color, 'LineWidth', LSb.LW, 'DisplayName', LSb.label);
            legend_set_a(4) = plot(axs(1), PF1.Re_s_alpha, PF1.alpha, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(5) = plot(axs(1), PFR.Re_s_alpha, PFR.alpha, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            %% NB alpha vs Res
            legend_set_b(1) = plot(axs(2), NB1.Re_s_alpha, NB1.alpha, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = plot(axs(2), NB2.Re_s_alpha, NB2.alpha, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = plot(axs(2), NB3.Re_s_alpha, NB3.alpha, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            %% cf vs Res
            legend_set_c(1) = fplot(axs(3), @(Re) 1./(obj.eta.*(Re)), [1e-1 5000], '-. k', 'Linewidth', 2, 'DisplayName', '$$\frac{1}{\eta Re_s}$$');
            legend_set_c(2) = plot(axs(3), PF1.Re_s, PF1.cf, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_c(3) = plot(axs(3), PFR.Re_s, PFR.cf, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            legend_set_c(4) = plot(axs(3), NB1.Re_s, NB1.cf, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_c(5) = plot(axs(3), NB2.Re_s, NB2.cf, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_c(6) = plot(axs(3), NB3.Re_s, NB3.cf, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            legend_set_c(7) = plot(axs(3), RK.Re_s, RK.cf, RK.specs,'Color', RK.color, 'LineWidth', RK.LW, 'DisplayName', RK.label);

            textbox_a = annotation('textbox', obj.textbox_pos21_a_SW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize',16);
            textbox_b = annotation('textbox', obj.textbox_pos21_b_SW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize',16);
            textbox_c = annotation('textbox', obj.textbox_pos21_c_SW,   'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize',16);

            ylabel(axs(1),'$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3),'$$c_f$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 3);

            set(axs, 'XScale', 'log')
            set(axs(3), 'YScale', 'log')

            axis(axs(1:2), obj.Res_alpha_range)
            axis(axs(3),obj.Res_cf_range)

            fig_out = AYfig_;
        end
        function fig_out = UB_XB_T_vs_omega(obj, AYfig_, UB1, UB2, XB1, XB2)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            %% handling pure fluid plot component
            legend_set_a(1) = errorbar(axs(1), UB1.omega, UB1.mu_torque, UB1.sigma_torque, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_a(2) = errorbar(axs(1), UB2.omega, UB2.mu_torque, UB2.sigma_torque, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            %% handling neutrally buoyant stuff
            legend_set_b(1) = errorbar(axs(2), XB1.omega, XB1.mu_torque, XB1.sigma_torque, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_b(2) = errorbar(axs(2), XB2.omega, XB2.mu_torque, XB2.sigma_torque, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            textbox_a = annotation('textbox', obj.textbox_pos2_a_NW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_NW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            ylabel(axs(1), '$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');

            set(axs, 'XScale', 'log', 'YScale', 'log');

            axis(axs,obj.omega_tau_range)

            fig_out = AYfig_;
        end
        function fig_out = UB_XB_cf_alpha_vs_Res(obj, AYfig_, UB1, UB2, XB1, XB2)
            AYfig_.ax_sub = set_alpha_cf_vs_Res_axes(AYfig_,obj.ax21);
            axs = AYfig_.ax_sub;
            hold(axs, 'on');
            box(axs,'on');

            %% PF alpha vs Res
            legend_set_a(1) = plot(axs(1), UB1.Re_s_alpha, UB1.alpha, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_a(2) = plot(axs(1), UB2.Re_s_alpha, UB2.alpha, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            %% NB alpha vs Res
            legend_set_b(1) = plot(axs(2), XB1.Re_s_alpha, XB1.alpha, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS);
            legend_set_b(2) = plot(axs(2), XB2.Re_s_alpha, XB2.alpha, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS);

            %% cf vs Res
            legend_set_c(1) = fplot(axs(3), @(Re) 1./(obj.eta.*(Re)), [1e-1 5000], '-. k', 'Linewidth', 2, 'DisplayName', '$$\frac{1}{\eta Re_s}$$');
            legend_set_c(2) = plot(axs(3), UB1.Re_s, UB1.cf, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_c(3) = plot(axs(3), UB2.Re_s, UB2.cf, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
            legend_set_c(4) = plot(axs(3), XB1.Re_s, XB1.cf, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_c(5) = plot(axs(3), XB2.Re_s, XB2.cf, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            textbox_a = annotation('textbox', obj.textbox_pos21_a_NW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize',16);
            textbox_b = annotation('textbox', obj.textbox_pos21_b_NW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize',16);
            textbox_c = annotation('textbox', obj.textbox_pos21_c_SW,   'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize',16);

            ylabel(axs(1),'$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3),'$$c_f$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'NorthEast', 'Interpreter', 'Latex');

            set(axs, 'XScale', 'log')
            set(axs(3), 'YScale', 'log')

            axis(axs(1:2), obj.Res_alpha_range)
            axis(axs(3),obj.Res_cf_range)

            fig_out = AYfig_;
        end
        function fig_out = ALL_G_vs_Res(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall)
            AYfig_.init_tiles([3,2]);
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
              legend_set_e(i) = plot(axs(5), FB1.exp(i).get_Re_s, FB1.exp(i).get_G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            fplot(axs(5), @(Re) (FB1.powerfit.b).*(Re).^(FB1.powerfit.m), [71 10000],'-', 'Color', FB1.color,'Linewidth', 2,'DisplayName', 'FB1 $$\beta Re_s^{\alpha}$$');

            for i=1:length(FB2.exp)
              % legend_set_f(i) = plot(axs(6), FB2.exp(i).Re_s, FB2.exp(i).G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
              legend_set_f(i) = plot(axs(6), FB2.exp(i).get_Re_s, FB2.exp(i).get_G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end
            fplot(axs(6), @(Re) (FB2.powerfit.b).*(Re).^(FB2.powerfit.m), [71 10000],'-', 'Color', FB2.color,'Linewidth', 2, 'DisplayName', 'FB2 $$\beta Re_s^{\alpha}$$');

            textbox_a = annotation('textbox', obj.textbox_pos222_a, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos222_b, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', obj.textbox_pos222_c, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', obj.textbox_pos222_d, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_e = annotation('textbox', obj.textbox_pos222_e, 'Interpreter', 'LaTeX', 'String', 'e)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_f = annotation('textbox', obj.textbox_pos222_f, 'Interpreter', 'LaTeX', 'String', 'f)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:2:5), '$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(5:6), '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(4), legend_set_d,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(5), legend_set_e,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(6), legend_set_f,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs,obj.Res_G_range)

            fig_out = AYfig_;
        end
        function fig_out = FB_T_vs_omega(obj, AYfig_, FB1, FB2, EC000, EC050, EC075, EC100)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            legend_set_a(1) = errorbar(axs(1), EC000.omega, EC000.mu_torque, EC000.sigma_torque, EC000.specs, 'Color', EC000.color, 'LineWidth', EC000.LW, 'MarkerSize', EC000.MS, 'DisplayName', EC000.label);
            legend_set_a(2) = errorbar(axs(1), EC050.omega, EC050.mu_torque, EC050.sigma_torque, EC050.specs, 'Color', EC050.color, 'LineWidth', EC050.LW, 'MarkerSize', EC050.MS, 'DisplayName', EC050.label);
            legend_set_a(3) = errorbar(axs(1), EC075.omega, EC075.mu_torque, EC075.sigma_torque, EC075.specs, 'Color', EC075.color, 'LineWidth', EC075.LW, 'MarkerSize', EC075.MS, 'DisplayName', EC075.label);
            legend_set_a(4) = errorbar(axs(1), EC100.omega, EC100.mu_torque, EC100.sigma_torque, EC100.specs, 'Color', EC100.color, 'LineWidth', EC100.LW, 'MarkerSize', EC100.MS, 'DisplayName', EC100.label);
            for i = 1:length(FB1.exp)
              legend_set_a(4+i) = errorbar(axs(1), FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.exp(i).sigma_torque, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_b(i) = errorbar(axs(2), FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.exp(i).sigma_torque, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', obj.textbox_pos2_a_SW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_SW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1), '$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs,obj.omega_tau_range)

            fig_out = AYfig_;
        end
        function fig_out = FB_appmu_vs_omega(obj, AYfig_, FB1, FB2, EC000, EC050, EC075, EC100)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i = 1:length(FB1.exp)
              legend_set_a(i) = plot(axs(1), FB1.exp(i).omega, FB1.exp(i).appmu, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_b(i) = plot(axs(2), FB2.exp(i).omega, FB2.exp(i).appmu, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', obj.textbox_pos2_a_NW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_NW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1), '$$\mu_{app}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs,obj.omega_appmu_range)

            fig_out = AYfig_;
        end
        function fig_out = FB_tauyrat_mup_vs_q(obj, AYfig_, FB1, FB2, EC000, EC050, EC075, EC100)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            % legend_set_a(1) = plot(axs(1), EC000.q, EC000.tau_y_rat, EC000.specs, 'Color', EC000.color, 'LineWidth', EC000.LW_L, 'MarkerSize', EC000.MS_L, 'DisplayName', EC000.label);
            % legend_set_a(2) = plot(axs(1), EC050.q, EC050.tau_y_rat, EC050.specs, 'Color', EC050.color, 'LineWidth', EC050.LW_L, 'MarkerSize', EC050.MS_L, 'DisplayName', EC050.label);
            % legend_set_a(3) = plot(axs(1), EC075.q, EC075.tau_y_rat, EC075.specs, 'Color', EC075.color, 'LineWidth', EC075.LW_L, 'MarkerSize', EC075.MS_L, 'DisplayName', EC075.label);
            % legend_set_a(4) = plot(axs(1), EC100.q, EC100.tau_y_rat, EC100.specs, 'Color', EC100.color, 'LineWidth', EC100.LW_L, 'MarkerSize', EC100.MS_L, 'DisplayName', EC100.label);
            for i = 1:length(FB1.exp)
                legend_set_a(i) = plot(axs(1), FB1.exp(i).q, (FB1.exp(i).tau_y)/(FB1.exp(i).tau_static), FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_b(i) = plot(axs(2), FB2.exp(i).q, (FB2.exp(i).tau_y)/(FB2.exp(i).tau_static), FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            for i = 1:length(FB1.exp)
                legend_set_c(i) = plot(axs(3), FB1.exp(i).q, FB1.exp(i).mu_p, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_d(i) = plot(axs(4), FB2.exp(i).q, FB2.exp(i).mu_p, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', obj.textbox_pos22_a_NE, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos22_b_NE, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', obj.textbox_pos22_c_NE, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', obj.textbox_pos22_d_NE, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);

            % set(axs(1:2),'YScale', 'log');
            set(axs,'YScale', 'log');

            set(axs(1:2), 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1])


            ylabel(axs(1), '$$\tau_y/\tau_{q = 0}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3), '$$\tilde{\mu}_{p}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:4), '$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            % legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            % legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs(1),[0 2 1e-5 1.0])
            axis(axs(2),[0 16 1e-5 1.0])
            axis(axs(3),[0 2 5e-2 1.2])
            axis(axs(4),[0 16 5e-2 1.2])

            fig_out = AYfig_;
        end
        function fig_out = FB_alpha_vs_omega_q(obj, AYfig_, FB1, FB2, EC000, EC050, EC075, EC100)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i = 1:length(FB1.exp)
                legend_set_a(i) = plot(axs(1), FB1.exp(i).omega, FB1.exp(i).alpha, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_b(i) = plot(axs(2), FB2.exp(i).omega, FB2.exp(i).alpha, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            for i = 1:length(FB1.exp)
                legend_set_c(i) = plot(axs(3), FB1.exp(i).q, FB1.exp(i).powerfit.m, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_d(i) = plot(axs(4), FB2.exp(i).q, FB2.exp(i).powerfit.m, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', obj.textbox_pos22_a_NE, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos22_b_NE, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', obj.textbox_pos22_c_NE, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', obj.textbox_pos22_d_NE, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs(1:2),'XScale', 'log');

            ylabel(axs(1), '$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3), '$$\alpha$$ fitted [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:4), '$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            % legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            % legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs(1:2), [1e-2, 2e2, -1, 3])
            axis(axs(3),[0 2 -1 3])
            axis(axs(4),[0 16 -1 3])

            fig_out = AYfig_;
        end
        function write_figures(obj, figs_, save_dir_, save_type_)
            if (obj.write_figs)
                if (obj.write_all_figs)
                    obj.figs_to_write = 1:length(figs_);
                end
                if (strcmp(save_type_, 'pdf'))
                    for i = obj.figs_to_write
                        exportgraphics(figs_(i).fig, [save_dir_ figs_(i).fig.Name '.pdf'], 'ContentType', 'vector')
                    end
                else
                    saveas(figs_(i).fig, [save_dir_ figs_(i).fig.Name], save_type_);
                end
            end
        end
    end
end

function ax_out = set_alpha_cf_vs_Res_axes(AYfig_, ax21_)
    ax_out = gobjects(3, 1);
    if (nargin==1)
        figure(AYfig_.fig.Number);
        axa = subplot(AYfig_.fig, 2,2,1);
        axb = subplot(AYfig_.fig, 2,2,2);
        axc = subplot(AYfig_.fig, 2,2,[3,4]);
    else
        axa = axes(AYfig_.fig, 'Position', ax21_(2, :));
        axb = axes(AYfig_.fig, 'Position', ax21_(3, :));
        axc = axes(AYfig_.fig, 'Position', ax21_(1, :));
    end
    ax_out(1) = axa;
    ax_out(2) = axb;
    ax_out(3) = axc;
end
