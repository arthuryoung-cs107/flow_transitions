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

        omega_tau_range = [1e-2 1e2 1e-9 1e-2];
        Res_G_range = [1e-2 1e4 1e0 1e8];
        Res_alpha_range = [1e0 5e5 -1 2.1];
        Res_cf_range = [1e-1 1e4 1e-4 1e4];
        Res_Grat_range = [1e-1 2e4 0.8 1e2];
        omega_appmu_range = [1e-2 2e2 6e-2 1e7];

        omega_range = [1e-2, 2e2];
        torque_range = [1e-7, 1e-2];
        Re_s_range = [1e-4, 2e4];
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
        posdimfull = [1 1 1728 1000];

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

            %% handling pure fluid plot component
            hold(axs(1), 'on');
            legend_set_a(1) = errorbar(axs(1), PF1.omega, PF1.mu_torque, PF1.sigma_torque, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = errorbar(axs(1), PFR.omega, PFR.mu_torque, PFR.sigma_torque, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            textbox_a = annotation('textbox', obj.textbox_pos2_a_NW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);

            %% handling neutrally buoyant stuff
            hold(axs(2), 'on');
            legend_set_b(1) = errorbar(axs(2), NB1.omega_full, NB1.mu_torque_full, NB1.sigma_torque_full, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = errorbar(axs(2), NB2.omega, NB2.mu_torque, NB2.sigma_torque, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = errorbar(axs(2), NB3.omega, NB3.mu_torque, NB3.sigma_torque, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_NW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            ylabel(axs(1), '$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)


            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');

            axis(axs,obj.omega_tau_range)
            set(axs, 'XScale', 'log', 'YScale', 'log');

            fig_out = AYfig_.fig;
        end
        function fig_out = PF_NB_Grat_vs_Res(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, RM10, RM20)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;

            %% pre KD-effective viscosity
            hold(axs(1), 'on');
            plot(axs(1), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
            plot(axs(1), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
            plot(axs(1), NB1.Re_s_noKD, NB1.G_rat_noKD, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
            plot(axs(1), NB2.Re_s_noKD, NB2.G_rat_noKD, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
            plot(axs(1), NB3.Re_s_noKD, NB3.G_rat_noKD, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
            textbox_a = annotation('textbox', obj.textbox_pos2_a_SW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);

            %% post KD-effective viscosity

            hold(axs(2), 'on');
            legend_set_b(1) = plot(axs(2), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_b(2) = plot(axs(2), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            legend_set_b(3) = plot(axs(2), NB1.Re_s, NB1.G_rat, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(4) = plot(axs(2), NB2.Re_s, NB2.G_rat, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(5) = plot(axs(2), NB3.Re_s, NB3.G_rat, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            legend_set_b(6) = plot(axs(2), RM10.Re_s, RM10.G_rat, RM10.specs, 'Color', RM10.color, 'LineWidth', RM10.LW, 'DisplayName', RM10.label);
            legend_set_b(7) = plot(axs(2), RM20.Re_s, RM20.G_rat, RM20.specs, 'Color', RM20.color, 'LineWidth', RM20.LW, 'DisplayName', RM20.label);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_SW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            ylabel(axs(1), '$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(1), '$$Re_s \textrm{ (PRE-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel('$$Re_s \textrm{ (POST-KD effective viscosity)}$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex');

            axis(axs,obj.Res_Grat_range)
            set(axs, 'XScale', 'log', 'YScale', 'log');

            fig_out = AYfig_.fig;
        end
        function fig_out = PF_NB_cf_alpha_vs_Res(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, RK, RV,LSa, LSb)
            AYfig_.ax_sub = set_alpha_cf_vs_Res_axes(AYfig_,obj.ax21);
            axs = AYfig_.ax_sub;

            %% PF alpha vs Res
            hold(axs(1), 'on');
            plot(axs(1), RV.Re_s, RV.alpha, RV.specs,'Color', RV.color, 'LineWidth', RV.LW, 'MarkerSize', RV.MS)
            plot(axs(1), LSa.Re_s, LSa.alpha, LSa.specs,'Color', LSa.color, 'LineWidth', LSa.LW)
            plot(axs(1), LSb.Re_s, LSb.alpha, LSb.specs,'Color', LSb.color, 'LineWidth', LSb.LW)
            plot(axs(1), PF1.Re_s_alpha, PF1.alpha, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
            plot(axs(1), PFR.Re_s_alpha, PFR.alpha, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
            textbox_a = annotation('textbox', obj.textbox_pos21_a_SW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize',16);
            labels_set_a = {RV.label, LSa.label, LSb.label, PF1.label, PFR.label};

            %% NB alpha vs Res
            hold(axs(2), 'on');
            plot(axs(2), NB1.Re_s_alpha, NB1.alpha, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
            plot(axs(2), NB2.Re_s_alpha, NB2.alpha, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
            plot(axs(2), NB3.Re_s_alpha, NB3.alpha, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
            textbox_b = annotation('textbox', obj.textbox_pos21_b_SW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize',16);


            %% cf vs Res
            hold(axs(3), 'on');
            fplot(axs(3), @(Re) 1./(obj.eta.*(Re)), [1e-1 5000], '-. k', 'Linewidth', 2)
            plot(axs(3), PF1.Re_s, PF1.cf, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS)
            plot(axs(3), PFR.Re_s, PFR.cf, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS)
            plot(axs(3), NB1.Re_s, NB1.cf, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS)
            plot(axs(3), NB2.Re_s, NB2.cf, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS)
            plot(axs(3), NB3.Re_s, NB3.cf, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS)
            plot(axs(3), RK.Re_s, RK.cf, RK.specs,'Color', RK.color, 'LineWidth', RK.LW)
            textbox_c = annotation('textbox', obj.textbox_pos21_c_SW,   'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize',16);
            labels_set_c = {'$$\frac{1}{\eta Re_s}$$', PF1.label, PFR.label, NB1.label, NB2.label, NB3.label, RK.label};

            ylabel(axs(1),'$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3),'$$c_f$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), labels_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(2), labels_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(3), labels_set_c,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 3);


            axis(axs,obj.Res_Grat_range)
            set(axs, 'XScale', 'log', 'YScale', 'log');

            fig_out = AYfig_.fig;

                axa = axes('Position', ax21(2, :));
                box on;
                set(gca, 'XScale', 'log')
                hold(axa, 'on');
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
                xlabel('$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)
                labels = {NB1.label, NB2.label, NB3.label};
                legend(labels,'Location', 'SouthEast', 'Interpreter', 'Latex');
                textbox_b = annotation('textbox', textbox_pos21_b_SW,   'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none');
                textbox_b.FontSize = 16;

                axis(axc,Res_cf_range)
                axis([axa axb], Res_alpha_range)
        end
        function write_figures(obj, figs_, save_dir_, save_type_)
            if (obj.write_figs)
              if (obj.write_all_figs)
                obj.figs_to_write = 1:length(figs_);
              end
              AYfig.save_figs(figs_, obj.figs_to_write, save_type_, save_dir_);
            end
        end
    end
end

function ax_out = make_alpha_cf_vs_Res_axes(AYfig_, ax21_)
    ax_out = gobjects(3, 1);
    if (nargin==1)
        figure(AYfig_.fig.Number);
        axa = subplot(AYfig_.fig, 2,2,1);
        axb = subplot(AYfig_.fig, 2,2,2);
        axc = subplot(AYfig_.fig, 2,2,[3,4]);
    else
        axa = axes(AYfig_.fig, 'Position', ax21(2, :));
        axb = axes(AYfig_.fig, 'Position', ax21(3, :));
        axc = axes(AYfig_.fig, 'Position', ax21(1, :));
    end
    ax_out(1) = axa;
    ax_out(2) = axb;
    ax_out(3) = axc;
end
