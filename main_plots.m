classdef main_plots
    properties (Constant)
        pos1_graphical_abstract = [40 0];
        dim1_graphical_abstract = [2.4 2.4];
        posdim1_graphical_abstract = [main_plots.pos1_graphical_abstract 4*main_plots.dim1_graphical_abstract];
    end
    properties
        write_figs;
        write_all_figs;
        figs_to_write;

        eta;
        G_obs_Res_slope;
        G_obs_Reb_slope;

        dim2;
        dim21;
        dim222;
        dim22;
        ax21;
        dim32 = [580 450];

        r_i = 0.01208;
        r_o = 0.025;
        L = 0.036;
        tau2T = 2*pi*(0.01208)*(0.01208)*(0.036);

        omega_tau_range = [1e-2 1e2 1e-9 1e-2];
        omega_T_range = [1e-2 1.26e2 1e-9 1e-2];
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
        % pos_spread = AYfig.get_pos_spread(1,[2 6]);

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
        dim22_med = [580 400];
        dim22_tall = [580 500];

        dim1 = [580 325];

        dim32_tall = [580 500];

        % posdimfull = [1 1 1438 796];
        % posdimfull = [1 1 1728 1000];
        posdimfull = [0 0 1728 1000];

        pos_top_row = [1 551 1728 460];
        pos_bottom_row = [0 1 1728 460];

        posdim_SW = [0 1 864 460];

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
            obj.G_obs_Reb_slope = (4*pi*obj.r_i*obj.r_o*obj.r_o)/((obj.r_o-obj.r_i)*(obj.r_o-obj.r_i)*(obj.r_i + obj.r_o));

            obj.dim2 = obj.dim2_short;
            obj.dim21 = obj.dim21_short;
            obj.dim222 = obj.dim222_tall;
            obj.dim22 = obj.dim22_short;

            obj.ax21 = obj.ax21_short;

            obj.pos_spread = AYfig.get_pos_spread(1,[2 6]);
        end
        function fig_out = FB_phi_vs_omegai_vs_q_slices(obj,AYfig_,FB1_phi_exp)
            axs = set_FB1_phi_vs_omegai_vs_q_slices_axes(AYfig_);
            box(axs,'on');
            % box(axs(2:end),'on');
            view(axs(1),[45 45 45]);

            o_slice = [0.01 10 100];

            poq = FB1_phi_exp.phi_o_q_mat;
            [phi_vec o_vec q_vec] = deal(poq(:,1), poq(:,2), poq(:,3));
            [poq1 poq2 poq3] = FB1_phi_exp.get_3_speeds(poq, o_slice);

            scatter3(axs(1),q_vec,o_vec,phi_vec,'*','LineWidth',0.75,'SizeData',50,'CData',phi_vec);
            scatter(axs(2),poq1(:,3),poq1(:,1),'*','LineWidth',0.75,'SizeData',50,'CData',poq1(:,1));
            scatter(axs(3),poq2(:,3),poq2(:,1),'*','LineWidth',0.75,'SizeData',50,'CData',poq2(:,1));
            scatter(axs(4),poq3(:,3),poq3(:,1),'*','LineWidth',0.75,'SizeData',50,'CData',poq3(:,1));
            set(axs(1),'YScale','log','YTick',[1e-2 1e-1 1e0 1e1 1e2])

            cbar = colorbar(axs(1));
            clims = cbar.Limits;
            cbar.Label.String = '$$ \phi $$';
            cbar.Label.Interpreter = 'Latex';
            cbar.Label.FontSize = 16;
            for iax = 1:length(axs)
                colormap(axs(iax),cool);
            end

            [qlim olim plim] = deal(axs(1).XLim, axs(1).YLim, axs(1).ZLim);

            q_mat = [qlim(2) qlim(2) qlim(1) qlim(1)];
            p_mat = [plim(1) plim(2) plim(2) plim(1)];

            fill3(axs(1), q_mat,o_slice(1)*ones(1,4),p_mat,clims(2),'LineStyle', ':', 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',2);
            fill3(axs(1), q_mat,o_slice(2)*ones(1,4),p_mat,clims(2),'LineStyle', ':', 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',2);
            fill3(axs(1), q_mat,o_slice(3)*ones(1,4),p_mat,clims(2),'LineStyle', ':', 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',2);

            axis(axs(2:end),[qlim plim])

            zlabel(axs(1),'$$\phi$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(1),'$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(1),'$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            ylabel(axs(2),'$$\phi$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(2:end),'$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            for iax = 2:length(axs)
                title(axs(iax), ['$$ \omega_i = ' num2str(o_slice(iax-1)) ' $$ rad/s'], 'Interpreter', 'LaTeX', 'FontSize', 14);
            end

            textbox_a = annotation('textbox', [0.09 0.88 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.27 0.27 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.59 0.27 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(c)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.92 0.27 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(d)}', 'LineStyle', 'none', 'FontSize', 16);

            fig_out = AYfig_;
        end
        function fig_out = PL_NUXB_Grat_vs_Res(obj, AYfig_, PFall, NBall, UBall, XBall, RMall)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            [PF1 PFR] = deal(PFall{1}, PFall{2});
            [NB1 NB2 NB3] = deal(NBall.NB1_in,NBall.NB2_in,NBall.NB3_in);
            [UB1 UB2] = deal(UBall.UB1_in,UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in,XBall.XB2_in);
            [RM10 RM20] = deal(RMall(1), RMall(2));

            Res_lims = [1e-1 2e4];
            Re_s_c_plot = 70;
            Res_CC_lims = [Res_lims(1) Re_s_c_plot];
            Res_TTV_lims = [Re_s_c_plot Res_lims(2)];
            Grat_lims = [8e-1 3e3];

            alpha_G_plot_R = 1.7;
            alpha_G_plot_S = 1.5;
            alpha_plot_R = alpha_G_plot_R-1;
            alpha_plot_S = alpha_G_plot_S-1;
            beta_plot_R = 1/(Re_s_c_plot^alpha_plot_R);
            beta_plot_S = 1/(Re_s_c_plot^alpha_plot_S);

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = plot(axs(1), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            plot(axs(1), Res_CC_lims, ones(size(Res_CC_lims)),'--','Color',[0 0 0], 'LineWidth',1.5,'DisplayName','$$ G/G_{cc} = 1 $$');
            fplot(axs(1), @(Re) beta_plot_R*(Re).^alpha_plot_R, Res_TTV_lims,':','Color',[0 0 0], 'LineWidth',2,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');
            fplot(axs(1), @(Re) beta_plot_S*(Re).^alpha_plot_S, Res_TTV_lims,'-.','Color',[0 0 0], 'LineWidth',1,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');

            legend_set_b(1) = plot(axs(2), NB1.Re_s, NB1.G_rat, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = plot(axs(2), NB2.Re_s, NB2.G_rat, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = plot(axs(2), NB3.Re_s, NB3.G_rat, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            legend_set_b(4) = plot(axs(2), RM10.Re_s, RM10.G_rat, RM10.specs,'Color', RM10.color, 'LineWidth', RM10.LW, 'MarkerSize', RM10.MS, 'DisplayName', RM10.label);
            legend_set_b(5) = plot(axs(2), RM20.Re_s, RM20.G_rat, RM20.specs,'Color', RM20.color, 'LineWidth', RM20.LW, 'MarkerSize', RM20.MS, 'DisplayName', RM20.label);
            plot(axs(2), Res_CC_lims, ones(size(Res_CC_lims)),'--','Color',[0 0 0], 'LineWidth',1.5,'DisplayName','$$ G/G_{cc} = 1 $$');
            fplot(axs(2), @(Re) beta_plot_R*(Re).^alpha_plot_R, Res_TTV_lims,':','Color',[0 0 0], 'LineWidth',2,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');
            fplot(axs(2), @(Re) beta_plot_S*(Re).^alpha_plot_S, Res_TTV_lims,'-.','Color',[0 0 0], 'LineWidth',1,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');

            legend_set_c(1) = plot(axs(3), UB1.Re_s, UB1.G_rat, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_c(2) = plot(axs(3), UB2.Re_s, UB2.G_rat, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
            plot(axs(3), Res_CC_lims, ones(size(Res_CC_lims)),'--','Color',[0 0 0], 'LineWidth',1.5,'DisplayName','$$ G/G_{cc} = 1 $$');
            fplot(axs(3), @(Re) beta_plot_R*(Re).^alpha_plot_R, Res_TTV_lims,':','Color',[0 0 0], 'LineWidth',2,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');
            fplot(axs(3), @(Re) beta_plot_S*(Re).^alpha_plot_S, Res_TTV_lims,'-.','Color',[0 0 0], 'LineWidth',1,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');

            % [Res_try_UB1 Grat_try_UB1] = UB1.comp_Grat_Res_KD(0.54);
            % [Res_try_UB2 Grat_try_UB2] = UB2.comp_Grat_Res_KD(0.54);
            % legend_set_c(3) = plot(axs(3), Res_try_UB1, Grat_try_UB1, 'p','Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            % legend_set_c(4) = plot(axs(3), Res_try_UB2, Grat_try_UB2, 'p','Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
            % plot(axs(3), Res_CC_lims, ones(size(Res_CC_lims)),'--','Color',[0 0 0], 'LineWidth',1.5,'DisplayName','$$ G/G_{cc} = 1 $$');

            legend_set_d(1) = plot(axs(4), XB1.Re_s, XB1.G_rat, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_d(2) = plot(axs(4), XB2.Re_s, XB2.G_rat, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);
            plot(axs(4), Res_CC_lims, ones(size(Res_CC_lims)),'--','Color',[0 0 0], 'LineWidth',1.5,'DisplayName','$$ G/G_{cc} = 1 $$');
            fplot(axs(4), @(Re) beta_plot_R*(Re).^alpha_plot_R, Res_TTV_lims,':','Color',[0 0 0], 'LineWidth',2,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');
            fplot(axs(4), @(Re) beta_plot_S*(Re).^alpha_plot_S, Res_TTV_lims,'-.','Color',[0 0 0], 'LineWidth',1,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');

            textbox_a = annotation('textbox', [0.09 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.56 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.09 0.38 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(c)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.56 0.38 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(d)}', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:2:3), '$$\mathcal{G}_{cc} = T_z/T_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:4), '$$Re_s$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'NorthEast', 'Interpreter', 'Latex');
            legend(axs(4), legend_set_d,'Location', 'NorthEast', 'Interpreter', 'Latex');

            axis(axs,[Res_lims Grat_lims])
            set(axs, 'XTick', [1e-1,1e0,1e1,1e2,1e3,1e4])
            set(axs, 'YTick', [1e0,1e1,1e2,1e3])

            fig_out = AYfig_;
        end
        function fig_out = PL_NUXB_T_vs_omega(obj, AYfig_, PFall, NBall, UBall, XBall)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            [PF1 PFR] = deal(PFall{1}, PFall{2});
            [NB1 NB2 NB3] = deal(NBall.NB1_in,NBall.NB2_in,NBall.NB3_in);
            [UB1 UB2] = deal(UBall.UB1_in,UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in,XBall.XB2_in);

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.omega, PF1.mu_torque, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = plot(axs(1), PFR.omega, PFR.mu_torque, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            legend_set_b(1) = plot(axs(2), NB1.omega, NB1.mu_torque, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = plot(axs(2), NB2.omega, NB2.mu_torque, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = plot(axs(2), NB3.omega, NB3.mu_torque, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            legend_set_c(1) = plot(axs(3), UB1.omega, UB1.mu_torque, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_c(2) = plot(axs(3), UB2.omega, UB2.mu_torque, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            legend_set_d(1) = plot(axs(4), XB1.omega, XB1.mu_torque, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_d(2) = plot(axs(4), XB2.omega, XB2.mu_torque, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            textbox_a = annotation('textbox', [0.09 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.56 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.09 0.38 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(c)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.56 0.38 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(d)}', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:2:3), '$$T_z$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:4), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(4), legend_set_d,'Location', 'SouthEast', 'Interpreter', 'Latex');

            axis(axs,obj.omega_T_range)
            set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2])
            set(axs, 'YTick', [1e-8,1e-6,1e-4,1e-2])

            fig_out = AYfig_;
        end
        function fig_out = UXB_FB_Grat_Res_compact(obj,AYfig_,UBall,XBall,FBall,FBext)
            axs_full = set_taustar_vs_Gamma_compact_axes(AYfig_);

            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);
            [FB1 FB2] = deal(FBall(1),FBall(2));

            UXB = {UB1; UB2; XB1; XB2};
            UXBlen = length(UXB);

            expALL = cell(4 + length(FB1.exp) + length(FB2.exp),1);
            iexp = 4;
            ileg = 0;
            for iFB = 1:length(FBall)
                FB = FBall(iFB);
                exp = FB.exp;
                for i = 1:length(exp)
                    expi = exp(i);
                    iexp = iexp + 1;
                    expALL{iexp} = expi;
                    if (expi.q < expi.qcrit_sgf)
                        [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                        [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                        [r_i r_o] = deal(fluid.r_i_def, fluid.r_o_def);
                        rc_full = expi.rc_Bingham;

                        % Re_try_full = (rc_full-r_i).*expi.Re_b_Bingham/(r_o-r_i);
                        Re_try_full = expi.Re_s_Bingham;

                        Re_try = Re_try_full(ifit);
                        Grat_full = expi.G_rat_Bingham;
                        Grat = Grat_full(ifit);

                        ileg = ileg +1;
                        legend_set_a(ileg) = plot(axs_full(1), Re_try_full, Grat_full, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        legend_set_b(ileg) = plot(axs_full(2), Re_try, Grat, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                    end
                end
            end

            for i = 1:length(FBext)
                expi = FBext{i};
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                [r_i r_o] = deal(fluid.r_i_def, fluid.r_o_def);
                rc_full = expi.rc_Bingham;

                % Re_try_full = (rc_full-r_i).*expi.Re_b_Bingham/(r_o-r_i);
                Re_try_full = expi.Re_s_Bingham;

                Re_try = Re_try_full(ifit);
                Grat_full = expi.G_rat_Bingham;
                Grat = Grat_full(ifit);

                ileg = ileg +1;
                legend_set_a(ileg) = plot(axs_full(1), Re_try_full, Grat_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Re_try, Grat, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end

            for i = 1:UXBlen
                expi = UXB{i};
                expALL{i} = expi;
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                [r_i r_o] = deal(fluid.r_i_def, fluid.r_o_def);
                rc_full = expi.rc_Bingham;

                % Re_try_full = (rc_full-r_i).*expi.Re_b_Bingham/(r_o-r_i);
                Re_try_full = expi.Re_s_Bingham;

                Re_try = Re_try_full(ifit);
                Grat_full = expi.G_rat_Bingham;
                Grat = Grat_full(ifit);

                ileg = ileg + 1;
                legend_set_a(ileg) = plot(axs_full(1), Re_try_full, Grat_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Re_try, Grat, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end
            set(axs_full,'XScale','log','YScale','log');
            [box_x box_y] = deal(axs_full(2).XLim, axs_full(2).YLim);
            box_a = [   box_x(1) box_y(1); ...
                        box_x(1) box_y(2); ...
                        box_x(2) box_y(2); ...
                        box_x(2) box_y(1)];

            patch(axs_full(1), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_a, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
            patch(axs_full(2), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_a, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);

            ylabel(axs_full(1), '$$ G_{rat} = G/G_{cc} $$ [dimensionless]', 'Interpreter', 'Latex')
            xlabel(axs_full(1), '$$ Re_s $$ [dimensionless]', 'Interpreter', 'Latex')

            ylabel(axs_full(2), '$$ G_{rat} $$', 'Interpreter', 'Latex')
            xlabel(axs_full(2), '$$ Re_s $$', 'Interpreter', 'Latex')

            legend_a = legend(axs_full(1), legend_set_a,'Location','NorthWest', 'Interpreter', 'Latex','NumColumns',1);

            % set(axs_full(2), 'XTick', [1e-3,1e-2,1e-1,1e-0])
            % ylim(axs_full(1), [1e-1 1e1])

            textbox_a = annotation('textbox', [0.15 0.125 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.425 0.775 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            fig_out = AYfig_;
        end
        function fig_out = UXB_FB_Grat_Rebstar_compact(obj,AYfig_,UBall,XBall,FBall,FBext)
            axs_full = set_taustar_vs_Gamma_compact_axes(AYfig_);

            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);
            [FB1 FB2] = deal(FBall(1),FBall(2));

            UXB = {UB1; UB2; XB1; XB2};
            UXBlen = length(UXB);

            subcrit_lims = [inf 0 inf 0];

            expALL = cell(4 + length(FB1.exp) + length(FB2.exp),1);
            iexp = 4;
            ileg = 0;

            for i = 1:length(FBext)
                expi = FBext{i};
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                [r_i r_o] = deal(fluid.r_i_def, fluid.r_o_def);
                rc_full = expi.rc_Bingham;

                Re_try_full = (rc_full-r_i).*expi.Re_b_Bingham/(r_o-r_i);
                % Re_try_full = expi.Re_b_Bingham;
                Re_try = Re_try_full(ifit);
                Grat_full = expi.G_rat_Bingham;
                Grat = Grat_full(ifit);

                tau_star_full = expi.tau_comp/expi.tau_y_Bingham;
                Gamma_full = (expi.mu_p_Bingham*expi.gamma_Bingham)/expi.tau_y_Bingham;
                tau_star = tau_star_full(ifit);
                Gamma = Gamma_full(ifit);

                subcrit_lims = update_lims(subcrit_lims, Re_try, Grat);

                ileg = ileg +1;
                legend_set_a(ileg) = plot(axs_full(1), Re_try_full, Grat_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end

            for iFB = 1:length(FBall)
                FB = FBall(iFB);
                exp = FB.exp;
                for i = 1:length(exp)
                    expi = exp(i);
                    iexp = iexp + 1;
                    expALL{iexp} = expi;
                    if (expi.q < expi.qcrit_sgf)
                        [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                        [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                        [r_i r_o] = deal(fluid.r_i_def, fluid.r_o_def);
                        rc_full = expi.rc_Bingham;

                        Re_try_full = (rc_full-r_i).*expi.Re_b_Bingham/(r_o-r_i);
                        % Re_try_full = expi.Re_b_Bingham;
                        Re_try = Re_try_full(ifit);
                        Grat_full = expi.G_rat_Bingham;
                        Grat = Grat_full(ifit);

                        tau_star_full = expi.tau/expi.tau_y_Bingham;
                        Gamma_full = (expi.mu_p_Bingham*expi.gamma_Bingham)/expi.tau_y_Bingham;
                        tau_star = tau_star_full(ifit);
                        Gamma = Gamma_full(ifit);

                        subcrit_lims = update_lims(subcrit_lims, Re_try, Grat);

                        ileg = ileg +1;
                        legend_set_a(ileg) = plot(axs_full(1), Re_try_full, Grat_full, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        legend_set_b(ileg) = plot(axs_full(2), Gamma, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                    end
                end
            end

            for i = 1:UXBlen
                expi = UXB{i};
                expALL{i} = expi;
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                [r_i r_o] = deal(fluid.r_i_def, fluid.r_o_def);
                rc_full = expi.rc_Bingham;

                Re_try_full = (rc_full-r_i).*expi.Re_b_Bingham/(r_o-r_i);
                % Re_try_full = expi.Re_b_Bingham;
                Re_try = Re_try_full(ifit);
                Grat_full = expi.G_rat_Bingham;
                Grat = Grat_full(ifit);

                tau_star_full = expi.tau_comp/expi.tau_y_Bingham;
                Gamma_full = (expi.mu_p_Bingham*expi.gamma_Bingham)/expi.tau_y_Bingham;
                tau_star = tau_star_full(ifit);
                Gamma = Gamma_full(ifit);

                subcrit_lims = update_lims(subcrit_lims, Re_try, Grat);

                ileg = ileg + 1;
                legend_set_a(ileg) = plot(axs_full(1), Re_try_full, Grat_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end
            set(axs_full,'XScale','log','YScale','log');
            [box_x_a box_y_a] = deal(subcrit_lims(1:2), subcrit_lims(3:4));
            [box_x_b box_y_b] = deal(axs_full(2).XLim, axs_full(2).YLim);
            box_a = [   box_x_a(1) box_y_a(1); ...
                        box_x_a(1) box_y_a(2); ...
                        box_x_a(2) box_y_a(2); ...
                        box_x_a(2) box_y_a(1)];
            box_b = [   box_x_b(1) box_y_b(1); ...
                        box_x_b(1) box_y_b(2); ...
                        box_x_b(2) box_y_b(2); ...
                        box_x_b(2) box_y_b(1)];

            Rebc = 50;
            alpha = 1.7;
            alpha_plot = alpha-1;
            beta_plot = Rebc^(-alpha_plot);

            legend_set_a(ileg+1) = fplot(axs_full(1), @(Re) beta_plot*Re.^alpha_plot, [Rebc 1e3],':','Color',[0 0 0], 'LineWidth',2.5,'DisplayName','$$ G_{rat} \propto Re^{\alpha-1}_{b} $$');
            legend_set_a(ileg+2) = fplot(axs_full(2), @(g) g+1, box_x_b,'--','Color',[0 0 0], 'LineWidth',1.5,'DisplayName','$$ \tau^* = \Gamma + 1 $$');

            patch(axs_full(1), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_a, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
            patch(axs_full(2), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_b, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);

            % ylabel(axs_full(1), '$$ G_{rat} = T_z / T_B(\omega_i) $$ [dimensionless]', 'Interpreter', 'Latex')
            ylabel(axs_full(1), '$$ \mathcal{G}_{B} = T_z / T_B(\omega_i) $$ [dimensionless]', 'Interpreter', 'Latex')
            xlabel(axs_full(1), '$$ Re_b^* = \rho_b \omega_i r_i (r_c - r_i) / \mu_p$$ [dimensionless]', 'Interpreter', 'Latex')

            ylabel(axs_full(2), '$$ \tau^* = \tau_i / \tau_y $$', 'Interpreter', 'Latex')
            xlabel(axs_full(2), '$$ \Gamma = \mu_p \cdot \dot{\gamma}_B (r_i) / \tau_y $$ [dimensionless]', 'Interpreter', 'Latex')

            legend_a = legend(axs_full(1), legend_set_a,'Location','NorthWest', 'Interpreter', 'Latex','NumColumns',1);

            axis(axs_full(2), [box_x_b box_y_b])
            axs_full(1).YLim(2) = 1e1;
            set(axs_full(2), 'XTick', [1e-3,1e-2,1e-1,1e-0])

            textbox_a = annotation('textbox', [0.150 0.125 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.425 0.775 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);

            fig_out = AYfig_;
        end
        function fig_out = UXB_FB_Tstar_vs_Tviscous_compact(obj,AYfig_,UBall,XBall,FBall,FBext)
            axs_full = set_taustar_vs_Gamma_compact_axes(AYfig_);

            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);
            [FB1 FB2] = deal(FBall(1),FBall(2));

            UXB = {UB1; UB2; XB1; XB2};
            UXBlen = length(UXB);

            expALL = cell(4 + length(FB1.exp) + length(FB2.exp),1);
            iexp = 4;
            ileg = 0;
            for iFB = 1:length(FBall)
                FB = FBall(iFB);
                exp = FB.exp;
                for i = 1:length(exp)
                    expi = exp(i);
                    iexp = iexp + 1;
                    expALL{iexp} = expi;

                    if (expi.q < expi.qcrit_sgf)
                        [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                        [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                        oi = expi.omega;
                        mup = expi.mu_p_Bingham;
                        r_i = fluid.r_i_def;
                        r_use = fluid.r_o_def;
                        L = fluid.h_def;
                        rho_rel = abs(expi.comp_rho_rel);
                        Tv = (mup*4*pi*L*r_i*r_i*(r_use.*r_use).*oi)./(r_use.*r_use - r_i*r_i);
                        Tqs_raw = pi*(r_use.*r_use)*rho_rel*9.81*L*L;
                        Tqs = (1-expi.q)*Tqs_raw;

                        T_star_full = expi.mu_torque./Tqs;
                        T_viscous_full = Tv./Tqs;
                        T_star = T_star_full(ifit);
                        T_viscous = T_viscous_full(ifit);

                        ileg = ileg +1;
                        legend_set_a(ileg) = plot(axs_full(1), T_viscous_full, T_star_full, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        legend_set_b(ileg) = plot(axs_full(2), T_viscous, T_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                    end
                end
            end

            for i = 1:length(FBext)
                expi = FBext{i};
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                oi = expi.omega;
                mup = expi.mu_p_Bingham;
                r_i = fluid.r_i_def;
                r_use = fluid.r_o_def;
                L = fluid.h_def;
                rho_rel = abs(expi.comp_rho_rel);
                Tv = (mup*4*pi*L*r_i*r_i*(r_use.*r_use).*oi)./(r_use.*r_use - r_i*r_i);
                Tqs_raw = pi*(r_use.*r_use)*rho_rel*9.81*L*L;
                Tqs = (1-expi.q)*Tqs_raw;

                T_star_full = expi.mu_torque./Tqs;
                T_viscous_full = Tv./Tqs;
                T_star = T_star_full(ifit);
                T_viscous = T_viscous_full(ifit);

                ileg = ileg +1;
                legend_set_a(ileg) = plot(axs_full(1), T_viscous_full, T_star_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), T_viscous, T_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end

            for i = 1:UXBlen
                expi = UXB{i};
                expALL{i} = expi;
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                oi = expi.omega;
                mup = expi.mu_p_Bingham;
                r_i = fluid.r_i_def;
                r_use = fluid.r_o_def;
                L = fluid.h_def;
                rho_rel = abs(expi.comp_rho_rel);
                Tv = (mup*4*pi*L*r_i*r_i*(r_use.*r_use).*oi)./(r_use.*r_use - r_i*r_i);
                Tqs_raw = pi*(r_use.*r_use)*rho_rel*9.81*L*L;
                Tqs = Tqs_raw;

                T_star_full = expi.mu_torque./Tqs;
                T_viscous_full = Tv./Tqs;
                T_star = T_star_full(ifit);
                T_viscous = T_viscous_full(ifit);

                ileg = ileg + 1;
                legend_set_a(ileg) = plot(axs_full(1), T_viscous_full, T_star_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), T_viscous, T_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end
            set(axs_full,'XScale','log','YScale','log');
            [box_x box_y] = deal(axs_full(2).XLim, axs_full(2).YLim);
            box_a = [   box_x(1) box_y(1); ...
                        box_x(1) box_y(2); ...
                        box_x(2) box_y(2); ...
                        box_x(2) box_y(1)];

            patch(axs_full(1), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_a, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
            patch(axs_full(2), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_a, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);

            ylabel(axs_full(1), '$$T^* = T_z / \pi r_o^2 \rho_{qs} g L^2 $$ [dimensionless]', 'Interpreter', 'Latex')
            xlabel(axs_full(1), '$$\hat{T} = T_v / \pi r_o^2 \rho_{qs} g L^2 $$ [dimensionless]', 'Interpreter', 'Latex')

            ylabel(axs_full(2), '$$T^*$$', 'Interpreter', 'Latex')
            xlabel(axs_full(2), '$$\hat{T}$$', 'Interpreter', 'Latex')

            legend_a = legend(axs_full(1), legend_set_a,'Location','NorthWest', 'Interpreter', 'Latex','NumColumns',1);

            % set(axs_full(2), 'XTick', [1e-3,1e-2,1e-1,1e-0])
            ylim(axs_full(1), [1e-1 1e1])

            textbox_a = annotation('textbox', [0.15 0.125 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.425 0.775 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            fig_out = AYfig_;
        end
        function fig_out = UXB_FB_taustar_vs_Gamma_compact(obj,AYfig_,UBall,XBall,FBall,FBext)
            axs_full = set_taustar_vs_Gamma_compact_axes(AYfig_);

            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);
            [FB1 FB2] = deal(FBall(1),FBall(2));

            UXB = {UB1; UB2; XB1; XB2};
            UXBlen = length(UXB);

            expALL = cell(4 + length(FB1.exp) + length(FB2.exp),1);
            iexp = 4;
            ileg = 0;
            for iFB = 1:length(FBall)
                FB = FBall(iFB);
                exp = FB.exp;
                for i = 1:length(exp)
                    expi = exp(i);
                    iexp = iexp + 1;
                    expALL{iexp} = expi;

                    if (expi.q < expi.qcrit_sgf)
                        [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                        [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));
                        tau_star = tfit/expi.tau_y_Bingham;
                        Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                        tau_star_full = expi.tau/expi.tau_y_Bingham;
                        Gamma_full = (expi.mu_p_Bingham*expi.gamma_Bingham)/expi.tau_y_Bingham;

                        ileg = ileg +1;
                        legend_set_a(ileg) = plot(axs_full(1), Gamma_full, tau_star_full, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        legend_set_b(ileg) = plot(axs_full(2), Gamma, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                    end
                end
            end

            for i = 1:length(FBext)
                expi = FBext{i};
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                tau_star = tfit/expi.tau_y_Bingham;
                Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                tau_star_full = expi.tau_comp/expi.tau_y_Bingham;
                Gamma_full = (expi.mu_p_Bingham*expi.gamma_Bingham)/expi.tau_y_Bingham;

                ileg = ileg +1;
                legend_set_a(ileg) = plot(axs_full(1), Gamma_full, tau_star_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end

            for i = 1:UXBlen
                expi = UXB{i};
                expALL{i} = expi;
                [mup ty] = deal(expi.mu_p_Bingham, expi.tau_y_Bingham);
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                tau_star = tfit/expi.tau_y_Bingham;
                Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                tau_star_full = expi.tau_comp/expi.tau_y_Bingham;
                Gamma_full = (expi.mu_p_Bingham*expi.gamma_Bingham)/expi.tau_y_Bingham;


                ileg = ileg + 1;
                legend_set_a(ileg) = plot(axs_full(1), Gamma_full, tau_star_full, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end
            set(axs_full,'XScale','log','YScale','log');
            [box_x box_y] = deal(axs_full(2).XLim, axs_full(2).YLim);
            box_a = [   box_x(1) box_y(1); ...
                        box_x(1) box_y(2); ...
                        box_x(2) box_y(2); ...
                        box_x(2) box_y(1)];

            patch(axs_full(1), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_a, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
            patch(axs_full(2), 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_a, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);

            ylabel(axs_full(1), '$$\tau^* = \tau/\tau_y$$ [dimensionless]', 'Interpreter', 'Latex')
            xlabel(axs_full(1), '$$\Gamma = \dot{\gamma}_B \mu_p/\tau_y$$ [dimensionless]', 'Interpreter', 'Latex')

            ylabel(axs_full(2), '$$\tau^*$$', 'Interpreter', 'Latex')
            xlabel(axs_full(2), '$$\Gamma$$', 'Interpreter', 'Latex')

            legend_a = legend(axs_full(1), legend_set_a,'Location','NorthWest', 'Interpreter', 'Latex','NumColumns',1);

            set(axs_full(2), 'XTick', [1e-3,1e-2,1e-1,1e-0])

            textbox_a = annotation('textbox', [0.15 0.125 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.425 0.775 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            fig_out = AYfig_;
        end
        function fig_out = FB_Bingham_fits_Grat_vs_Reb(obj, AYfig_, FB1, FB2, PFR)
            axs = set_FB_Bingham_Grat_vs_Reb_axes(AYfig_);

            for i = 1:length(FB1.exp)
                legend_set_a(i) = plot(axs(1), FB1.exp(i).Re_b_Bingham, FB1.exp(i).G_rat_Bingham, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_b(i) = plot(axs(2), FB2.exp(i).Re_b_Bingham, FB2.exp(i).G_rat_Bingham, FB2.exp(i).specs,'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            for i = 1:length(FB1.exp)
                legend_set_c(i) = plot(axs(3), FB1.exp(i).q, (FB1.exp(i).tau_y_Bingham)/(FB1.exp(i).tau_static), FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_d(i) = plot(axs(4), FB2.exp(i).q, (FB2.exp(i).tau_y_Bingham)/(FB2.exp(i).tau_static), FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            for i = 1:length(FB1.exp)
                legend_set_e(i) = plot(axs(5), FB1.exp(i).q, FB1.exp(i).mu_p_Bingham, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_f(i) = plot(axs(6), FB2.exp(i).q, FB2.exp(i).mu_p_Bingham, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            if (nargin==5)
                legend_set_a(length(FB1.exp) + 1) = plot(axs(1), PFR.Re_b_comp, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
                legend_set_B(length(FB2.exp) + 1) = plot(axs(2), PFR.Re_b_comp, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            end
            % textbox_a = annotation('textbox', obj.textbox_pos22_a_NE, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_b = annotation('textbox', obj.textbox_pos22_b_NE, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_c = annotation('textbox', obj.textbox_pos22_c_NE, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            % textbox_d = annotation('textbox', obj.textbox_pos22_d_NE, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs(1:2),'YScale', 'log','XScale','log');

            % set(axs(1:2), 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1])

            ylabel(axs(1), '$$G_{rat} = G/G_B$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(1:2), '$$Re_b$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            ylabel(axs(3), '$$\tau_y/\tau_{q = 0}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(5), '$$\tilde{\mu}_{p}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(5:6), '$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            % legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            % legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs(1:2),[1e-2 1e4 1e-2 1e1])
            % axis(axs(3),[0 2 1e-5 1.0])
            % axis(axs(4),[0 16 1e-5 1.0])
            % axis(axs(5),[0 2 1e-1 1.2])
            % axis(axs(6),[0 16 1e-2 1.2])

            fig_out = AYfig_;
        end
        function fig_out = ALL_G_vs_Reb(obj, AYfig_, FBall, PFall, NBall, UBall, XBall, FBext)
            AYfig_.init_tiles([3,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            [FB1 FB2] = deal(FBall(1), FBall(2));
            [PF1 PFR] = deal(PFall{1}, PFall{2});
            [NB1 NB2 NB3] = deal(NBall.NB1_in,NBall.NB2_in,NBall.NB3_in);
            [UB1 UB2] = deal(UBall.UB1_in,UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in,XBall.XB2_in);

            G_lims = [1e-1 4e6];
            Reb_lims = [1e-2 1e4];
            Reb_c_plot = 50;
            Reb_CC_lims = [Reb_lims(1) Reb_c_plot];
            Reb_TTV_lims = [Reb_c_plot Reb_lims(2)];

            alpha_plot_R = 1.7;
            alpha_plot_S = 1.5;
            alpha_plot_try = 2.0;
            beta_plot_R = (obj.G_obs_Reb_slope*Reb_c_plot)/(Reb_c_plot^alpha_plot_R);
            beta_plot_S = (obj.G_obs_Reb_slope*Reb_c_plot)/(Reb_c_plot^alpha_plot_S);
            beta_plot_try = (obj.G_obs_Reb_slope*Reb_c_plot)/(Reb_c_plot^alpha_plot_try);

            for i = 1:6
                fplot(axs(i), @(Re) obj.G_obs_Reb_slope*(Re), Reb_CC_lims,'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{4 \pi r_i r_o^2}{(r_o + r_i) \cdot (r_o-r_i)^2} Re_b $$')
                fplot(axs(i), @(Re) beta_plot_R*(Re).^alpha_plot_R, Reb_TTV_lims,':','Color',[0 0 0], 'LineWidth',2,'DisplayName','$$ G/G_{cc} \propto Re_b^{0.7} $$');
                fplot(axs(i), @(Re) beta_plot_S*(Re).^alpha_plot_S, Reb_TTV_lims,'-.','Color',[0 0 0], 'LineWidth',1,'DisplayName','$$ G/G_{cc} \propto Re_b^{0.7} $$');
                % fplot(axs(i), @(Re) beta_plot_try*(Re).^alpha_plot_try, Reb_TTV_lims,'-','Color',[0 0 0], 'LineWidth',1,'DisplayName','$$ G/G_{cc} \propto Re_b^{0.7} $$');
            end

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.Re_b_comp, PF1.G, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = plot(axs(1), PFR.Re_b_comp, PFR.G, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            legend_set_b(1) = plot(axs(2), NB1.Re_b_comp, NB1.G, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = plot(axs(2), NB2.Re_b_comp, NB2.G, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = plot(axs(2), NB3.Re_b_comp, NB3.G, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            legend_set_c(1) = plot(axs(3), UB1.Re_b_comp, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_c(2) = plot(axs(3), UB2.Re_b_comp, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            legend_set_d(1) = plot(axs(4), XB1.Re_b_comp, XB1.G, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_d(2) = plot(axs(4), XB2.Re_b_comp, XB2.G, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            if (nargin==8)
                legend_set_e(1) = plot(axs(5), FBext{1}.Re_b_Bingham, FBext{1}.G_b_Bingham, FBext{1}.specs, 'Color', FBext{1}.color, 'LineWidth', FBext{1}.LW, 'MarkerSize', FBext{1}.MS, 'DisplayName', FBext{1}.label);
                legend_set_e(2) = plot(axs(5), FBext{2}.Re_b_Bingham, FBext{2}.G_b_Bingham, FBext{2}.specs, 'Color', FBext{2}.color, 'LineWidth', FBext{2}.LW, 'MarkerSize', FBext{2}.MS, 'DisplayName', FBext{2}.label);
                legend_set_e(3) = plot(axs(5), FBext{3}.Re_b_Bingham, FBext{3}.G_b_Bingham, FBext{3}.specs, 'Color', FBext{3}.color, 'LineWidth', FBext{3}.LW, 'MarkerSize', FBext{3}.MS, 'DisplayName', FBext{3}.label);
                legend_set_e(4) = plot(axs(5), FBext{4}.Re_b_Bingham, FBext{4}.G_b_Bingham, FBext{4}.specs, 'Color', FBext{4}.color, 'LineWidth', FBext{4}.LW, 'MarkerSize', FBext{4}.MS, 'DisplayName', FBext{4}.label);
                ie_start = 4;
            else
                ie_start=0;
            end

            for i=1:length(FB1.exp)
              legend_set_e(ie_start + i) = plot(axs(5), FB1.exp(i).get_Re_b, FB1.exp(i).get_G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end

            for i=1:length(FB2.exp)
              legend_set_f(i) = plot(axs(6), FB2.exp(i).get_Re_b, FB2.exp(i).get_G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', [0.09 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.56 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.09 0.60 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(c)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.56 0.60 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(d)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_e = annotation('textbox', [0.09 0.33 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(e)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_f = annotation('textbox', [0.56 0.33 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(f)}', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:2:5), '$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(5:6), '$$Re_b = \rho \omega_i r_i (r_o-r_i) / \mu$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2,1e3,1e4])
            set(axs, 'YTick', [1e0,1e1,1e2,1e3,1e4,1e5,1e6])

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(4), legend_set_d,'Location', 'SouthEast', 'Interpreter', 'Latex');

            leg_comb = legend(axs(5), [legend_set_e legend_set_f],'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 4);
            leg_comb.Layout.Tile = 'south';

            axis(axs,[Reb_lims G_lims])

            fig_out = AYfig_;
        end
        function fig_out = ALL_G_vs_Reb_tall(obj, AYfig_, FBall, PFall, NBall, UBall, XBall, FBext)
            AYfig_.init_tiles([3,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            [FB1 FB2] = deal(FBall(1), FBall(2));
            [PF1 PFR] = deal(PFall{1}, PFall{2});
            [NB1 NB2 NB3] = deal(NBall.NB1_in,NBall.NB2_in,NBall.NB3_in);
            [UB1 UB2] = deal(UBall.UB1_in,UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in,XBall.XB2_in);

            G_lims = [3e-4 4e6];
            Reb_lims = [1e-2 1e4];
            Reb_c_plot = 50;
            Reb_CC_lims = [Reb_lims(1) Reb_c_plot];
            Reb_TTV_lims = [Reb_c_plot Reb_lims(2)];

            alpha_plot_R = 1.7;
            alpha_plot_S = 1.5;
            beta_plot_R = (obj.G_obs_Reb_slope*Reb_c_plot)/(Reb_c_plot^alpha_plot_R);
            beta_plot_S = (obj.G_obs_Reb_slope*Reb_c_plot)/(Reb_c_plot^alpha_plot_S);

            for i = 1:6
                fplot(axs(i), @(Re) obj.G_obs_Reb_slope*(Re), Reb_CC_lims,'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', '$$ \frac{4 \pi r_i r_o^2}{(r_o + r_i) \cdot (r_o-r_i)^2} Re_b $$')
                fplot(axs(i), @(Re) beta_plot_R*(Re).^alpha_plot_R, Reb_TTV_lims,':','Color',[0 0 0], 'LineWidth',2,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');
                fplot(axs(i), @(Re) beta_plot_S*(Re).^alpha_plot_S, Reb_TTV_lims,'-.','Color',[0 0 0], 'LineWidth',1,'DisplayName','$$ G/G_{cc} \propto Re_s^{0.7} $$');
            end

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.Re_b_comp, PF1.G, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = plot(axs(1), PFR.Re_b_comp, PFR.G, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            legend_set_b(1) = plot(axs(2), NB1.Re_b_comp, NB1.G, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = plot(axs(2), NB2.Re_b_comp, NB2.G, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = plot(axs(2), NB3.Re_b_comp, NB3.G, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            legend_set_c(1) = plot(axs(3), UB1.Re_b_comp, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_c(2) = plot(axs(3), UB2.Re_b_comp, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            legend_set_d(1) = plot(axs(4), XB1.Re_b_comp, XB1.G, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_d(2) = plot(axs(4), XB2.Re_b_comp, XB2.G, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            if (nargin==8)
                legend_set_e(1) = plot(axs(5), FBext{1}.Re_b_Bingham, FBext{1}.G_b_Bingham, FBext{1}.specs, 'Color', FBext{1}.color, 'LineWidth', FBext{1}.LW, 'MarkerSize', FBext{1}.MS, 'DisplayName', FBext{1}.label);
                legend_set_e(2) = plot(axs(5), FBext{2}.Re_b_Bingham, FBext{2}.G_b_Bingham, FBext{2}.specs, 'Color', FBext{2}.color, 'LineWidth', FBext{2}.LW, 'MarkerSize', FBext{2}.MS, 'DisplayName', FBext{2}.label);
                legend_set_e(3) = plot(axs(5), FBext{3}.Re_b_Bingham, FBext{3}.G_b_Bingham, FBext{3}.specs, 'Color', FBext{3}.color, 'LineWidth', FBext{3}.LW, 'MarkerSize', FBext{3}.MS, 'DisplayName', FBext{3}.label);
                legend_set_e(4) = plot(axs(5), FBext{4}.Re_b_Bingham, FBext{4}.G_b_Bingham, FBext{4}.specs, 'Color', FBext{4}.color, 'LineWidth', FBext{4}.LW, 'MarkerSize', FBext{4}.MS, 'DisplayName', FBext{4}.label);
                ie_start = 4;
            else
                ie_start=0;
            end

            for i=1:length(FB1.exp)
              legend_set_e(ie_start + i) = plot(axs(5), FB1.exp(i).get_Re_b, FB1.exp(i).get_G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end

            for i=1:length(FB2.exp)
              legend_set_f(i) = plot(axs(6), FB2.exp(i).get_Re_b, FB2.exp(i).get_G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', [0.09 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.56 0.86 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.09 0.54 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.56 0.54 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_e = annotation('textbox', [0.09 0.21 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'e)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_f = annotation('textbox', [0.56 0.21 0.1 0.1], 'Interpreter', 'LaTeX', 'String', 'f)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:2:5), '$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(5:6), '$$Re_b = \rho \omega_i r_i (r_o-r_i) / \mu$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2,1e3,1e4])

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(4), legend_set_d,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(5), legend_set_e,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(6), legend_set_f,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs,[Reb_lims G_lims])

            fig_out = AYfig_;
        end
        function fig_out = UXB_Grat_vs_Res_NB_comparison(obj, AYfig_, UBall, XBall, NBall)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            [UB1 UB2] = deal(UBall.UB1_in,UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in,XBall.XB2_in);
            [NB1 NB2 NB3] = deal(NBall.NB1_in,NBall.NB2_in,NBall.NB3_in);

            [beta_G, alpha_G] = deal(NBall.powerfit.b, NBall.powerfit.m);
            [beta_Gr, alpha_Gr] = deal(NBall.powerfit_Grat_Res.b, NBall.powerfit_Grat_Res.m);

            Res_lims = [1e2 2e3];

            Gcomp = @(Res,G) (G-(beta_G*Res.^(alpha_G)))./(beta_G*Res.^(alpha_G));
            Grcomp = @(Res,Grat) Grat./(beta_Gr*Res.^(alpha_Gr));

            legend_set_a(1) = plot(axs(1), UB1.Re_s, Gcomp(UB1.Re_s,UB1.G), UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'Displayname', UB1.label);
            legend_set_a(2) = plot(axs(1), UB2.Re_s, Gcomp(UB2.Re_s,UB2.G), UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'Displayname', UB2.label);
            % legend_set_a(1) = plot(axs(1), UB1.Re_s, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'Displayname', UB1.label);
            % legend_set_a(2) = plot(axs(1), UB2.Re_s, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'Displayname', UB2.label);
            % legend_set_a(3) = fplot(axs(1), @(Re) beta_G*(Re).^alpha_G, Res_lims,'-','Color',NB3.color, 'LineWidth',NB3.LW);


            legend_set_b(1) = plot(axs(2), UB1.Re_s, Grcomp(UB1.Re_s,UB1.G_rat), UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'Displayname', UB1.label);
            legend_set_b(2) = plot(axs(2), UB2.Re_s, Grcomp(UB2.Re_s,UB2.G_rat), UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'Displayname', UB2.label);

            textbox_a = annotation('textbox', obj.textbox_pos2_a_SW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_SW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            % ylabel(axs(1), '$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            % ylabel(axs(2), '$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            % legend(axs(1),legend_set_a,'Location', 'Northeast', 'Interpreter', 'Latex');
            % legend(axs(2), legend_set_b,'Location', 'Northeast', 'Interpreter', 'Latex');

            % set(axs(1), 'XScale', 'log', 'YScale', 'log');
            set(axs(1), 'XScale', 'log', 'YScale', 'linear');
            set(axs(2), 'XScale', 'log', 'YScale', 'linear');

            xlim(axs, Res_lims)

            fig_out = AYfig_;
        end
        function fig_out = UXB_NB_PL_Grat_vs_Res(obj, AYfig_, UB1, UB2, XB1, XB2, PF1, PFR, NB1, NB2, NB3)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            [bGratPF1 aGratPF1 bGratPFR aGratPFR] = deal(PF1.powerfit_Grat_Res.b,PF1.powerfit_Grat_Res.m,PFR.powerfit_Grat_Res.b,PFR.powerfit_Grat_Res.m);

            legend_set_a(1) = plot(axs(1), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'Displayname', PF1.label);
            legend_set_a(2) = plot(axs(1), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'Displayname', PFR.label);
            legend_set_a(3) = plot(axs(1), NB1.Re_s, NB1.G_rat, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'Displayname', NB1.label);
            legend_set_a(4) = plot(axs(1), NB2.Re_s, NB2.G_rat, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'Displayname', NB2.label);
            legend_set_a(5) = plot(axs(1), NB3.Re_s, NB3.G_rat, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'Displayname', NB3.label);
            legend_set_a(6) = plot(axs(1), UB1.Re_s, UB1.G_rat, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'Displayname', UB1.label);
            legend_set_a(7) = plot(axs(1), UB2.Re_s, UB2.G_rat, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'Displayname', UB2.label);

            legend_set_b(1) = plot(axs(2), PF1.Re_s, PF1.G_rat, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'Displayname', PF1.label);
            legend_set_b(2) = plot(axs(2), PFR.Re_s, PFR.G_rat, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'Displayname', PFR.label);
            legend_set_b(3) = plot(axs(2), NB1.Re_s, NB1.G_rat, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'Displayname', NB1.label);
            legend_set_b(4) = plot(axs(2), NB2.Re_s, NB2.G_rat, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'Displayname', NB2.label);
            legend_set_b(5) = plot(axs(2), NB3.Re_s, NB3.G_rat, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'Displayname', NB3.label);
            legend_set_b(6) = plot(axs(2), XB1.Re_s, XB1.G_rat, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'Displayname', XB1.label);
            legend_set_b(7) = plot(axs(2), XB2.Re_s, XB2.G_rat, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'Displayname', XB2.label);

            textbox_a = annotation('textbox', obj.textbox_pos2_a_SW,   'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_SW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            ylabel(axs(1), '$$G_{rat} = G/G_{cc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'Northeast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'Northeast', 'Interpreter', 'Latex');

            set(axs, 'XScale', 'log', 'YScale', 'log');

            axis(axs,[1e-1 2e4 8e-1 3e3])

            fig_out = AYfig_;
        end
        function fig_out = FB_phi_vs_omegai_vs_q(obj,AYfig_,FB1_phi_exp)
            AYfig_.init_tiles([1,1]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            poq = FB1_phi_exp.phi_o_q_mat;
            [phi_vec o_vec q_vec] = deal(poq(:,1), poq(:,2), poq(:,3));

            scatter3(axs(1),q_vec,o_vec,phi_vec,'*','LineWidth',0.75,'SizeData',50,'CData',phi_vec);

            % set(axs(1),'YScale','log');
            % ylim(axs(1), [1e-3, 1.1e2]);

            zlabel('$$\phi$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel('$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            set(axs(1),'YScale','log','YTick',[1e-2 1e-1 1e0 1e1 1e2])
            cbar = colorbar(axs(1));
            cbar.Label.String = '$$ \phi $$';
            cbar.Label.Interpreter = 'Latex';
            cbar.Label.FontSize = 16;
            colormap(axs(1),cool);
            view([45 45 45]);
            fig_out = AYfig_;
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

            ylabel(axs(1), '$$T_{z}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');

            set(axs, 'XScale', 'log', 'YScale', 'log');

            % axis(axs,obj.omega_tau_range)
            axis(axs,obj.omega_T_range)
            set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2])

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

            ylabel(axs(1), '$$T_{z}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex');

            set(axs, 'XScale', 'log', 'YScale', 'log');

            % axis(axs,obj.omega_tau_range)
            axis(axs,obj.omega_T_range)
            set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2])


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

            textbox_a = annotation('textbox', obj.textbox_pos2_a_SW, 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_SW, 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1), '$$T_{z}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs,obj.omega_T_range)
            set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2])
            set(axs, 'YTick', [1e-8,1e-6,1e-4,1e-2])

            fig_out = AYfig_;
        end
        function fig_out = FB_appmu_vs_omegai(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i = 1:length(FB1.exp)
                appmup = glass_particles.compute_appmu_Bingham(FB1.exp(i).tau_y_Bingham, FB1.exp(i).omega, FB1.exp(i).tau);
                legend_set_a(i) = plot(axs(1), FB1.exp(i).omega, appmup, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                appmup = glass_particles.compute_appmu_Bingham(FB2.exp(i).tau_y_Bingham, FB2.exp(i).omega, FB2.exp(i).tau);
                legend_set_b(i) = plot(axs(2), FB2.exp(i).omega, appmup, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', obj.textbox_pos2_a_NW, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos2_b_NW, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1), '$$\mu_{app}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            % axis(axs,obj.omega_appmu_range)

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
        function fig_out = FB_tauyrat_vs_q(obj, AYfig_, FB1, FB2, EC000, EC050, EC075, EC100)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            % legend_set_a(1) = plot(axs(1), EC000.q, EC000.tau_y_rat, EC000.specs, 'Color', EC000.color, 'LineWidth', EC000.LW_L, 'MarkerSize', EC000.MS_L, 'DisplayName', EC000.label);
            % legend_set_a(2) = plot(axs(1), EC050.q, EC050.tau_y_rat, EC050.specs, 'Color', EC050.color, 'LineWidth', EC050.LW_L, 'MarkerSize', EC050.MS_L, 'DisplayName', EC050.label);
            % legend_set_a(3) = plot(axs(1), EC075.q, EC075.tau_y_rat, EC075.specs, 'Color', EC075.color, 'LineWidth', EC075.LW_L, 'MarkerSize', EC075.MS_L, 'DisplayName', EC075.label);
            % legend_set_a(4) = plot(axs(1), EC100.q, EC100.tau_y_rat, EC100.specs, 'Color', EC100.color, 'LineWidth', EC100.LW_L, 'MarkerSize', EC100.MS_L, 'DisplayName', EC100.label);
            for i = 1:length(FB1.exp)
                legend_set_a(i) = plot(axs(1), FB1.exp(i).q, (FB1.exp(i).tau_y_Bingham)/(FB1.exp(i).tau_static), FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_b(i) = plot(axs(2), FB2.exp(i).q, (FB2.exp(i).tau_y_Bingham)/(FB2.exp(i).tau_static), FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            for i = 1:length(FB1.exp)
                legend_set_c(i) = plot(axs(3), FB1.exp(i).q, FB1.exp(i).mu_p_Bingham, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_d(i) = plot(axs(4), FB2.exp(i).q, FB2.exp(i).mu_p_Bingham, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
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
            axis(axs(3),[0 2 1e-1 1.2])
            axis(axs(4),[0 16 1e-2 1.2])

            fig_out = AYfig_;
        end
        function fig_out = FB_tauyrat_mup_vs_q(obj, AYfig_, FB1, FB2, EC000, EC050, EC075, EC100)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            qlims_FB1 = [0 2.05];
            qlims_FB2 = [0 15.2];

            tyr_lims = [1e-4 1e0];
            mup_lims = [1e-2 3e-1];

            % legend_set_a(1) = plot(axs(1), EC000.q, EC000.tau_y_rat, EC000.specs, 'Color', EC000.color, 'LineWidth', EC000.LW_L, 'MarkerSize', EC000.MS_L, 'DisplayName', EC000.label);
            % legend_set_a(2) = plot(axs(1), EC050.q, EC050.tau_y_rat, EC050.specs, 'Color', EC050.color, 'LineWidth', EC050.LW_L, 'MarkerSize', EC050.MS_L, 'DisplayName', EC050.label);
            % legend_set_a(3) = plot(axs(1), EC075.q, EC075.tau_y_rat, EC075.specs, 'Color', EC075.color, 'LineWidth', EC075.LW_L, 'MarkerSize', EC075.MS_L, 'DisplayName', EC075.label);
            % legend_set_a(4) = plot(axs(1), EC100.q, EC100.tau_y_rat, EC100.specs, 'Color', EC100.color, 'LineWidth', EC100.LW_L, 'MarkerSize', EC100.MS_L, 'DisplayName', EC100.label);
            for ax_active = [axs(1) axs(2)]
                for i = 1:length(FB1.exp)
                    plot(ax_active, FB1.exp(i).q, (FB1.exp(i).tau_y_Bingham)/(FB1.exp(i).tau_static), FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
                end
                for i = 1:length(FB2.exp)
                    plot(ax_active, FB2.exp(i).q, (FB2.exp(i).tau_y_Bingham)/(FB2.exp(i).tau_static), FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
                end
            end

            for ax_active = [axs(3) axs(4)]
                for i = 1:length(FB1.exp)
                    plot(ax_active, FB1.exp(i).q, FB1.exp(i).mu_p_Bingham, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
                end
                for i = 1:length(FB2.exp)
                    plot(ax_active, FB2.exp(i).q, FB2.exp(i).mu_p_Bingham, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
                end
            end

            textbox_a = annotation('textbox', obj.textbox_pos22_a_NE, 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', obj.textbox_pos22_b_NE, 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', obj.textbox_pos22_c_NE, 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', obj.textbox_pos22_d_NE, 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs(1:2),'YScale', 'log');
            % set(axs,'YScale', 'log');

            set(axs(1:2), 'YTick', [1e-4,1e-3,1e-2,1e-1,1e0])

            ylabel(axs(1), '$$\tau_y/\tau_{q = 0}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3), '$$\mu_p$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:4), '$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            % legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            % legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            ylim([axs(1) axs(2)], tyr_lims);
            ylim([axs(3) axs(4)], mup_lims);
            xlim([axs(1) axs(3)], qlims_FB1);
            xlim([axs(2) axs(4)], qlims_FB2);

            fig_out = AYfig_;
        end
        function fig_out = FB_tauyrat_mup_vs_q_compact(obj, AYfig_, FB1, FB2, FBext)
            axs = set_mup_tyr_vs_q_compact_axes(AYfig_);
            axs_FB1_qlim = [axs(1) axs(3)];
            axs_FB2_qlim = [axs(2) axs(4)];
            axs_tyr = [axs(1) axs(2)];
            axs_mup = [axs(3) axs(4)];

            qlims_FB1 = [0 2.05];
            qlims_FB2 = [0 15.2];

            tyr_lims = [1e-4 1e0];
            % mup_lims = [1e-2 3e-1];
            mup_lims = [0 3e-1];

            box_tyr_vs_q = [    qlims_FB1(1) tyr_lims(1); ...
                                qlims_FB1(1) tyr_lims(2); ...
                                qlims_FB1(2) tyr_lims(2); ...
                                qlims_FB1(2) tyr_lims(1)];

            box_mup_vs_q = [    qlims_FB1(1) mup_lims(1); ...
                                qlims_FB1(1) mup_lims(2); ...
                                qlims_FB1(2) mup_lims(2); ...
                                qlims_FB1(2) mup_lims(1)];

            for ax_active = axs_tyr
                patch(ax_active, 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_tyr_vs_q, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
                if (nargin>4)
                    for i = 1:length(FBext)
                        expi = FBext{i};
                        plot(ax_active, expi.q, (expi.tau_y_Bingham)/(expi.tau_static), expi.specs,'Color', expi.color, 'LineWidth', expi.LW_L, 'MarkerSize', expi.MS_L, 'DisplayName', expi.label);
                    end
                end
                for i = 1:length(FB1.exp)
                    plot(ax_active, FB1.exp(i).q, (FB1.exp(i).tau_y_Bingham)/(FB1.exp(i).tau_static), FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
                end
                for i = 1:length(FB2.exp)
                    plot(ax_active, FB2.exp(i).q, (FB2.exp(i).tau_y_Bingham)/(FB2.exp(i).tau_static), FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
                end
            end

            for ax_active = axs_mup
                patch(ax_active , 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_mup_vs_q, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
                if (nargin>4)
                    for i = 1:length(FBext)
                        expi = FBext{i};
                        plot(ax_active, expi.q, expi.mu_p_Bingham, expi.specs,'Color', expi.color, 'LineWidth', expi.LW_L, 'MarkerSize', expi.MS_L, 'DisplayName', expi.label);
                    end
                end
                for i = 1:length(FB1.exp)
                    plot(ax_active, FB1.exp(i).q, FB1.exp(i).mu_p_Bingham, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
                end
                for i = 1:length(FB2.exp)
                    plot(ax_active, FB2.exp(i).q, FB2.exp(i).mu_p_Bingham, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
                end
            end

            fplot(axs_tyr(1), @(q) 1-q, [0 0.9], '--', 'Color', [0 0 0], 'LineWidth', 1, 'DisplayName', '$$ 1 - q $$');
            fplot(axs_tyr(2), @(q) 1-q, [0 0.9], '--', 'Color', [0 0 0], 'LineWidth', 1, 'DisplayName', '$$ 1 - q $$');

            textbox_a = annotation('textbox', [0.425 0.20 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.25 0.7 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.925 0.20 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(c)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.75 0.85 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(d)}', 'LineStyle', 'none', 'FontSize', 16);

            set(axs_tyr,'YScale', 'log');
            set(axs_tyr(1), 'YTick', [1e-4,1e-2,1e0])
            set(axs_tyr(2), 'YTick', [1e-4,1e-3,1e-2,1e-1,1e0])

            ylabel(axs_FB2_qlim(1), '$$\tau_y/\tau_{q = 0}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs_FB2_qlim(2), '$$\mu_p$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_FB2_qlim, '$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            ylabel(axs_FB1_qlim(1), '$$\tau_y/\tau_{q = 0}$$', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs_FB1_qlim(2), '$$\mu_p$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_FB1_qlim, '$$q = Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',12)

            ylim(axs_tyr, tyr_lims);
            ylim(axs_mup, mup_lims);
            xlim(axs_FB1_qlim, qlims_FB1);
            xlim(axs_FB2_qlim, qlims_FB2);

            fig_out = AYfig_;
        end
        function fig_out = FB_Carreau_par_vs_q(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([3,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i = 1:length(FB1.exp)
                legend_set_a(i) = plot(axs(1), FB1.exp(i).q, FB1.exp(i).mu0_Carreau, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_b(i) = plot(axs(2), FB2.exp(i).q, FB2.exp(i).mu0_Carreau, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            for i = 1:length(FB1.exp)
                legend_set_c(i) = plot(axs(3), FB1.exp(i).q, FB1.exp(i).lambda_Carreau, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_d(i) = plot(axs(4), FB2.exp(i).q, FB2.exp(i).lambda_Carreau, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            for i = 1:length(FB1.exp)
                legend_set_c(i) = plot(axs(5), FB1.exp(i).q, FB1.exp(i).n_Carreau, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end

            for i = 1:length(FB2.exp)
                legend_set_d(i) = plot(axs(6), FB2.exp(i).q, FB2.exp(i).n_Carreau, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            textbox_a = annotation('textbox', [0.1, 0.675, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'a)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.575, 0.675, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'b)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.1, 0.35, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'c)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.575, 0.35, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'd)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_e = annotation('textbox', [0.1, 0.2, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'e)', 'LineStyle', 'none', 'FontSize', 16);
            textbox_f = annotation('textbox', [0.575, 0.2, 0.1, 0.1], 'Interpreter', 'LaTeX', 'String', 'f)', 'LineStyle', 'none', 'FontSize', 16);

            set(axs(1:4),'YScale', 'log');

            % set(axs(1:2), 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1])

            ylabel(axs(1), '$$\mu_0$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3), '$$\lambda$$ [s]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(5), '$$n$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(5:6), '$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            % legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            % legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            qrange_FB1 = [0 2];
            qrange_FB2 = [0 16];

            axis(axs(1),[qrange_FB1 4e-1 1e4])
            axis(axs(2),[qrange_FB2 4e-1 1e4])
            axis(axs(3),[qrange_FB1 1e-1 1e3])
            axis(axs(4),[qrange_FB2 1e-1 1e3])
            axis(axs(5),[qrange_FB1 0 1e0])
            axis(axs(6),[qrange_FB2 0 1e0])
            fig_out = AYfig_;
        end
        function fig_out = FB_Carreau_par_vs_q_compact(obj, AYfig_, FB1, FB2)
            axs = set_Carreau_params_vs_q_compact_axes(AYfig_);
            axs_m = [axs(1) axs(2)];
            axs_l = [axs(3) axs(4)];
            axs_n = [axs(5) axs(6)];

            axs_FB2_qlims = [axs(1) axs(3) axs(5)];
            axs_FB1_qlims = [axs(2) axs(4) axs(6)];

            qlims_FB1 = [0 2.05];
            qlims_FB2 = [0 15.2];

            m_lims = [4e-1 1e4];
            l_lims = [1e-1 1e3];
            n_lims = [0 1];

            box_m_vs_q = [  qlims_FB1(1) m_lims(1); ...
                            qlims_FB1(1) m_lims(2); ...
                            qlims_FB1(2) m_lims(2); ...
                            qlims_FB1(2) m_lims(1)];

            box_l_vs_q = [  qlims_FB1(1) l_lims(1); ...
                            qlims_FB1(1) l_lims(2); ...
                            qlims_FB1(2) l_lims(2); ...
                            qlims_FB1(2) l_lims(1)];

            box_n_vs_q = [  qlims_FB1(1) n_lims(1); ...
                            qlims_FB1(1) n_lims(2); ...
                            qlims_FB1(2) n_lims(2); ...
                            qlims_FB1(2) n_lims(1)];

            for ax_active = axs_m
                patch(ax_active, 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_m_vs_q, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
                for i = 1:length(FB1.exp)
                    legend_set_a(i) = plot(ax_active, FB1.exp(i).q, FB1.exp(i).mu0_Carreau, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
                end
                for i = 1:length(FB2.exp)
                    legend_set_b(i) = plot(ax_active, FB2.exp(i).q, FB2.exp(i).mu0_Carreau, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
                end
            end

            for ax_active = axs_l
                patch(ax_active, 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_l_vs_q, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
                for i = 1:length(FB1.exp)
                    legend_set_c(i) = plot(ax_active, FB1.exp(i).q, FB1.exp(i).lambda_Carreau, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
                end
                for i = 1:length(FB2.exp)
                    legend_set_d(i) = plot(ax_active, FB2.exp(i).q, FB2.exp(i).lambda_Carreau, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
                end
            end

            for ax_active = axs_n
                patch(ax_active, 'Faces', 1:4, 'LineStyle', ':', 'Vertices', box_n_vs_q, 'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',1.5);
                for i = 1:length(FB1.exp)
                    legend_set_c(i) = plot(ax_active, FB1.exp(i).q, FB1.exp(i).n_Carreau, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW_L, 'MarkerSize', FB1.MS_L, 'DisplayName', FB1.exp(i).label);
                end
                for i = 1:length(FB2.exp)
                    legend_set_d(i) = plot(ax_active, FB2.exp(i).q, FB2.exp(i).n_Carreau, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW_L, 'MarkerSize', FB2.MS_L, 'DisplayName', FB2.exp(i).label);
                end
            end

            textbox_a = annotation('textbox', [0.28 0.28 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(a)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_b = annotation('textbox', [0.19 0.73 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(b)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_c = annotation('textbox', [0.60 0.28 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(c)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_d = annotation('textbox', [0.52 0.73 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(d)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_e = annotation('textbox', [0.93 0.28 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(e)}', 'LineStyle', 'none', 'FontSize', 16);
            textbox_f = annotation('textbox', [0.84 0.73 0.1 0.1], 'Interpreter', 'LaTeX', 'String', '\textit{(f)}', 'LineStyle', 'none', 'FontSize', 16);

            set([axs_m axs_l],'YScale', 'log');

            set(axs_m(1), 'YTick', [1e0,1e1,1e2,1e3,1e4])
            set(axs_m(2), 'YTick', [1e0,1e2,1e4])
            set(axs_l(1), 'YTick', [1e-1,1e0,1e1,1e2,1e3])
            set(axs_l(2), 'YTick', [1e-1,1e1,1e3])

            ylabel(axs_m, '$$\mu_0$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs_l, '$$\lambda$$ [s]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs_n(1), '$$n$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs_n(2), '$$n$$', 'Interpreter', 'LaTeX','FontSize',12)

            xlabel(axs_FB2_qlims, '$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_FB1_qlims, '$$q = Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',12)

            % legend(axs(1), legend_set_a,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);
            % legend(axs(2), legend_set_b,'Location', 'NorthEast', 'Interpreter', 'Latex', 'NumColumns', 2);

            % qrange_FB1 = [0 2];
            % qrange_FB2 = [0 16];
            %
            ylim(axs_m, m_lims);
            ylim(axs_l, l_lims);
            ylim(axs_n, n_lims);
            xlim(axs_FB1_qlims, qlims_FB1);
            xlim(axs_FB2_qlims, qlims_FB2);

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
        function fig_out = FB1_FB2_Bingham_and_Carreau_fluid_params(obj,AYfig_,FB1,FB2)
            [tdim1,tdim2] = deal(2,2);
            axs=set_FB_Bingham_Carreau_axes(AYfig_);
            ax_B=0;
            ax_C=4;

            inds_full = 1:length(axs);
            FB1_inds = [1 3 5 6 7];
            FB2_inds = [2 4 8 9 10];

            axit=0;
            for FB = [FB1 FB2]
                explen=length(FB.exp);
                [tauy_vec mup_vec q_vec] = deal(nan(explen,1));
                [mu0_vec lambda_vec n_vec] = deal(nan(explen,1));
                color_mat=nan(explen,3);
                for i=1:(explen)
                    expi=FB.exp(i);

                    [tauy_vec(i) mup_vec(i)] = deal(expi.tau_y_Bingham/expi.tau_static, expi.mu_p_Bingham);
                    [mu0_vec(i) lambda_vec(i) n_vec(i)] = deal(expi.mu0_Carreau, expi.lambda_Carreau, expi.n_Carreau);

                    [color_mat(i,:) q_vec(i)] = deal(FB.exp(i).color, FB.exp(i).q);
                end
                scatter(axs(1+axit), q_vec, tauy_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',1*FB.MS_L*FB.MS_L);
                scatter(axs(3+axit), q_vec, mup_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',1*FB.MS_L*FB.MS_L);
                scatter(axs(1+ax_C), q_vec, mu0_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',1*FB.MS_L*FB.MS_L);
                scatter(axs(2+ax_C), q_vec, lambda_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',1*FB.MS_L*FB.MS_L);
                scatter(axs(3+ax_C), q_vec, n_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',1*FB.MS_L*FB.MS_L);

                ylabel(axs(1+ax_C), '$$\mu_0$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',14)
                ylabel(axs(2+ax_C), '$$\lambda$$ [s]', 'Interpreter', 'LaTeX','FontSize',14)
                ylabel(axs(3+ax_C), '$$n$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',14)

                ax_B=ax_B+2;
                ax_C=ax_C+3;
                axit = axit + 1;
            end
            set(axs,'YScale','log');
            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(1),'$$\tau_y/\tau_{q=0}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(3),'$$\tilde{\mu}_p$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',14)

            ylim([axs(1) axs(2)], [1e-4 1]); %% set yield stress ratio range
            ylim([axs(3) axs(4)], [1e-2 3e-1]); %% set plastic viscosity range
            set([axs(7) axs(10)],'YScale','linear');
            ylim([axs(7) axs(10)], [0 0.6]); %% set thinning index range

            xlim(axs(FB1_inds), [0 2.1]);
            xlim(axs(FB2_inds), [0 15.2]);

            fig_out=AYfig_;
        end
        function fig_out = FB1_dimensional_regime_plot(obj,AYfig_,FB1, label_flag_)
            if (nargin==3)
                label_flag = 'names_only';
            else
                label_flag = label_flag_;
            end

            run figure_properties.m
            color_ttv = [255 220 0]/255;
            color_fgm = [0 190 0]/255;
            color_sgf = [60 255 255]/255;
            color_dgf = [0 160 255]/255;
            color_frs = [0 228 116]/255;
            ncolor = 30;

            cgrad_FGM_2_TTV = make_color_gradient(color_fgm,color_ttv,ncolor);
            cgrad_FRS_2_TTV = make_color_gradient(color_frs,color_ttv,ncolor);
            cgrad_DGF_2_FRS = make_color_gradient(color_dgf,color_frs,ncolor);
            cgrad_DGF_2_SGF = make_color_gradient(color_dgf,color_sgf,ncolor);
            cgrad_SGF_2_FGM = make_color_gradient(color_sgf,color_fgm,ncolor);

            % AYfig_.init_tiles([1,2]);
            % axs = AYfig_.ax_tile;
            axs = AYfig_.ax;

            % set(axs, 'Units', 'centimeters');
            % axs.Position = [axs.Position(1:2) obj.dim1_graphical_abstract]

            axs.Position = [axs.Position(1:2) 0.75 0.625];

            hold(axs, 'on');
            box(axs,'on');
            olims = [1e-2 2e3];
            % qlims_mat = [0 3.0; 0 3.0];
            % qlims_mat_plot = [0 3.0; 0 3.0];
            qlims_mat = [0 2.0; 0 2.0];
            qlims_mat_plot = [0 2.0; 0 2.0];


            txFB1_row1 = [  0.275 0.525 0.1 0.1; ...
                            0.275 0.350 0.1 0.1; ...
                            0.275 0.150 0.1 0.1; ...
                            0.630 0.500 0.1 0.1; ...
                            0.675 0.250 0.1 0.1];
            txFB2_row1 = [  0.565 0.750 0.1 0.1; ...
                            0.565 0.500 0.1 0.1; ...
                            0.565 0.200 0.1 0.1; ...
                            0.795 0.650 0.1 0.1; ...
                            0.810 0.300 0.1 0.1];

            tx_row1_tens = nan(size(txFB1_row1,1),size(txFB1_row1,2),2);
            tx_row1_tens(:,:,1) = txFB1_row1;
            tx_row1_tens(:,:,2) = txFB2_row1;

            f_s = 20;

            iFB=1;
            % for FB = [FB1 FB2]
            for FB = FB1

                tx_row1 = tx_row1_tens(:,:,iFB);
                tx_row2 = tx_row1 - [0 0.05 0 0];
                tx_row2(end,:) = tx_row2(end,:) + [0.01 0 0 0];

                tx_row3 = tx_row2 - [0 0.05 0 0];
                tx_row3(end,:) = tx_row3(end,:) + [0.01 0 0 0];

                tx_row4 = tx_row3 - [0 0.05 0 0];
                tx_row4(end,:) = tx_row4(end,:) + [0.01 0 0 0];

                qlims = qlims_mat(iFB,:);
                qcrit_sgf = FB.exp(iFB).qcrit_sgf;
                qcrit_ttv = FB.exp(iFB).qcrit_ttv;
                qcrit_fgm = FB.exp(iFB).qcrit_fgm;

                ocrit_dgf_ql = FB.exp(iFB).ocrit_dgf_ql;
                ocrit_dgf_qh = FB.exp(iFB).ocrit_dgf_qh;
                ocrit_fgm_ql = FB.exp(iFB).ocrit_fgm_ql;
                ocrit_fgm_qh = FB.exp(iFB).ocrit_fgm_qh;

                ocrit_3 = ((ocrit_fgm_ql - ocrit_dgf_qh)/(qcrit_fgm-qcrit_sgf))*(qcrit_ttv-qcrit_sgf) + ocrit_dgf_qh;
                qcrit_3_dgf_2_sgf = qcrit_sgf - (0.5*(qcrit_ttv-qcrit_sgf));
                ocrit_3_dgf = ((ocrit_dgf_qh - ocrit_dgf_ql)/(qcrit_sgf-qlims(1)))*(qcrit_3_dgf_2_sgf-qlims(1)) + ocrit_dgf_ql;

                dgf = [ olims(1) qlims(1); ...
                        olims(1) qcrit_sgf; ...
                        ocrit_dgf_qh qcrit_sgf; ...
                        ocrit_dgf_ql, qlims(1)];
                sgf = [ dgf(2,:); ...
                        olims(1) qcrit_fgm; ...
                        ocrit_fgm_ql qcrit_fgm; ...
                        ocrit_3 qcrit_ttv; ...
                        dgf(3,:)];
                fgm = [ sgf(2,:); ...
                        olims(1) qlims(2); ...
                        ocrit_fgm_qh qlims(2); ...
                        sgf(3,:)];

                frs = [ olims(2) qcrit_ttv; ...
                        olims(2) qlims(1); ...
                        ocrit_dgf_ql qlims(1); ...
                        ocrit_dgf_qh qcrit_sgf; ...
                        ocrit_3 qcrit_ttv];
                ttv = [ olims(2) qlims(2); ...
                        olims(2) qcrit_ttv; ...
                        ocrit_3 qcrit_ttv; ...
                        ocrit_fgm_ql qcrit_fgm; ...
                        ocrit_fgm_qh qlims(2)];

                ofac = 2;

                fgm_2_ttv = [ocrit_fgm_qh/ofac qlims(2); ...
                             ocrit_fgm_ql/ofac qcrit_fgm; ...
                             ocrit_fgm_ql qcrit_fgm; ...
                             ocrit_fgm_qh qlims(2)];
                sgf_2_fgm = [olims(1) qcrit_fgm; ...
                             ocrit_fgm_ql/ofac qcrit_fgm; ...
                             ocrit_fgm_ql qcrit_fgm; ...
                             ocrit_3 qcrit_ttv; ...
                             ocrit_dgf_qh qcrit_sgf; ...
                             ocrit_fgm_ql qcrit_ttv; ...
                             ocrit_fgm_ql/ofac qcrit_ttv; ...
                             olims(1) qcrit_ttv];
                frs_2_ttv = [ocrit_3 qcrit_ttv; ...
                             ocrit_dgf_qh qcrit_ttv; ...
                             olims(2) qcrit_ttv; ...
                             % olims(2) 0.5*(qcrit_sgf+qcrit_ttv); ...
                             % olims(2) qcrit_sgf; ...
                             olims(2) qcrit_3_dgf_2_sgf; ...
                             olims(2) qlims(1); ...
                             ocrit_dgf_ql qlims(1); ...
                             ocrit_dgf_qh qcrit_sgf];
                dgf_2_sgf = [olims(1) qcrit_sgf; ...
                             ocrit_3 qcrit_sgf; ...
                             ocrit_dgf_qh qcrit_sgf; ...
                             ocrit_3_dgf qcrit_3_dgf_2_sgf; ...
                             ocrit_dgf_ql qlims(1); ...
                             ocrit_dgf_qh qcrit_3_dgf_2_sgf; ...
                             ocrit_3 qcrit_3_dgf_2_sgf; ...
                             olims(1) qcrit_3_dgf_2_sgf];

                c_fgm_2_ttv = [ color_fgm; ...
                                color_fgm; ...
                                color_ttv; ...
                                color_ttv];
                c_sgf_2_fgm = [ color_fgm; ...
                                color_fgm; ...
                                color_ttv; ...
                                color_ttv; ...
                                color_ttv; ... % color_frs; ...
                                color_sgf; ...
                                color_sgf; ...
                                color_sgf];
                c_frs_2_ttv = [ color_ttv; ...
                                color_ttv; ...
                                color_ttv; ...
                                % color_ttv; ...
                                % color_ttv; ...
                                color_ttv; ... % color_frs; ...
                                color_frs; ...
                                color_frs; ...
                                color_ttv]; % color_frs];
                c_dgf_2_sgf = [ color_sgf; ...
                                color_sgf; ...
                                color_ttv; ... % color_frs; ...
                                color_frs; ...
                                color_dgf; ...
                                color_dgf; ...
                                color_dgf; ...
                                color_dgf];

                abb = 0;
                fa = 1;

                patch(axs(iFB), 'Faces', 1:size(dgf,1), 'Vertices', dgf, 'FaceColor',color_dgf,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(sgf,1), 'Vertices', sgf, 'FaceColor',color_sgf,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(fgm,1), 'Vertices', fgm, 'FaceColor',color_fgm,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(frs,1), 'Vertices', frs, 'FaceColor',color_frs,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(ttv,1), 'Vertices', ttv, 'FaceColor',color_ttv,'EdgeAlpha',abb,'FaceAlpha',fa);

                ab = 0;
                mrk='none';
                mrk2='o';

                patch(axs(iFB), 'Faces', 1:size(fgm_2_ttv,1), 'Vertices', fgm_2_ttv, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_fgm_2_ttv,'Marker',mrk);
                patch(axs(iFB), 'Faces', 1:size(sgf_2_fgm,1), 'Vertices', sgf_2_fgm, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_sgf_2_fgm,'Marker',mrk);
                patch(axs(iFB), 'Faces', 1:size(frs_2_ttv,1), 'Vertices', frs_2_ttv, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_frs_2_ttv,'Marker',mrk);
                patch(axs(iFB), 'Faces', 1:size(dgf_2_sgf,1), 'Vertices', dgf_2_sgf, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_dgf_2_sgf,'Marker',mrk);

                l_s = ':';
                abl = 1;
                l_w = 1.5;

                patch(axs(iFB), 'Faces', 1:size(frs,1), 'LineStyle', l_s, 'LineWidth', l_w, 'Vertices', frs, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(ttv,1), 'LineStyle', l_s, 'LineWidth', l_w, 'Vertices', ttv, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(dgf,1), 'LineStyle', l_s, 'LineWidth', l_w, 'Vertices', dgf, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(sgf,1), 'LineStyle', l_s, 'LineWidth', l_w, 'Vertices', sgf, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(fgm,1), 'LineStyle', l_s, 'LineWidth', l_w, 'Vertices', fgm, 'FaceAlpha',0,'EdgeAlpha',abl);

                if (strcmp(label_flag,'names_only'))
                    textbox_FGM_row1 = annotation('textbox',tx_row1(1,:), 'Interpreter', 'LaTeX', ...
                    'String', 'FGM', 'LineStyle', 'none', 'FontSize', f_s);
                    textbox_SGF_row1 = annotation('textbox',tx_row1(2,:), 'Interpreter', 'LaTeX', ...
                    'String', 'SGF', 'LineStyle', 'none', 'FontSize', f_s);
                    textbox_DGF_row1 = annotation('textbox',tx_row1(3,:), 'Interpreter', 'LaTeX', ...
                    'String', 'DGF', 'LineStyle', 'none', 'FontSize', f_s);
                    textbox_TTV_row1 = annotation('textbox',tx_row1(4,:), 'Interpreter', 'LaTeX', ...
                    'String', 'TTV', 'LineStyle', 'none', 'FontSize', f_s);
                    textbox_FRS_row1 = annotation('textbox',tx_row1(5,:), 'Interpreter', 'LaTeX', ...
                    'String', 'FRS', 'LineStyle', 'none', 'FontSize', f_s);

                    % textbox_FGM_row2 = annotation('textbox',tx_row2(1,:), 'Interpreter', 'LaTeX', ...
                    % 'String', 'material (FGM)', 'LineStyle', 'none', 'FontSize', 12);
                    % textbox_SGF_row2 = annotation('textbox',tx_row2(2,:), 'Interpreter', 'LaTeX', ...
                    % 'String', 'flow (SGF)', 'LineStyle', 'none', 'FontSize', 12);
                    % textbox_DGF_row2 = annotation('textbox',tx_row2(3,:), 'Interpreter', 'LaTeX', ...
                    % 'String', 'flow (DGF)', 'LineStyle', 'none', 'FontSize', 12);
                    % textbox_TTV_row2 = annotation('textbox',tx_row2(4,:), 'Interpreter', 'LaTeX', ...
                    % 'String', 'toroidal vortices', 'LineStyle', 'none', 'FontSize', 12);
                    % textbox_frs_row2 = annotation('textbox',tx_row2(5,:), 'Interpreter', 'LaTeX', ...
                    % 'String', 'suspension', 'LineStyle', 'none', 'FontSize', 12);
                    %
                    % textbox_TTV_row3 = annotation('textbox',tx_row3(4,:), 'Interpreter', 'LaTeX', ...
                    % 'String', '(TTV)', 'LineStyle', 'none', 'FontSize', 12);
                    % textbox_FRS_row3 = annotation('textbox',tx_row3(5,:), 'Interpreter', 'LaTeX', ...
                    % 'String', '(FRS)', 'LineStyle', 'none', 'FontSize', 12);

                else
                    textbox_FGM_row1 = annotation('textbox',tx_row1(1,:), 'Interpreter', 'LaTeX', ...
                    'String', 'FGM: $$ \tau \propto \omega_i^{\alpha}, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_SGF_row1 = annotation('textbox',tx_row1(2,:), 'Interpreter', 'LaTeX', ...
                    'String', 'SGF: $$ \tau \neq f(\omega_i), $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_DGF_row1 = annotation('textbox',tx_row1(3,:), 'Interpreter', 'LaTeX', ...
                    'String', 'DGF: $$ \tau \approx \tau_{q = 0} \cdot ( 1 - q ), $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row1 = annotation('textbox',tx_row1(4,:), 'Interpreter', 'LaTeX', ...
                    'String', 'TTV: $$ \tau \propto \omega_i^{\alpha}, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_FRS_row1 = annotation('textbox',tx_row1(5,:), 'Interpreter', 'LaTeX', ...
                    'String', 'FRS: $$ \tau \propto \omega_i^{\alpha}, $$', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_FGM_row2 = annotation('textbox',tx_row2(1,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0 < | \alpha | < 1, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_SGF_row2 = annotation('textbox',tx_row2(2,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.55 < \phi < 0.57 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_DGF_row2 = annotation('textbox',tx_row2(3,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.57 < \phi < 0.58 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row2 = annotation('textbox',tx_row2(4,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ \alpha > 1.5, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_frs_row2 = annotation('textbox',tx_row2(5,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0 < \alpha < 1, $$', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_DGF_row3 = annotation('textbox',tx_row3(1,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.51 < \phi < 0.55 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row3 = annotation('textbox',tx_row3(4,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.51 < \phi < 0.55 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_FRS_row3 = annotation('textbox',tx_row3(5,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ \phi > 0.55, $$', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_FRS_row4 = annotation('textbox',tx_row4(5,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ \phi < 0.57 $$', 'LineStyle', 'none', 'FontSize', 12);
                end

                iFB = iFB + 1;
            end

            % title(axs(1), '\textit{FB1}, $$ 113 $$ $$\mu$$m', 'Interpreter', 'LaTeX', 'FontSize', 14);
            % title(axs(2), '\textit{FB2}, $$ 49 $$ $$\mu$$m', 'Interpreter', 'LaTeX', 'FontSize', 14);

            set(axs,'YScale','linear','XScale','log');
            % xlabel(axs, '$$\omega_i$$', 'Interpreter', 'LaTeX','FontSize',f_s)
            % ylabel(axs(1), '$$q$$', 'Interpreter', 'LaTeX','FontSize',f_s)
            xlim(axs, olims);
            % set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2,1e3])
            set(axs, 'XTick', [])
            set(axs, 'XTickLabel', [])
            set(axs, 'YTick', [])
            set(axs, 'YTickLabel', [])
            ylim(axs(1), qlims_mat_plot(1,:));
            % ylim(axs(2), qlims_mat_plot(2,:));

            fig_out = AYfig_;
        end
        function fig_out = FB_dimensional_regime_plot(obj,AYfig_,FB1, FB2, label_flag_)
            if (nargin==4)
                label_flag = 'names_only';
            else
                label_flag = label_flag_;
            end

            run figure_properties.m
            color_ttv = [255 220 0]/255;
            color_fgm = [0 190 0]/255;
            color_sgf = [60 255 255]/255;
            color_dgf = [0 160 255]/255;
            color_frs = [0 228 116]/255;
            ncolor = 30;

            cgrad_FGM_2_TTV = make_color_gradient(color_fgm,color_ttv,ncolor);
            cgrad_FRS_2_TTV = make_color_gradient(color_frs,color_ttv,ncolor);
            cgrad_DGF_2_FRS = make_color_gradient(color_dgf,color_frs,ncolor);
            cgrad_DGF_2_SGF = make_color_gradient(color_dgf,color_sgf,ncolor);
            cgrad_SGF_2_FGM = make_color_gradient(color_sgf,color_fgm,ncolor);

            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');
            olims = [1e-2 2e3];
            qlims_mat = [0 3.0; 0 3.0];
            qlims_mat_plot = [0 3.0; 0 3.0];


            txFB1_row1 = [  0.090 0.650 0.1 0.1; ...
                            0.090 0.350 0.1 0.1; ...
                            0.090 0.200 0.1 0.1; ...
                            0.320 0.650 0.1 0.1; ...
                            0.340 0.300 0.1 0.1];
            txFB2_row1 = [  0.565 0.750 0.1 0.1; ...
                            0.565 0.500 0.1 0.1; ...
                            0.565 0.200 0.1 0.1; ...
                            0.795 0.650 0.1 0.1; ...
                            0.810 0.300 0.1 0.1];

            tx_row1_tens = nan(size(txFB1_row1,1),size(txFB1_row1,2),2);
            tx_row1_tens(:,:,1) = txFB1_row1;
            tx_row1_tens(:,:,2) = txFB2_row1;

            iFB=1;
            for FB = [FB1 FB2]

                tx_row1 = tx_row1_tens(:,:,iFB);
                tx_row2 = tx_row1 - [0 0.05 0 0];
                tx_row2(end,:) = tx_row2(end,:) + [0.01 0 0 0];

                tx_row3 = tx_row2 - [0 0.05 0 0];
                tx_row3(end,:) = tx_row3(end,:) + [0.01 0 0 0];

                tx_row4 = tx_row3 - [0 0.05 0 0];
                tx_row4(end,:) = tx_row4(end,:) + [0.01 0 0 0];

                qlims = qlims_mat(iFB,:);
                qcrit_sgf = FB.exp(iFB).qcrit_sgf;
                qcrit_ttv = FB.exp(iFB).qcrit_ttv;
                qcrit_fgm = FB.exp(iFB).qcrit_fgm;

                ocrit_dgf_ql = FB.exp(iFB).ocrit_dgf_ql;
                ocrit_dgf_qh = FB.exp(iFB).ocrit_dgf_qh;
                ocrit_fgm_ql = FB.exp(iFB).ocrit_fgm_ql;
                ocrit_fgm_qh = FB.exp(iFB).ocrit_fgm_qh;

                ocrit_3 = ((ocrit_fgm_ql - ocrit_dgf_qh)/(qcrit_fgm-qcrit_sgf))*(qcrit_ttv-qcrit_sgf) + ocrit_dgf_qh;
                qcrit_3_dgf_2_sgf = qcrit_sgf - (0.5*(qcrit_ttv-qcrit_sgf));
                ocrit_3_dgf = ((ocrit_dgf_qh - ocrit_dgf_ql)/(qcrit_sgf-qlims(1)))*(qcrit_3_dgf_2_sgf-qlims(1)) + ocrit_dgf_ql;

                dgf = [ olims(1) qlims(1); ...
                        olims(1) qcrit_sgf; ...
                        ocrit_dgf_qh qcrit_sgf; ...
                        ocrit_dgf_ql, qlims(1)];
                sgf = [ dgf(2,:); ...
                        olims(1) qcrit_fgm; ...
                        ocrit_fgm_ql qcrit_fgm; ...
                        ocrit_3 qcrit_ttv; ...
                        dgf(3,:)];
                fgm = [ sgf(2,:); ...
                        olims(1) qlims(2); ...
                        ocrit_fgm_qh qlims(2); ...
                        sgf(3,:)];

                frs = [ olims(2) qcrit_ttv; ...
                        olims(2) qlims(1); ...
                        ocrit_dgf_ql qlims(1); ...
                        ocrit_dgf_qh qcrit_sgf; ...
                        ocrit_3 qcrit_ttv];
                ttv = [ olims(2) qlims(2); ...
                        olims(2) qcrit_ttv; ...
                        ocrit_3 qcrit_ttv; ...
                        ocrit_fgm_ql qcrit_fgm; ...
                        ocrit_fgm_qh qlims(2)];

                ofac = 2;

                fgm_2_ttv = [ocrit_fgm_qh/ofac qlims(2); ...
                             ocrit_fgm_ql/ofac qcrit_fgm; ...
                             ocrit_fgm_ql qcrit_fgm; ...
                             ocrit_fgm_qh qlims(2)];
                sgf_2_fgm = [olims(1) qcrit_fgm; ...
                             ocrit_fgm_ql/ofac qcrit_fgm; ...
                             ocrit_fgm_ql qcrit_fgm; ...
                             ocrit_3 qcrit_ttv; ...
                             ocrit_dgf_qh qcrit_sgf; ...
                             ocrit_fgm_ql qcrit_ttv; ...
                             ocrit_fgm_ql/ofac qcrit_ttv; ...
                             olims(1) qcrit_ttv];
                frs_2_ttv = [ocrit_3 qcrit_ttv; ...
                             ocrit_dgf_qh qcrit_ttv; ...
                             olims(2) qcrit_ttv; ...
                             % olims(2) 0.5*(qcrit_sgf+qcrit_ttv); ...
                             % olims(2) qcrit_sgf; ...
                             olims(2) qcrit_3_dgf_2_sgf; ...
                             olims(2) qlims(1); ...
                             ocrit_dgf_ql qlims(1); ...
                             ocrit_dgf_qh qcrit_sgf];
                dgf_2_sgf = [olims(1) qcrit_sgf; ...
                             ocrit_3 qcrit_sgf; ...
                             ocrit_dgf_qh qcrit_sgf; ...
                             ocrit_3_dgf qcrit_3_dgf_2_sgf; ...
                             ocrit_dgf_ql qlims(1); ...
                             ocrit_dgf_qh qcrit_3_dgf_2_sgf; ...
                             ocrit_3 qcrit_3_dgf_2_sgf; ...
                             olims(1) qcrit_3_dgf_2_sgf];

                c_fgm_2_ttv = [ color_fgm; ...
                                color_fgm; ...
                                color_ttv; ...
                                color_ttv];
                c_sgf_2_fgm = [ color_fgm; ...
                                color_fgm; ...
                                color_ttv; ...
                                color_ttv; ...
                                color_ttv; ... % color_frs; ...
                                color_sgf; ...
                                color_sgf; ...
                                color_sgf];
                c_frs_2_ttv = [ color_ttv; ...
                                color_ttv; ...
                                color_ttv; ...
                                % color_ttv; ...
                                % color_ttv; ...
                                color_ttv; ... % color_frs; ...
                                color_frs; ...
                                color_frs; ...
                                color_ttv]; % color_frs];
                c_dgf_2_sgf = [ color_sgf; ...
                                color_sgf; ...
                                color_ttv; ... % color_frs; ...
                                color_frs; ...
                                color_dgf; ...
                                color_dgf; ...
                                color_dgf; ...
                                color_dgf];

                abb = 0;
                fa = 1;

                patch(axs(iFB), 'Faces', 1:size(dgf,1), 'Vertices', dgf, 'FaceColor',color_dgf,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(sgf,1), 'Vertices', sgf, 'FaceColor',color_sgf,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(fgm,1), 'Vertices', fgm, 'FaceColor',color_fgm,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(frs,1), 'Vertices', frs, 'FaceColor',color_frs,'EdgeAlpha',abb,'FaceAlpha',fa);
                patch(axs(iFB), 'Faces', 1:size(ttv,1), 'Vertices', ttv, 'FaceColor',color_ttv,'EdgeAlpha',abb,'FaceAlpha',fa);

                ab = 0;
                mrk='none';
                mrk2='o';

                patch(axs(iFB), 'Faces', 1:size(fgm_2_ttv,1), 'Vertices', fgm_2_ttv, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_fgm_2_ttv,'Marker',mrk);
                patch(axs(iFB), 'Faces', 1:size(sgf_2_fgm,1), 'Vertices', sgf_2_fgm, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_sgf_2_fgm,'Marker',mrk);
                patch(axs(iFB), 'Faces', 1:size(frs_2_ttv,1), 'Vertices', frs_2_ttv, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_frs_2_ttv,'Marker',mrk);
                patch(axs(iFB), 'Faces', 1:size(dgf_2_sgf,1), 'Vertices', dgf_2_sgf, ...
                'FaceColor','interp','EdgeAlpha',ab,'FaceVertexCData',c_dgf_2_sgf,'Marker',mrk);

                l_s = ':';
                abl = 1;

                patch(axs(iFB), 'Faces', 1:size(frs,1), 'LineStyle', l_s, 'Vertices', frs, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(ttv,1), 'LineStyle', l_s, 'Vertices', ttv, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(dgf,1), 'LineStyle', l_s, 'Vertices', dgf, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(sgf,1), 'LineStyle', l_s, 'Vertices', sgf, 'FaceAlpha',0,'EdgeAlpha',abl);
                patch(axs(iFB), 'Faces', 1:size(fgm,1), 'LineStyle', l_s, 'Vertices', fgm, 'FaceAlpha',0,'EdgeAlpha',abl);

                if (strcmp(label_flag,'names_only'))
                    textbox_FGM_row1 = annotation('textbox',tx_row1(1,:), 'Interpreter', 'LaTeX', ...
                    'String', 'Fluidized granular', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_SGF_row1 = annotation('textbox',tx_row1(2,:), 'Interpreter', 'LaTeX', ...
                    'String', 'Sheared granular', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_DGF_row1 = annotation('textbox',tx_row1(3,:), 'Interpreter', 'LaTeX', ...
                    'String', 'Dense granular', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row1 = annotation('textbox',tx_row1(4,:), 'Interpreter', 'LaTeX', ...
                    'String', 'Turbulent', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_FRS_row1 = annotation('textbox',tx_row1(5,:), 'Interpreter', 'LaTeX', ...
                    'String', 'Frictional', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_FGM_row2 = annotation('textbox',tx_row2(1,:), 'Interpreter', 'LaTeX', ...
                    'String', 'material (FGM)', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_SGF_row2 = annotation('textbox',tx_row2(2,:), 'Interpreter', 'LaTeX', ...
                    'String', 'flow (SGF)', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_DGF_row2 = annotation('textbox',tx_row2(3,:), 'Interpreter', 'LaTeX', ...
                    'String', 'flow (DGF)', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row2 = annotation('textbox',tx_row2(4,:), 'Interpreter', 'LaTeX', ...
                    'String', 'toroidal vortices', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_frs_row2 = annotation('textbox',tx_row2(5,:), 'Interpreter', 'LaTeX', ...
                    'String', 'suspension', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_TTV_row3 = annotation('textbox',tx_row3(4,:), 'Interpreter', 'LaTeX', ...
                    'String', '(TTV)', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_FRS_row3 = annotation('textbox',tx_row3(5,:), 'Interpreter', 'LaTeX', ...
                    'String', '(FRS)', 'LineStyle', 'none', 'FontSize', 12);

                else
                    textbox_FGM_row1 = annotation('textbox',tx_row1(1,:), 'Interpreter', 'LaTeX', ...
                    'String', 'FGM: $$ \tau \propto \omega_i^{\alpha}, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_SGF_row1 = annotation('textbox',tx_row1(2,:), 'Interpreter', 'LaTeX', ...
                    'String', 'SGF: $$ \tau \neq f(\omega_i), $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_DGF_row1 = annotation('textbox',tx_row1(3,:), 'Interpreter', 'LaTeX', ...
                    'String', 'DGF: $$ \tau \approx \tau_{q = 0} \cdot ( 1 - q ), $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row1 = annotation('textbox',tx_row1(4,:), 'Interpreter', 'LaTeX', ...
                    'String', 'TTV: $$ \tau \propto \omega_i^{\alpha}, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_FRS_row1 = annotation('textbox',tx_row1(5,:), 'Interpreter', 'LaTeX', ...
                    'String', 'FRS: $$ \tau \propto \omega_i^{\alpha}, $$', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_FGM_row2 = annotation('textbox',tx_row2(1,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0 < | \alpha | < 1, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_SGF_row2 = annotation('textbox',tx_row2(2,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.55 < \phi < 0.57 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_DGF_row2 = annotation('textbox',tx_row2(3,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.57 < \phi < 0.58 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row2 = annotation('textbox',tx_row2(4,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ \alpha > 1.5, $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_frs_row2 = annotation('textbox',tx_row2(5,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0 < \alpha < 1, $$', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_DGF_row3 = annotation('textbox',tx_row3(1,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.51 < \phi < 0.55 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_TTV_row3 = annotation('textbox',tx_row3(4,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ 0.51 < \phi < 0.55 $$', 'LineStyle', 'none', 'FontSize', 12);
                    textbox_FRS_row3 = annotation('textbox',tx_row3(5,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ \phi > 0.55, $$', 'LineStyle', 'none', 'FontSize', 12);

                    textbox_FRS_row4 = annotation('textbox',tx_row4(5,:), 'Interpreter', 'LaTeX', ...
                    'String', '$$ \phi < 0.57 $$', 'LineStyle', 'none', 'FontSize', 12);
                end

                iFB = iFB + 1;
            end

            title(axs(1), '\textit{FB1}, $$ 113 $$ $$\mu$$m', 'Interpreter', 'LaTeX', 'FontSize', 14);
            title(axs(2), '\textit{FB2}, $$ 49 $$ $$\mu$$m', 'Interpreter', 'LaTeX', 'FontSize', 14);

            set(axs,'YScale','linear','XScale','log');
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(1), '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            xlim(axs, olims);
            set(axs, 'XTick', [1e-2,1e-1,1e0,1e1,1e2,1e3])
            ylim(axs(1), qlims_mat_plot(1,:));
            ylim(axs(2), qlims_mat_plot(2,:));

            fig_out = AYfig_;
        end
        function fig_out = FB_dimensionless_regime_plot(obj,AYfig_,FB1, FB2)
            run figure_properties.m

            color_ttv = [255 220 0]/255;
            color_fgm = [0 190 0]/255;
            color_sgf = [60 255 255]/255;
            color_dgf = [0 160 255]/255;
            color_frs = [0 228 116]/255;
            ncolor = 30;

            AYfig_.init_tiles([1,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');
            % olims = [1e-2 1.5e2];
            olims = [1e-2 1e3];
            qlims_mat = [0 3.0; 0 3.0];
            qlims_mat_plot = [0 3.0; 0 3.0];

            iFB=1;
            for FB = [FB1 FB2]
                qlims = qlims_mat(iFB,:);
                qcrit_sgf = FB.exp(iFB).qcrit_sgf;
                qcrit_ttv = FB.exp(iFB).qcrit_ttv;
                qcrit_fgm = FB.exp(iFB).qcrit_fgm;

                ocrit_dgf_ql = FB.exp(iFB).ocrit_dgf_ql;
                ocrit_dgf_qh = FB.exp(iFB).ocrit_dgf_qh;
                ocrit_fgm_ql = FB.exp(iFB).ocrit_fgm_ql;
                ocrit_fgm_qh = FB.exp(iFB).ocrit_fgm_qh;

                ocrit_3 = ((ocrit_fgm_ql - ocrit_dgf_qh)/(qcrit_fgm-qcrit_sgf))*(qcrit_ttv-qcrit_sgf) + ocrit_dgf_qh;

                dgf = [ olims(1) qlims(1); ...
                        olims(1) qcrit_sgf; ...
                        ocrit_dgf_qh qcrit_sgf; ...
                        ocrit_dgf_ql, qlims(1)];
                sgf = [ dgf(2,:); ...
                        olims(1) qcrit_fgm; ...
                        ocrit_fgm_ql qcrit_fgm; ...
                        ocrit_3 qcrit_ttv; ...
                        dgf(3,:)];
                fgm = [ sgf(2,:); ...
                        olims(1) qlims(2); ...
                        ocrit_fgm_qh qlims(2); ...
                        sgf(3,:)];

                frs = [ olims(2) qcrit_ttv; ...
                        olims(2) qlims(1); ...
                        ocrit_dgf_ql qlims(1); ...
                        ocrit_dgf_qh qcrit_sgf; ...
                        ocrit_3 qcrit_ttv];
                ttv = [ olims(2) qlims(2); ...
                        olims(2) qcrit_ttv; ...
                        ocrit_3 qcrit_ttv; ...
                        ocrit_fgm_ql qcrit_fgm; ...
                        ocrit_fgm_qh qlims(2)];

                patch(axs(iFB), 'Faces', [1 2 3 4], 'Vertices', dgf, 'FaceColor',color_dgf,'EdgeAlpha',0);
                patch(axs(iFB), 'Faces', [1 2 3 4 5], 'Vertices', sgf, 'FaceColor',color_sgf,'EdgeAlpha',0);
                patch(axs(iFB), 'Faces', [1 2 3 4], 'Vertices', fgm, 'FaceColor',color_fgm,'EdgeAlpha',0);
                patch(axs(iFB), 'Faces', [1 2 3 4 5], 'Vertices', frs, 'FaceColor',color_frs,'EdgeAlpha',0);
                patch(axs(iFB), 'Faces', [1 2 3 4 5], 'Vertices', ttv, 'FaceColor',color_ttv,'EdgeAlpha',0);

                iFB = iFB + 1;
            end

            title(axs(1), '$$ 113 $$ micron', 'Interpreter', 'LaTeX', 'FontSize', 14);
            title(axs(2), '$$ 49 $$ micron', 'Interpreter', 'LaTeX', 'FontSize', 14);

            set(axs,'YScale','linear','XScale','log');
            xlabel(axs, '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(1), '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            xlim(axs, olims);
            ylim(axs(1), qlims_mat_plot(1,:));
            ylim(axs(2), qlims_mat_plot(2,:));

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
function ax_out = set_FB1_phi_vs_omegai_vs_q_slices_axes(AYfig_)
    ax_out = gobjects(4, 1);
    figure(AYfig_.fig.Number);
    ax_out(1) = axes('Position', [0.100 0.525 0.800 0.470]);
    ax_out(2) = axes('Position', [0.075 0.100 0.250 0.270]);
    ax_out(3) = axes('Position', [0.400 0.100 0.250 0.270]);
    ax_out(4) = axes('Position', [0.730 0.100 0.250 0.270]);
    hold(ax_out, 'on');
end
function ax_out = set_taustar_vs_Gamma_compact_axes(AYfig_)
    ax_out = gobjects(2,1);
    figure(AYfig_.fig.Number);
    ax_out(1) = AYfig_.ax;
    ax_out(2) = axes('Position', [0.4 0.6 0.3 0.3]);
    box(ax_out(1), 'on');
    hold(ax_out, 'on');
end
function ax_out = set_Carreau_params_vs_q_compact_axes(AYfig_)
    ax_out = gobjects(6, 1);
    figure(AYfig_.fig.Number);
    ax_out(1) = axes('Position', [0.075 0.150 0.255 0.625]);
    ax_out(3) = axes('Position', [0.400 0.150 0.255 0.625]);
    ax_out(5) = axes('Position', [0.730 0.150 0.255 0.625]);
    ax_out(2) = axes('Position', [0.180 0.725 0.160 0.250]);
    ax_out(4) = axes('Position', [0.510 0.725 0.160 0.250]);
    ax_out(6) = axes('Position', [0.830 0.725 0.160 0.250]);
    hold(ax_out, 'on');
end
function ax_out = set_mup_tyr_vs_q_compact_axes(AYfig_)
    ax_out = gobjects(4, 1);
    figure(AYfig_.fig.Number);
    ax_out(2) = axes('Position', [0.075 0.15 0.4 0.8]);
    ax_out(4) = axes('Position', [0.575 0.15 0.4 0.8]);
    ax_out(1) = axes('Position', [0.2225 0.675 0.270 0.30]);
    ax_out(3) = axes('Position', [0.7225 0.675 0.270 0.30]);
    hold(ax_out, 'on');
end
function ax_out = set_FB_Bingham_Grat_vs_Reb_axes(AYfig_)
    ax_out = gobjects(6, 1);
    figure(AYfig_.fig.Number);
    ax_out(1) = subplot(2,2,1);
    ax_out(2) = subplot(2,2,2);
    ax_out(3) = subplot(4,2,5);
    ax_out(4) = subplot(4,2,6);
    ax_out(5) = subplot(4,2,7);
    ax_out(6) = subplot(4,2,8);
    hold(ax_out, 'on');
    box(ax_out, 'on');
end
function ax_out = set_FB_Bingham_Carreau_axes(AYfig_)
    ax_out = gobjects(10, 1);
    figure(AYfig_.fig.Number);
    ax_out(1) = subplot(4,2,1);
    ax_out(2) = subplot(4,2,2);
    ax_out(3) = subplot(4,2,3);
    ax_out(4) = subplot(4,2,4);
    ax_out(5) = subplot(4,3,7);
    ax_out(6) = subplot(4,3,8);
    ax_out(7) = subplot(4,3,9);
    ax_out(8) = subplot(4,3,10);
    ax_out(9) = subplot(4,3,11);
    ax_out(10) = subplot(4,3,12);
    hold(ax_out, 'on');
    box(ax_out, 'on');
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

function colors_out = make_color_gradient(color1_, color2_, n_)
    colors_out = [linspace(color1_(1),color2_(1),n_); linspace(color1_(2),color2_(2),n_); linspace(color1_(3),color2_(3),n_)]';
end

function new_lims = update_lims(old_lims, x_, y_)
    new_lims = old_lims;
    new_lims(1) = min([old_lims(1) min(x_)]);
    new_lims(2) = max([old_lims(2) max(x_)]);
    new_lims(3) = min([old_lims(3) min(y_)]);
    new_lims(4) = max([old_lims(4) max(y_)]);
end
