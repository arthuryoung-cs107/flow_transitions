classdef FB_transition_plots < main_plots
    properties

    end
    methods
        function obj = FB_transition_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@main_plots(write_figs_, write_all_figs_, figs_to_write_);
        end
        function fig_out = FB1_FB2_tau_vs_gamma(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([4,6]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');
            len1 = length(FB1.exp);

            pLW = 2*FB1.LW;
            pMS = 2*FB1.MS;

            for i=1:length(FB1.exp)
                [gammai,taui] = FB1.exp(i).gamma_tau_analytical([-2 2]);
                plot(axs(i), gammai, taui, '-', 'Color', FB1.exp(i).color, 'LineWidth', 2, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
              % legend(axs(i), 'Show', 'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            end
            for i=1:length(FB2.exp)
                plot(axs(i+len1), FB2.exp(i).gamma_analytical, FB2.exp(i).tau_analytical, '-', 'Color', FB2.exp(i).color, 'LineWidth', 2, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
              % legend(axs(i+len1), 'Show', 'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            end

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:6:19), '$$\tau$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(19:24), '$$\dot{\gamma}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)

            % axis(axs,obj.Res_G_range)

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_alpha_comp_vs_shear(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([2,3]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            pLW = 2*FB1.LW;
            pMS = 2*FB1.MS;

            for i=1:length(FB1.exp)
              legend_set_a(i) = plot(axs(1), FB1.exp(i).omega, FB1.exp(i).alpha_T, [FB1.specs '-'], 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB1.exp)
              legend_set_b(i) = plot(axs(2), FB1.exp(i).Re_s, FB1.exp(i).alpha_cf, [FB1.specs '-'], 'Color', FB1.exp(i).color,'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB1.exp)
              legend_set_c(i) = plot(axs(3), FB1.exp(i).Re_s, FB1.exp(i).alpha_G, [FB1.specs '-'], 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
              legend_set_d(i) = plot(axs(4), FB2.exp(i).omega, FB2.exp(i).alpha_T, [FB2.specs '-'], 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end
            for i=1:length(FB2.exp)
              legend_set_e(i) = plot(axs(5), FB2.exp(i).Re_s, FB2.exp(i).alpha_cf, [FB2.specs '-'], 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end
            for i=1:length(FB2.exp)
              legend_set_f(i) = plot(axs(6), FB2.exp(i).Re_s, FB2.exp(i).alpha_G, [FB2.specs '-'], 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            set(axs,'XScale', 'log');

            ylim(axs, [-0.5 obj.alpha_range(2)])
            xlim(axs(1:3:4), obj.omega_range)
            xlim([axs(2:3) axs(5:end)], obj.Re_s_range)

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_alphaT_vs_omega(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([4,6]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');
            len1 = length(FB1.exp);

            pLW = 2*FB1.LW;
            pMS = 2*FB1.MS;

            for i=1:length(FB1.exp)
                alpha = FB1.exp(i).alpha_T;
                plot(axs(i), FB1.exp(i).omega, alpha, [FB1.specs '-'], 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
                plot(axs(i), FB1.exp(i).omega_TV_Ta, alpha(FB1.exp(i).TV_range_Ta), FB1.specs, 'Color', [0 0 0], 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', [FB1.exp(i).label ', transitioned']);
                plot(axs(i), FB1.exp(i).omega_c3_Ta, FB1.exp(i).alpha_tol, '|', 'Color', [0 0 0], 'LineWidth', pLW, 'MarkerSize', 2*pMS, 'DisplayName', [FB1.exp(i).label '$$, \omega_{i,c,3}$$']);
                legend(axs(i), 'Show', 'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            end
            for i=1:length(FB2.exp)
                alpha = FB2.exp(i).alpha_T;
                plot(axs(i+len1), FB2.exp(i).omega, FB2.exp(i).alpha_T, [FB2.specs '-'], 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
                plot(axs(i+len1), FB2.exp(i).omega_TV_Ta, alpha(FB2.exp(i).TV_range_Ta), FB2.specs, 'Color', [0 0 0], 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', [FB2.exp(i).label ', transitioned']);
                plot(axs(i+len1), FB2.exp(i).omega_c3_Ta, FB2.exp(i).alpha_tol, '|', 'Color', [0 0 0], 'LineWidth', pLW, 'MarkerSize', 2*pMS, 'DisplayName', [FB2.exp(i).label '$$, \omega_{i,c,3}$$']);
                legend(axs(i+len1), 'Show', 'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            end

            set(axs, 'XScale', 'log');

            ylabel(axs(1:6:19), '$$\alpha_{T}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(19:24), '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            axis(axs, [obj.omega_range -0.5 obj.alpha_range(2)])

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_T_vs_omega_trans(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([4,6]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');
            len1 = length(FB1.exp);

            pLW = 2*FB1.LW;
            pMS = 2*FB1.MS;

            for i=1:length(FB1.exp)
              plot(axs(i), FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
              plot(axs(i), FB1.exp(i).omega_TV_Ta, FB1.exp(i).T_TV_Ta, FB1.specs, 'Color',[0 0 0], 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
              plot(axs(i), obj.omega_range, FB1.exp(i).tau_c*obj.tau2T*ones(1,2), '-', 'Color',[0 0 0], 'LineWidth', 1, 'MarkerSize', 1, 'DisplayName', [FB1.exp(i).label ' T_c']);
            end
            for i=1:length(FB2.exp)
              plot(axs(i+len1), FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
              plot(axs(i+len1), FB2.exp(i).omega_TV_Ta, FB2.exp(i).T_TV_Ta, FB2.specs, 'Color',[0 0 0], 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
              plot(axs(i+len1), obj.omega_range, FB2.exp(i).tau_c*obj.tau2T*ones(1,2), '-', 'Color',[0 0 0], 'LineWidth', 1, 'MarkerSize', 1, 'DisplayName', [FB2.exp(i).label ' T_c']);
            end

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:6:19), '$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(19:24), '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            axis(axs,[obj.omega_range obj.torque_range])

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_G_vs_Res_trans(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([4,6]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');
            len1 = length(FB1.exp);

            pLW = 2*FB1.LW;
            pMS = 2*FB1.MS;

            for i = 1:length(axs)
                fplot(axs(i), @(Re) obj.G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 1, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
            end
            for i=1:length(FB1.exp)
                plot(axs(i), FB1.exp(i).Re_s, FB1.exp(i).G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
                plot(axs(i), FB1.exp(i).Re_s_TV, FB1.exp(i).G_TV, FB1.specs, 'Color',[0 0 0], 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', [FB1.exp(i).label ', transitioned']);
                plot(axs(i), FB1.exp(i).Re_sc2, FB1.exp(i).G_c2, 'p', 'Color',[0 0 0], 'LineWidth', pLW, 'MarkerSize', pMS, 'DisplayName', [FB1.exp(i).label '$$, Re_{s,c,2}$$']);
                plot(axs(i), FB1.exp(i).Re_sc3, FB1.exp(i).G_c3, '|', 'Color',[0 0 0], 'LineWidth', pLW, 'MarkerSize', 2*pMS, 'DisplayName', [FB1.exp(i).label '$$, Re_{s,c,3}$$']);
                fplot(axs(i), @(Re) (FB1.exp(i).powerfit.b).*(Re).^(FB1.exp(i).powerfit.m), [71 10000],'-', 'Color', FB1.exp(i).color,'Linewidth', 1, 'DisplayName', [FB1.exp(i).label ' $$\beta Re_s^{\alpha}$$']);
                legend(axs(i), 'Show', 'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            end
            for i=1:length(FB2.exp)
                plot(axs(i+len1), FB2.exp(i).Re_s, FB2.exp(i).G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
                plot(axs(i+len1), FB2.exp(i).Re_s_TV, FB2.exp(i).G_TV, FB2.specs, 'Color',[0 0 0], 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', [FB2.exp(i).label ', transitioned']);
                plot(axs(i+len1), FB2.exp(i).Re_sc2, FB2.exp(i).G_c2, 'p', 'Color',[0 0 0], 'LineWidth', pLW, 'MarkerSize', pMS, 'DisplayName', [FB2.exp(i).label '$$, Re_{s,c,2}$$']);
                plot(axs(i+len1), FB2.exp(i).Re_sc3, FB2.exp(i).G_c3, '|', 'Color',[0 0 0], 'LineWidth', pLW, 'MarkerSize', 2*pMS, 'DisplayName', [FB2.exp(i).label '$$, Re_{s,c,3}$$']);
                fplot(axs(i+len1), @(Re) (FB2.exp(i).powerfit.b).*(Re).^(FB2.exp(i).powerfit.m), [71 10000],'-', 'Color', FB2.exp(i).color,'Linewidth', 1, 'DisplayName', [FB2.exp(i).label ' $$\beta Re_s^{\alpha}$$']);
                legend(axs(i+len1), 'Show', 'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            end
            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:6:19), '$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(19:24), '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            axis(axs,obj.Res_G_range)

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Re_sc_vs_q(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([2,3]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i=1:length(FB1.exp)
              plot(axs(1), FB1.exp(i).q, FB1.exp(i).Re_sc1, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', 2*FB1.LW_L, 'MarkerSize', 2*FB1.MS_L, 'DisplayName', FB1.exp(i).label);
              plot(axs(2), FB1.exp(i).q, FB1.exp(i).Re_sc2, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', 2*FB1.LW_L, 'MarkerSize', 2*FB1.MS_L, 'DisplayName', FB1.exp(i).label);
              plot(axs(3), FB1.exp(i).q, FB1.exp(i).Re_sc3, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', 2*FB1.LW_L, 'MarkerSize', 2*FB1.MS_L, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
              plot(axs(4), FB2.exp(i).q, FB2.exp(i).Re_sc1, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', 2*FB2.LW_L, 'MarkerSize', 2*FB2.MS_L, 'DisplayName', FB2.exp(i).label);
              plot(axs(5), FB2.exp(i).q, FB2.exp(i).Re_sc2, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', 2*FB2.LW_L, 'MarkerSize', 2*FB2.MS_L, 'DisplayName', FB2.exp(i).label);
              plot(axs(6), FB2.exp(i).q, FB2.exp(i).Re_sc3, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', 2*FB2.LW_L, 'MarkerSize', 2*FB2.MS_L, 'DisplayName', FB2.exp(i).label);
            end

            ylabel(axs(1:3:4), '$$Re_{s,c1}$$ ', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(2:3:5), '$$Re_{s,c2}$$ ', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3:3:6), '$$Re_{s,c3}$$ ', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$q$$', 'Interpreter', 'LaTeX','FontSize',12)

            ylim(axs, [0,120])

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_T_alphaT_vs_omega(obj, AYfig_, FB1, FB2)
            AYfig_.init_tiles([2,2]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            pLW = 2*FB1.LW;
            pMS = 2*FB1.MS;

            for i=1:length(FB1.exp)
              legend_set_a(i) = plot(axs(1), FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
              plot(axs(1), FB1.exp(i).omega_c3_Ta, FB1.exp(i).T_c3_Ta, '|', 'Color', [0 0 0], 'LineWidth', pLW, 'MarkerSize', pMS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
              legend_set_b(i) = plot(axs(2), FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
              plot(axs(2), FB2.exp(i).omega_c3_Ta, FB2.exp(i).T_c3_Ta, '|', 'Color', [0 0 0],'LineWidth', pLW, 'MarkerSize', pMS, 'DisplayName', FB2.exp(i).label);
            end
            for i=1:length(FB1.exp)
              legend_set_c(i) = plot(axs(3), FB1.exp(i).omega, FB1.exp(i).alpha_T, [FB1.specs '-'], 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
              plot(axs(3), FB1.exp(i).omega_c3_Ta, FB1.exp(i).alpha_tol, '|', 'Color', [0 0 0], 'LineWidth', pLW, 'MarkerSize', 2*pMS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
              legend_set_d(i) = plot(axs(4), FB2.exp(i).omega, FB2.exp(i).alpha_T, [FB2.specs '-'], 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
              plot(axs(4), FB2.exp(i).omega_c3_Ta, FB2.exp(i).alpha_tol, '|', 'Color', [0 0 0], 'LineWidth', pLW, 'MarkerSize', 2*pMS, 'DisplayName', FB2.exp(i).label);
            end

            set(axs,'XScale', 'log');
            set(axs(1:2),'YScale', 'log');

            ylabel(axs(1), '$$T$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3), '$$\alpha_T$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:4), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns',2);
            legend(axs(2), legend_set_b,'Location', 'SouthEast', 'Interpreter', 'Latex', 'NumColumns',2);
            legend(axs(3), legend_set_c,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',2);
            legend(axs(4), legend_set_d,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',2);

            axis(axs(1:2), [obj.omega_range obj.torque_range])
            axis(axs(3:4), [obj.omega_range -0.5 obj.alpha_range(2)])

            fig_out = AYfig_;
        end
    end
end
