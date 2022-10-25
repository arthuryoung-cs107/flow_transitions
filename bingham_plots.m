classdef bingham_plots < main_plots
    properties

    end
    methods
        function obj = bingham_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@main_plots(write_figs_, write_all_figs_, figs_to_write_);
        end
        function fig_out = FB1_FB2_Bingham_fluid_params(obj,AYfig_,FB1,FB2)
            [tdim1,tdim2] = deal(2,2);
            axs=prep_tiles(AYfig_,[tdim1 tdim2]);
            axi=0;

            for FB = [FB1 FB2]
                explen=length(FB.exp);
                [tauy_vec mup_vec q_vec] = deal(nan(explen,1));
                color_mat=nan(explen,3);
                for i=1:(explen)
                    expi=FB.exp(i);

                    [tauy_vec(i) mup_vec(i)] = deal(expi.tau_y_Bingham, expi.mu_p_Bingham);
                    [color_mat(i,:) q_vec(i)] = deal(FB.exp(i).color, FB.exp(i).q);
                end
                scatter(axs(1+axi), q_vec, tauy_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                scatter(axs(2+axi), q_vec, mup_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                axi=axi+tdim2;
            end
            ylabel([axs(1) axs(1+tdim2)], '$$\tau_y$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(2) axs(2+tdim2)], '$$\mu_p$$', 'Interpreter', 'LaTeX','FontSize',14)

            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            set(axs,'YScale','log')

            fig_out=AYfig_;
        end
        function fig_out = FB_tau_vs_omega_linscale(obj,AYfig_,FB,tile_dims_)
            axs=prep_tiles(AYfig_,tile_dims_);
            for i=1:length(FB.exp)
                % plot(axs(i), FB.exp(i).omega, FB.exp(i).tau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 2*FB.MS, 'DisplayName', FB.exp(i).label);

                inds = FB.exp(i).omega<40;
                plot(axs(i), FB.exp(i).omega(inds), FB.exp(i).tau(inds), FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 2*FB.MS, 'DisplayName', FB.exp(i).label);

                title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
            end
            % set(axs,'YScale', 'log', 'XScale', 'log');
            ylabel(axs, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)
            fig_out = AYfig_;
        end
        function fig_out = FB_tau_vs_omega_tauy_fit_range(obj,AYfig_,FB,tile_dims_)
            axs=prep_tiles(AYfig_,tile_dims_);
            xinds = 1:10;
            for i=1:length(FB.exp)
                plot(axs(i), FB.exp(i).omega(xinds), FB.exp(i).tau(xinds), FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 2*FB.MS, 'DisplayName', FB.exp(i).label);
                title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
            end
            ylabel(axs, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)
            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_T_vs_rpm_tau_vs_omega(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[2,3]);

            FB=FB1;
            for i=1:length(FB.exp)
              plot(axs(1), FB.exp(i).mu_rpm, FB.exp(i).mu_torque, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end
            for i=1:length(FB.exp)
              plot(axs(2), FB.exp(i).omega, FB.exp(i).tau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
              plot(axs(2), 9.5e-3, FB.exp(i).tau_y,'p','Color',FB.exp(i).color,'LineWidth',2)
            end
            for i=1:length(FB.exp)
              plot(axs(3), FB.exp(i).q, FB.exp(i).tau_y, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 3*FB.MS, 'DisplayName', FB.exp(i).label);
            end

            FB=FB2;
            for i=1:length(FB.exp)
              plot(axs(4), FB.exp(i).mu_rpm, FB.exp(i).mu_torque, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end
            for i=1:length(FB.exp)
              plot(axs(5), FB.exp(i).omega, FB.exp(i).tau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
              plot(axs(5), 9.5e-3, FB.exp(i).tau_y,'p','Color',FB.exp(i).color,'LineWidth',2)
            end
            for i=1:length(FB.exp)
              plot(axs(6), FB.exp(i).q, FB.exp(i).tau_y, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 3*FB.MS, 'DisplayName', FB.exp(i).label);
            end

            set(axs,'YScale', 'log', 'XScale', 'log');
            set(axs(3:3:6),'XScale', 'linear');
            % set(axs(3:3:6),'YScale', 'linear', 'XScale', 'linear');

            ylabel(axs(1:3:4), '$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(1:3:4), 'rpm', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(2:3:5), '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(2:3:5), '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3:3:6), '$$\tau_{y}$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(3:3:6), '$$q=\frac{Q}{Q_{inc}}$$', 'Interpreter', 'LaTeX','FontSize',12)

            title(axs(1), 'FB1: $$rpm$$ vs. $$T_{avg}$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(2), 'FB1: $$\omega_i$$ vs. $$\tau_i$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(3), 'FB1: $$q$$ vs. $$\tau_y$$', 'Interpreter', 'Latex', 'Fontsize', 14);

            title(axs(4), 'FB2: $$rpm$$ vs. $$T_{avg}$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(5), 'FB2: $$\omega_i$$ vs. $$\tau_i$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(6), 'FB2: $$q$$ vs. $$\tau_y$$', 'Interpreter', 'Latex', 'Fontsize', 14);

            ylim(axs(1:3:4),obj.torque_range);
            ylim(axs(2:3:5),[1e-3,1e3]);

            xlim(axs(1:3:4),[1e-1 1e4]);
            xlim(axs(2:3:5),[9e-3 2e2]);

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_bingham_visuals(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[3,4]);
            omegai_logspace = logspace(-2,2,100);
            omegai_plot = reshape(omegai_logspace,[],1);

            FB1_gammai_ro_mat = nan(length(omegai_logspace),length(FB1.exp));
            FB1_gammai_rc_mat = nan(length(omegai_logspace),length(FB1.exp));
            FB1_gammai_rx_mat = nan(length(omegai_logspace),length(FB1.exp));
            FB1_taui_ro_mat = nan(length(omegai_logspace),length(FB1.exp));
            FB1_taui_rc_mat = nan(length(omegai_logspace),length(FB1.exp));
            FB1_taui_rx_mat = nan(length(omegai_logspace),length(FB1.exp));
            FB1_plug_ogt = nan(3,length(FB1.exp));

            FB2_gammai_ro_mat = nan(length(omegai_logspace),length(FB2.exp));
            FB2_gammai_rc_mat = nan(length(omegai_logspace),length(FB2.exp));
            FB2_gammai_rx_mat = nan(length(omegai_logspace),length(FB2.exp));
            FB2_taui_ro_mat = nan(length(omegai_logspace),length(FB2.exp));
            FB2_taui_rc_mat = nan(length(omegai_logspace),length(FB2.exp));
            FB2_taui_rx_mat = nan(length(omegai_logspace),length(FB2.exp));
            FB2_plug_ogt = nan(3,length(FB2.exp));

            for i=1:length(FB1.exp)
                [gammai_ro,taui_ro] = FB1.exp(i).comp_gammai_taui_ro(omegai_logspace);
                [gammai_rc,taui_rc] = FB1.exp(i).comp_gammai_taui_rc(omegai_logspace);
                [gammai_rx,taui_rx,ogt] = FB1.exp(i).comp_gammai_taui(omegai_logspace);

                FB1_gammai_ro_mat(:,i)=gammai_ro;
                FB1_gammai_rc_mat(:,i)=gammai_rc;
                FB1_gammai_rx_mat(:,i)=gammai_rx;
                FB1_taui_ro_mat(:,i)=taui_ro;
                FB1_taui_rc_mat(:,i)=taui_rc;
                FB1_taui_rx_mat(:,i)=taui_rx;
                FB1_plug_ogt(:,i)=ogt;
            end
            [FB1_o_plug_sorted,FB1_ogt_order]=sort(FB1_plug_ogt(1,:));
            FB1_g_plug_sorted=FB1_plug_ogt(2,FB1_ogt_order);
            FB1_t_plug_sorted=FB1_plug_ogt(3,FB1_ogt_order);

            for i=1:length(FB2.exp)
                [gammai_ro,taui_ro] = FB2.exp(i).comp_gammai_taui_ro(omegai_logspace);
                [gammai_rc,taui_rc] = FB2.exp(i).comp_gammai_taui_rc(omegai_logspace);
                [gammai_rx,taui_rx,ogt] = FB2.exp(i).comp_gammai_taui(omegai_logspace);

                FB2_gammai_ro_mat(:,i)=gammai_ro;
                FB2_gammai_rc_mat(:,i)=gammai_rc;
                FB2_gammai_rx_mat(:,i)=gammai_rx;
                FB2_taui_ro_mat(:,i)=taui_ro;
                FB2_taui_rc_mat(:,i)=taui_rc;
                FB2_taui_rx_mat(:,i)=taui_rx;
                FB2_plug_ogt(:,i)=ogt;
            end
            [FB2_o_plug_sorted,FB2_ogt_order]=sort(FB2_plug_ogt(1,:));
            FB2_g_plug_sorted=FB2_plug_ogt(2,FB2_ogt_order);
            FB2_t_plug_sorted=FB2_plug_ogt(3,FB2_ogt_order);

            %% gamma plots

            plot_FB_data(axs(1),FB1,omegai_plot,abs(FB1_gammai_ro_mat));
            plot_FB_data(axs(5),FB1,omegai_plot,abs(FB1_gammai_rc_mat));
            plot_FB_data(axs(9),FB1,omegai_plot,abs(FB1_gammai_rx_mat));

            plot_FB_data(axs(3),FB2,omegai_plot,abs(FB2_gammai_ro_mat));
            plot_FB_data(axs(7),FB2,omegai_plot,abs(FB2_gammai_rc_mat));
            plot_FB_data(axs(11),FB2,omegai_plot,abs(FB2_gammai_rx_mat));

            plot_FB_data(axs(2),FB1,omegai_plot,FB1_taui_ro_mat);
            plot_FB_data(axs(6),FB1,omegai_plot,FB1_taui_rc_mat);
            plot_FB_data(axs(10),FB1,omegai_plot,FB1_taui_rx_mat);

            plot_FB_data(axs(4),FB2,omegai_plot,FB2_taui_ro_mat);
            plot_FB_data(axs(8),FB2,omegai_plot,FB2_taui_rc_mat);
            plot_FB_data(axs(12),FB2,omegai_plot,FB2_taui_rx_mat);


            plot(axs(1), FB1_o_plug_sorted, abs(FB1_g_plug_sorted), 'k - p', 'LineWidth',1);
            plot(axs(5), FB1_o_plug_sorted, abs(FB1_g_plug_sorted), 'k - p', 'LineWidth',1);
            plot(axs(9), FB1_o_plug_sorted, abs(FB1_g_plug_sorted), 'k - p', 'LineWidth',1);

            plot(axs(3), FB2_o_plug_sorted, abs(FB2_g_plug_sorted), 'k - p', 'LineWidth',1);
            plot(axs(7), FB2_o_plug_sorted, abs(FB2_g_plug_sorted), 'k - p', 'LineWidth',1);
            plot(axs(11), FB2_o_plug_sorted, abs(FB2_g_plug_sorted), 'k - p', 'LineWidth',1);

            plot(axs(2), FB1_o_plug_sorted, FB1_t_plug_sorted, 'k - p', 'LineWidth',1);
            plot(axs(6), FB1_o_plug_sorted, FB1_t_plug_sorted, 'k - p', 'LineWidth',1);
            plot(axs(10), FB1_o_plug_sorted, FB1_t_plug_sorted, 'k - p', 'LineWidth',1);

            plot(axs(4), FB2_o_plug_sorted, FB2_t_plug_sorted, 'k - p', 'LineWidth',1);
            plot(axs(8), FB2_o_plug_sorted, FB2_t_plug_sorted, 'k - p', 'LineWidth',1);
            plot(axs(12), FB2_o_plug_sorted, FB2_t_plug_sorted, 'k - p', 'LineWidth',1);


            title(axs(1), 'FB1: $$\omega_i$$ vs. $$ \dot{\gamma}_i$$, $$r_c=r_o$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(5), 'FB1: $$\omega_i$$ vs. $$ \dot{\gamma}_i$$, $$r_c=r_i \sqrt{\tau_i/\tau_y}$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(9), 'FB1: $$\omega_i$$ vs. $$ \dot{\gamma}_i$$, $$r_c$$ conditionally evaluated', 'Interpreter', 'Latex', 'Fontsize', 14);

            title(axs(3), 'FB2: $$\omega_i$$ vs. $$ \dot{\gamma}_i$$, $$r_c=r_o$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(7), 'FB2: $$\omega_i$$ vs. $$ \dot{\gamma}_i$$, $$r_c=r_i \sqrt{\tau_i/\tau_y}$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(11), 'FB2: $$\omega_i$$ vs. $$ \dot{\gamma}_i$$, $$r_c$$ conditionally evaluated', 'Interpreter', 'Latex', 'Fontsize', 14);

            title(axs(2), 'FB1: $$\omega_i$$ vs. $$\tau_i$$, $$r_c=r_o$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(6), 'FB1: $$\omega_i$$ vs. $$\tau_i$$, $$r_c=r_i \sqrt{\tau_i/\tau_y}$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(10), 'FB1: $$\omega_i$$ vs. $$\tau_i$$, $$r_c$$ conditionally evaluated', 'Interpreter', 'Latex', 'Fontsize', 14);

            title(axs(4), 'FB2: $$\omega_i$$ vs. $$\tau_i$$, $$r_c=r_o$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(8), 'FB2: $$\omega_i$$ vs. $$\tau_i$$, $$r_c=r_i \sqrt{\tau_i/\tau_y}$$', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(12), 'FB2: $$\omega_i$$ vs. $$\tau_i$$, $$r_c$$ conditionally evaluated', 'Interpreter', 'Latex', 'Fontsize', 14);


            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:4:9), '$$\dot{\gamma}_{i}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)
            ylim(axs(1:4:9),[1e-3 1e3]);

            ylabel(axs(3:4:11), '$$\dot{\gamma}_{i}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)
            ylim(axs(3:4:11),[1e-3 1e3]);

            ylabel(axs(2:4:10), '$$\tau_{i}$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            ylim(axs(2:4:10),[1e-3 1e3]);

            ylabel(axs(4:4:12), '$$\tau_{i}$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            ylim(axs(4:4:12),[1e-3 1e3]);

            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_radial_shear(obj,AYfig_,FB1,FB2)
            omega_count=20;
            r_count=50;

            run AYfigprops.m
            axs=prep_tiles(AYfig_,[2,2]);

            omegai_logspace = logspace(-2,2,omega_count);

            i = 1;

            FB=FB1;
            [r_it o_it g_it t_it]=FB.exp(i).comp_radial_shear(omegai_logspace,r_count);
            r_plot=r_it(1:(r_count-1),:);
            o_plot=o_it(1:(r_count-1),:);
            t_plot=t_it(1:(r_count-1),:);
            g_plot=abs(g_it(1:(r_count-1),:));

            plot3(axs(1),r_plot,o_plot,g_plot,['- ' FB.specs],'Color',FB.exp(i).color,'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            xlabel(axs(1), '$$r$$ [m]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(1), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
            zlabel(axs(1), '$$\dot{\gamma}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)
            set(axs(1),'ZScale', 'log', 'YScale', 'log');
            view(axs(1),view_mat(1,:));

            plot3(axs(2),r_plot,g_plot,t_plot,['- ' FB.specs],'Color',FB.exp(i).color,'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            xlabel(axs(2), '$$r$$ [m]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(2), '$$\dot{\gamma}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)
            zlabel(axs(2), '$$\tau$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            set(axs(2),'ZScale', 'log', 'YScale', 'log');
            view(axs(2),view_mat(1,:));

            FB=FB2;
            [r_it o_it g_it t_it]=FB.exp(i).comp_radial_shear(omegai_logspace,r_count);
            r_plot=r_it(1:(r_count-1),:);
            o_plot=o_it(1:(r_count-1),:);
            t_plot=t_it(1:(r_count-1),:);
            g_plot=abs(g_it(1:(r_count-1),:));

            plot3(axs(3),r_plot,o_plot,g_plot,['- ' FB.specs],'Color',FB.exp(i).color,'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            xlabel(axs(3), '$$r$$ [m]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(3), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
            zlabel(axs(3), '$$\dot{\gamma}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)
            set(axs(3),'ZScale', 'log', 'YScale', 'log');
            view(axs(3),view_mat(1,:));

            plot3(axs(4),r_plot,g_plot,t_plot,['- ' FB.specs],'Color',FB.exp(i).color,'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            xlabel(axs(4), '$$r$$ [m]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs(4), '$$\dot{\gamma}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)
            zlabel(axs(4), '$$\tau$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            set(axs(4),'ZScale', 'log', 'YScale', 'log');
            view(axs(4),view_mat(1,:));

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_analytical_vs_empirical_shear(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[2,3]);
            omegai_logspace = logspace(-2,2,100);

            FB=FB1;
            for i=1:length(FB.exp)
                plot(axs(1), FB.exp(i).omega, FB.exp(i).tau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end
            FB=FB2;
            for i=1:length(FB.exp)
                plot(axs(4), FB.exp(i).omega, FB.exp(i).tau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end

            FB=FB1;
            for i=1:length(FB.exp)
                [gamma_plot,tau_plot] = FB.exp(i).comp_gammai_taui(omegai_logspace);
                plot(axs(2), omegai_logspace, tau_plot, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end
            FB=FB2;
            for i=1:length(FB.exp)
                [gamma_plot,tau_plot] = FB.exp(i).comp_gammai_taui(omegai_logspace);
                plot(axs(5), omegai_logspace, tau_plot, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end

            FB=FB1;
            for i=1:length(FB.exp)
                [gamma_plot,tau_plot] = FB.exp(i).comp_gammai_taui(FB.exp(i).omega);
                plot(axs(3), FB.exp(i).omega, FB.exp(i).tau./tau_plot, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end
            FB=FB2;
            for i=1:length(FB.exp)
                [gamma_plot,tau_plot] = FB.exp(i).comp_gammai_taui(FB.exp(i).omega);
                plot(axs(6), FB.exp(i).omega, FB.exp(i).tau./tau_plot, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            end

            title(axs(1), 'FB1: $$\omega_i$$ vs. $$\tau_i$$ (measured)', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(2), 'FB1: $$\omega_i$$ vs. $$\hat{\tau}_i$$ (analytical)', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(3), 'FB1: $$\omega_i$$ vs. $$\tau_i/\hat{\tau}_i$$', 'Interpreter', 'Latex', 'Fontsize', 14);

            title(axs(4), 'FB2: $$\omega_i$$ vs. $$\tau_i$$ (measured)', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(5), 'FB2: $$\omega_i$$ vs. $$\hat{\tau}_i$$ (analytical)', 'Interpreter', 'Latex', 'Fontsize', 14);
            title(axs(6), 'FB2: $$\omega_i$$ vs. $$\tau_i/\hat{\tau}_i$$', 'Interpreter', 'Latex', 'Fontsize', 14);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel([axs(1), axs(2), axs(4), axs(5)], '$$\tau_i$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            ylabel(axs(3:3:6), '$$G_{rat} = \tau_i/\hat{\tau}_i$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

            ylim([axs(1), axs(2), axs(4), axs(5)],obj.tau_range);
            % ylim(axs(3:3:6),obj.Grat_range);
            ylim(axs(3:3:6),[3e-1 1e1]);
            xlim(axs,obj.omega_range);

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Bingham_fluid_fits_full(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    omega_low = min(FB.exp(i).omega);
                    omega_cap = max(FB.exp(i).omega);
                    omega_linspace=linspace(omega_low,omega_cap,100);

                    [gi,ti] = FB.exp(i).comp_gammai_taui(omega_linspace);

                    plot(axs(i), FB.exp(i).omega, FB.exp(i).tau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), omega_linspace, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            ylabel(axs_full, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full, 'YScale', 'log', 'XScale', 'log');

            fig_out = AYfig_;
        end

    end
end

function plot_FB_data(ax_,FB_,x_,Y_)
    for i=1:length(FB_.exp)
        plot(ax_, x_,Y_(:,i), FB_.specs, 'Color', FB_.exp(i).color, 'LineWidth', FB_.LW, 'MarkerSize', FB_.MS, 'DisplayName', FB_.exp(i).label);
    end
end

function ax_out = prep_tiles(AYfig_,dims_)
    AYfig_.init_tiles(dims_);
    ax_out=AYfig_.ax_tile;
    hold(ax_out,'on');
    box(ax_out, 'on');
end
