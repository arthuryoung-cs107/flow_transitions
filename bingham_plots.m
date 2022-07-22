classdef bingham_plots < main_plots
    properties

    end
    methods
        function obj = bingham_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@main_plots(write_figs_, write_all_figs_, figs_to_write_);
        end
        function fig_out = FB1_FB2_T_vs_omega(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[1,2]);

            for i=1:length(FB1.exp)
              plot(axs(1), FB1.exp(i).omega, FB1.exp(i).mu_torque, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
              plot(axs(2), FB2.exp(i).omega, FB2.exp(i).mu_torque, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs, '$$T_{avg}$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            axis(axs,[obj.omega_range obj.torque_range])

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_gammai_vs_omega(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[1,2]);

            omegai_logspace = logspace(-2,2,100);

            for i=1:length(FB1.exp)
                plot(axs(1), omegai_logspace, FB1.exp(i).comp_gammai_ro(omegai_logspace), FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            for i=1:length(FB2.exp)
                plot(axs(2), omegai_logspace, FB2.exp(i).comp_gammai_ro(omegai_logspace), FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs, '$$\dot{\gamma}_{i}$$ [1/s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            % axis(axs,[obj.omega_range obj.torque_range])

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_bingham_visuals(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[2,3]);

            fig_out=AYfig_;
        end
    end
end

function ax_out = prep_tiles(AYfig_,dims_)
    AYfig_.init_tiles(dims_);
    ax_out=AYfig_.ax_tile;
    hold(ax_out,'on');
    box(ax_out, 'on');
end
