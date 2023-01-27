classdef FB_raw_data_plots
    properties (Constant)
        r_i = 0.01208;
        r_o = 0.025;
        L = 0.036;
        tau2T = 2*pi*(0.01208)*(0.01208)*(0.036);

        posdim_full = [0 0 1728 1000];
        posdim_toprow = [1 551 1728 460];
        posdim_bottomrow = [0 1 1728 460];
    end
    properties
        write_figs;
        write_all_figs;
        figs_to_write;
    end
    methods
        function obj = FB_raw_data_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj.write_figs = write_figs_;
            obj.write_all_figs = write_all_figs_;
            obj.figs_to_write = figs_to_write_;
        end
        function fig_out_ = FB1_torque_vs_rpm_raw_plot(obj,AYfig_,FB)
            axs=prep_tiles(AYfig_,[4,5]);

            for i = 1:length(FB.exp)
                exp = FB.exp(i);
                plot(axs(i), exp.nf_raw, exp.Tf_raw, exp.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                plot(axs(i), exp.mu_rpm, exp.mu_torque, '- p', 'Color', [0 0 0], 'LineWidth', 2*FB.LW, 'MarkerSize', 3*FB.MS, 'DisplayName', FB.exp(i).label);

                plot(axs(i+10), exp.nf_clean, exp.Tf_clean, exp.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                plot(axs(i+10), exp.mu_rpm, exp.mu_torque, '- p', 'Color', [0 0 0], 'LineWidth', 2*FB.LW, 'MarkerSize', 3*FB.MS, 'DisplayName', FB.exp(i).label);
            end
            fig_out_=AYfig_;
            % set(axs, 'YScale', 'log');
            set(axs, 'YScale', 'log', 'XScale', 'log');
        end
    end
    methods (Static)
        function AYfig_out = make_fig(name_,posdim_)
            AYfig_out = AYfig(AYfig.specs_gen(name_,posdim_),false);
        end
    end
end
function ax_out = prep_tiles(AYfig_,dims_)
    AYfig_.init_tiles(dims_);
    ax_out=AYfig_.ax_tile;
    hold(ax_out,'on');
    box(ax_out, 'on');
end
