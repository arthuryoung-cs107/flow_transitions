classdef bingham_plots < main_plots
    properties

    end
    methods
        function obj = bingham_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@main_plots(write_figs_, write_all_figs_, figs_to_write_);
        end
        function fig_out = FBext_Bingham_fits_summary(obj,AYfig_,FBall)
            axs_full = prep_tiles(AYfig_, [4,4]);

            [EC000 EC050 EC075 EC100 FB1 FB2] = deal(FBall{1},FBall{2},FBall{3},FBall{4},FBall{5},FBall{6});

            iEC = 0;
            for EC = [EC000 EC050 EC075 EC100]
                iEC = iEC+1;

                plot(axs_full(iEC), EC.omega, EC.tau_comp, EC.specs, 'Color', EC.color, 'LineWidth', EC.LW, 'MarkerSize', EC.MS, 'DisplayName', EC.label);
                plot(axs_full(iEC+4), EC.omega, EC.tau_comp, EC.specs, 'Color', EC.color, 'LineWidth', EC.LW, 'MarkerSize', EC.MS, 'DisplayName', EC.label);
                plot(axs_full(iEC), EC.omega_fit_Bingham, EC.tau_fit_Bingham, 'o', 'Color', EC.color, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', EC.label);
                plot(axs_full(iEC+4), EC.omega_fit_Bingham, EC.tau_fit_Bingham, 'o', 'Color', EC.color, 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', EC.label);

                [mup ty] = deal(EC.mu_p_Bingham, EC.tau_y_Bingham);
                otest_log = logspace(-2,2,100);
                otest_lin = linspace(min(EC.omega),max(EC.omega),100);
                otest_fit = linspace(min(EC.omega_fit_Bingham),max(EC.omega_fit_Bingham),100);
                ti_log = glass_particles.taui_pred_Bingham(mup,ty,otest_log);
                ti_lin = glass_particles.taui_pred_Bingham(mup,ty,otest_lin);
                ti_fit = glass_particles.taui_pred_Bingham(mup,ty,otest_fit);

                plot(axs_full(iEC), otest_lin, ti_lin, '-', 'Color', EC.color, 'LineWidth', EC.LW, 'MarkerSize', EC.MS, 'DisplayName', EC.label);
                plot(axs_full(iEC+4), otest_log, ti_log, '-', 'Color', EC.color, 'LineWidth', EC.LW, 'MarkerSize', EC.MS, 'DisplayName', EC.label);

                plot(axs_full(9), EC.omega, EC.tau_comp, EC.specs, 'Color', EC.color, 'LineWidth', EC.LW, 'MarkerSize', EC.MS, 'DisplayName', EC.label);
                plot(axs_full(9), otest_log, ti_log, '-', 'Color', EC.color, 'LineWidth', EC.LW, 'MarkerSize', EC.MS, 'DisplayName', EC.label);

                plot(axs_full(13), EC.q, EC.tau_y_Bingham, EC.specs, 'Color', EC.color, 'LineWidth', 2, 'MarkerSize', 9, 'DisplayName', EC.label);
                plot(axs_full(15), EC.q, EC.mu_p_Bingham, EC.specs, 'Color', EC.color, 'LineWidth', 2, 'MarkerSize', 9, 'DisplayName', EC.label);
            end

            FB_orig = [FB1 FB2];

            for iFB = 1:length(FB_orig)
                FB = FB_orig(iFB);
                exp = FB.exp;
                for i = 1:length(exp)
                    expi = exp(i);
                    [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                    [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                    plot(axs_full(8+iFB), expi.omega, expi.tau, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);

                    plot(axs_full(13+iFB-1), expi.q, expi.tau_y_Bingham, FB.specs, 'Color', expi.color, 'LineWidth', 2, 'MarkerSize', 9, 'DisplayName', expi.label);
                    plot(axs_full(15+iFB-1), expi.q, expi.mu_p_Bingham, FB.specs, 'Color', expi.color, 'LineWidth', 2, 'MarkerSize', 9, 'DisplayName', expi.label);
                end
            end

            set(axs_full(5:10), 'XScale', 'log', 'YScale', 'log');

            fig_out = AYfig_;
        end
        function fig_out = UXB_FB_taustar_vs_Gamma_combined(obj,AYfig_,UBall,XBall,FBall,FBext)
            axs_full = prep_tiles(AYfig_, [1,2]);

            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);
            [FB1 FB2] = deal(FBall(1),FBall(2));

            UXB = {UB1; UB2; XB1; XB2};
            UXBlen = length(UXB);
            UXBlabels = cell(UXBlen,1);

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
                    % if (expi.q < expi.qcrit_Bingham2Carreau)
                        [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                        [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));
                        tau_star = tfit/expi.tau_y_Bingham;
                        Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                        Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit./(rfit-fluid.r_i_def)))/expi.tau_y_Bingham;
                        % Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit/(fluid.r_o_def-  fluid.r_i_def)))/expi.tau_y_Bingham;

                        ileg = ileg +1;
                        legend_set_a(ileg) = plot(axs_full(1), Gamma, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        legend_set_b(ileg) = plot(axs_full(2), Gamma_omega, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                    end
                end
            end

            for i = 1:length(FBext)
                expi = FBext{i};
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                tau_star = tfit/expi.tau_y_Bingham;
                Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit./(rfit-fluid.r_i_def)))/expi.tau_y_Bingham;


                ileg = ileg +1;
                legend_set_a(ileg) = plot(axs_full(1), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Gamma_omega, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end

            for i = 1:UXBlen
                expi = UXB{i};
                expALL{i} = expi;
                [mup ty] = deal(expi.mu_p_Bingham, expi.tau_y_Bingham);
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                tau_star = tfit/expi.tau_y_Bingham;
                Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                % Gamma_omega = (expi.mu_eff*fluid.r_i_def*(ofit./(rfit-fluid.r_i_def)))/expi.tau_qs;
                Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit./(rfit-fluid.r_i_def)))/expi.tau_y_Bingham;
                % Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit/(fluid.r_o_def-fluid.r_i_def)))/expi.tau_y_Bingham;

                ileg = ileg + 1;
                legend_set_a(ileg) = plot(axs_full(1), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                legend_set_b(ileg) = plot(axs_full(2), Gamma_omega, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
            end
            set(axs_full,'XScale','log','YScale','log');

            ylabel(axs_full, '$$\tau^* = \tau/\tau_y$$ [dimensionless]', 'Interpreter', 'Latex')
            xlabel(axs_full(1), '$$\Gamma = \dot{\gamma}_B \mu_p/\tau_y$$ [dimensionless]', 'Interpreter', 'Latex')
            xlabel(axs_full(2), '$$\Gamma = [\omega_i r_i / (r_c-r_i) ] \mu_p /\tau_y$$ [dimensionless]', 'Interpreter', 'Latex')

            legend(axs_full(1), legend_set_a,'Location', 'NorthWest', 'Interpreter', 'Latex','NumColumns',2);
            legend(axs_full(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex','NumColumns',2);

            fig_out = AYfig_;
        end
        function fig_out = UXB_FB_taustar_vs_Gamma(obj,AYfig_,UBall,XBall,FBall)
            axs_full = prep_tiles(AYfig_, [4,4]);

            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);
            [FB1 FB2] = deal(FBall(1),FBall(2));

            UXB = {UB1; UB2; XB1; XB2};
            UXBlen = length(UXB);
            UXBlabels = cell(UXBlen,1);

            expALL = cell(4 + length(FB1.exp) + length(FB2.exp),1);
            iexp = 4;
            for iFB = 1:length(FBall)
                FB = FBall(iFB);
                exp = FB.exp;
                for i = 1:length(exp)
                    expi = exp(i);
                    iexp = iexp + 1;
                    expALL{iexp} = expi;
                    [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                    [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                    tau_star = tfit/expi.tau_y_Bingham;
                    Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                    Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit./(rfit-fluid.r_i_def)))/expi.tau_y_Bingham;
                    % Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit/(fluid.r_o_def-fluid.r_i_def)))/expi.tau_y_Bingham;

                    % if (expi.q < expi.qcrit_Bingham2Carreau)
                    if (expi.q < expi.qcrit_sgf)
                        plot(axs_full(2+iFB), Gamma, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        plot(axs_full(2+iFB + 4), Gamma, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);

                        plot(axs_full(9), Gamma, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        plot(axs_full(13), Gamma, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);

                        plot(axs_full(11), Gamma_omega, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                        plot(axs_full(15), Gamma_omega, tau_star, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                    end
                end
            end

            for i = 1:UXBlen
                expi = UXB{i};
                expALL{i} = expi;
                [mup ty] = deal(expi.mu_p_Bingham, expi.tau_y_Bingham);
                [ofit tfit ifit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                [gfit rfit] = deal(expi.gamma_Bingham(ifit), expi.rc_Bingham(ifit));

                tau_star = tfit/expi.tau_y_Bingham;
                Gamma = (expi.mu_p_Bingham*gfit)/expi.tau_y_Bingham;
                % Gamma_omega = (expi.mu_eff*fluid.r_i_def*(ofit./(rfit-fluid.r_i_def)))/expi.tau_qs;
                Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit./(rfit-fluid.r_i_def)))/expi.tau_y_Bingham;
                % Gamma_omega = (expi.mu_p_Bingham*fluid.r_i_def*(ofit/(fluid.r_o_def-fluid.r_i_def)))/expi.tau_y_Bingham;

                plot(axs_full(1+floor((i-1)/2)), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(1+floor((i-1)/2) + 4), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(9), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(10), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(13), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(14), Gamma, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);

                plot(axs_full(11), Gamma_omega, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(12), Gamma_omega, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(15), Gamma_omega, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs_full(16), Gamma_omega, tau_star, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);

            end

            set(axs_full(1:4),'XScale','log','YScale','log');

            set(axs_full(13:end),'XScale','log','YScale','log');

            fig_out = AYfig_;
        end
        function fig_out = UXB_Bingham_fit_summary(obj,AYfig_,UBall,XBall,FBall)
            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);
            [FB1 FB2] = deal(FBall(1),FBall(2));

            UXB = {UB1; UB2; XB1; XB2};
            UXBlen = length(UXB);
            UXBlabels = cell(UXBlen,1);

            expALL = cell(4 + length(FB1.exp) + length(FB2.exp),1);

            axs_full = set_UXB_Bingham_summary_axes(AYfig_);
            axFB = axs_full((length(axs_full)-1):end);
            iexp = 4;
            for FB = FBall
                exp = FB.exp;
                for i = 1:length(exp)
                    expi = exp(i);
                    iexp = iexp + 1;
                    expALL{iexp} = expi;
                    plot(axFB(1), expi.q, expi.tau_y_Bingham, expi.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                    plot(axFB(2), expi.q, expi.mu_p_Bingham, FB.specs, 'Color', expi.color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', expi.label);
                end
            end

            mup_bounds = nan(UXBlen,3);
            tauy_bounds = nan(UXBlen,3);
            for i = 1:UXBlen
                axs = axs_full(i:UXBlen:(length(axs_full)-(UXBlen-i)));

                expi = UXB{i};
                expALL{i} = expi;
                [mup ty] = deal(expi.mu_p_Bingham, expi.tau_y_Bingham);

                otest_log = logspace(-2,2,100);
                otest_lin = linspace(min(expi.omega),max(expi.omega),100);
                otest_fit = linspace(min(expi.omega_fit_Bingham),max(expi.omega_fit_Bingham),100);
                ti_log = glass_particles.taui_pred_Bingham(mup,ty,otest_log);
                ti_lin = glass_particles.taui_pred_Bingham(mup,ty,otest_lin);
                ti_fit = glass_particles.taui_pred_Bingham(mup,ty,otest_fit);

                plot(axs(1), expi.omega, expi.tau_comp, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs(1), otest_log, ti_log, '-', 'Color', [0 0 0], 'LineWidth', 2);

                plot(axs(2), expi.omega, expi.tau_comp, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs(2), otest_lin, ti_lin, '-', 'Color', [0 0 0], 'LineWidth', 2);

                plot(axs(3), expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.specs, 'Color', expi.color, 'LineWidth', expi.LW, 'MarkerSize', expi.MS, 'DisplayName', expi.label);
                plot(axs(3), otest_fit, ti_fit, '-', 'Color', [0 0 0], 'LineWidth', 2);

                plot(axs_full(length(axs_full)-3), i, expi.tau_y_Bingham, 'p', 'Color', expi.color, 'LineWidth', 1.5, 'MarkerSize', 5);
                plot(axs_full(length(axs_full)-3), i, expi.tau_qs, 's', 'Color', expi.color, 'LineWidth', 1.5, 'MarkerSize', 5);
                plot(axs_full(length(axs_full)-2), i, expi.mu_p_Bingham, 'p', 'Color', expi.color, 'LineWidth', 1.5, 'MarkerSize', 5);
                plot(axs_full(length(axs_full)-2), i, expi.mu_eff, 's', 'Color', expi.color, 'LineWidth', 1.5, 'MarkerSize', 5);

                tauy_bounds(i,:) = expi.Bingham_tauy_bounds;
                mup_bounds(i,:) = expi.Bingham_mup_bounds;
            end

            plot(axs_full(length(axs_full)-3), 1:UXBlen, tauy_bounds(:,1), 'v', 'Color', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 5);
            plot(axs_full(length(axs_full)-3), 1:UXBlen, tauy_bounds(:,2), 'o', 'Color', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 5);
            plot(axs_full(length(axs_full)-3), 1:UXBlen, tauy_bounds(:,3), '^', 'Color', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 5);

            plot(axs_full(length(axs_full)-2), 1:UXBlen, mup_bounds(:,1), 'v', 'Color', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 5);
            plot(axs_full(length(axs_full)-2), 1:UXBlen, mup_bounds(:,2), 'o', 'Color', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 5);
            plot(axs_full(length(axs_full)-2), 1:UXBlen, mup_bounds(:,3), '^', 'Color', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 5);

            set(axs_full(1:4),'XScale','log','YScale','log');

            set(axs_full((length(axs_full)-3):(length(axs_full)-2)),'YScale','log');

            fig_out = AYfig_;
        end
        function fig_out = FB_Bingham_err_vs_num_spread(obj,AYfig_,FB1,FB2,FBext,log_flag)
            if (log_flag)
                axscale = 'log';
                % xlim_set = [1e-2 1e2];
                % ylim_set = [1e-4 1e3];
            else
                axscale = 'linear';
                % xlim_set = [0 130];
                % ylim_set = [0 200];
            end

            dims = [4 6];
            axs_full=prep_tiles(AYfig_,dims);

            B_coeff = (fluid.r_o_def-fluid.r_i_def)/(fluid.r_i_def);

            ql_count = 0;
            qh_count = 0;
            Bl_sum = 0;
            Bh_sum = 0;
            Bl_set = [];
            Bh_set = [];

            for i=1:length(FBext)
                expi = FBext{i};

                [omega_fit tau_fit i_fit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham, expi.i_fit_Bingham);
                omega_nfit = expi.omega(~i_fit);

                B_fit = B_coeff*(expi.tau_y_Bingham/expi.mu_p_Bingham)./(omega_fit);

                ql_count = ql_count +1;
                Bl_sum = Bl_sum + sum(B_fit);
                Bl_set = [Bl_set; B_fit];
            end

            axs = axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    expi = FB.exp(i);
                    par = [expi.mu_p_Bingham expi.tau_y_Bingham];
                    ti = glass_particles.taui_pred_Bingham(par(1),par(2),expi.omega);
                    tau_err = (abs(expi.tau - ti))./expi.tau;

                    [omega_fit tau_fit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham);
                    [omega_not_fit i_nf] = setdiff(expi.omega, omega_fit);
                    tau_not_fit = expi.tau(i_nf);
                    i_f = logical(ones(size(expi.omega)));
                    i_f(i_nf) = 0;
                    [tau_err_fit tau_err_not_fit] = deal(abs(tau_fit-ti(i_f))./tau_fit,abs(tau_not_fit-ti(i_nf))./tau_not_fit);

                    B_all = B_coeff*(expi.tau_y_Bingham/expi.mu_p_Bingham)./expi.omega;
                    B_fit = B_coeff*(expi.tau_y_Bingham/expi.mu_p_Bingham)./(expi.omega(i_f));
                    B_not_fit = B_coeff*(expi.tau_y_Bingham/expi.mu_p_Bingham)./(expi.omega(i_nf));

                    plot(axs(i), B_all, tau_err, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), B_fit, tau_err_fit, 'o', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), B_not_fit, tau_err_not_fit, 'x', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', FB.exp(i).label);

                    if (expi.q<1)
                        ql_count = ql_count +1;
                        Bl_sum = Bl_sum + sum(B_fit);
                        Bl_set = [Bl_set; B_fit];
                    else
                        qh_count = qh_count +1;
                        Bh_sum = Bh_sum + sum(B_fit);
                        Bh_set = [Bh_set; B_fit];
                    end

                    title(axs(i), expi.label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+1:end);
            end

            ylabel(axs_full(1:dims(2):22), '$$| \tau_i - \tau_B |/\tau_i $$', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full(22:end), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full,'XScale',axscale,'YScale',axscale);

            fprintf('Mean q<1 fitted Bingham number: %e, median :%e. Mean q>1 fitted Bingham number: %e, median: %e. \n', Bl_sum/ql_count, median(Bl_set),  Bh_sum/qh_count, median(Bh_set));

            fig_out = AYfig_;
        end
        function fig_out = FB_taui_rat_Bingham_vs_omegai_spread(obj,AYfig_,FB1,FB2,log_flag_)
            if (nargin==4)
                log_flag = true;
            else
                log_flag = log_flag_;
            end
            if (log_flag)
                axscale = 'log';
                xlim_set = [1e-2 1e2];
                ylim_set = [1e-4 1e3];
            else
                axscale = 'linear';
                xlim_set = [0 130];
                ylim_set = [0 200];
            end

            dims = [4 6];
            axs_full=prep_tiles(AYfig_,dims);

            axs = axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    expi = FB.exp(i);
                    par = [expi.mu_p_Bingham expi.tau_y_Bingham];
                    ti = glass_particles.taui_pred_Bingham(par(1),par(2),expi.omega);
                    tau_rat = expi.tau./ti;

                    [omega_fit tau_fit] = deal(expi.omega_fit_Bingham, expi.tau_fit_Bingham);
                    [omega_not_fit i_nf] = setdiff(expi.omega, omega_fit);
                    tau_not_fit = expi.tau(i_nf);
                    i_f = logical(ones(size(expi.omega)));
                    i_f(i_nf) = 0;
                    [tau_rat_fit tau_rat_not_fit] = deal(tau_fit./ti(i_f),tau_not_fit./ti(i_nf));

                    [oc tc] = glass_particles.get_plug_crit(par(1),par(2));
                    plot(axs(i), expi.omega, tau_rat, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), omega_fit, tau_rat_fit, 'o', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), omega_not_fit, tau_rat_not_fit, 'x', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', FB.exp(i).label);

                    plot(axs(i), [min(expi.omega) oc oc], [1 1 min(tau_rat)], '--', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 1);

                    title(axs(i), expi.label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+1:end);
            end

            ylabel(axs_full(1:dims(2):22), '$$\tau_i$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full(22:end), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full,'XScale',axscale,'YScale',axscale);

            fig_out = AYfig_;
        end
        function fig_out = ALL_tau_vs_omegai_spread(obj,AYfig_,FBall,PFall,NBall,UBall,XBall,log_flag_)
            if (nargin==7)
                log_flag = true;
            else
                log_flag = log_flag_;
            end
            if (log_flag)
                axscale = 'log';
                xlim_set = [1e-2 1e2];
                ylim_set = [1e-4 1e3];
            else
                axscale = 'linear';
                % xlim_set = [0 130];
                % ylim_set = [0 250];
                xlim_set = [0 130];
                ylim_set = [0 200];
            end

            [PF1 PFR] = deal(PFall{1}, PFall{2});
            [NB1 NB2 NB3] = deal(NBall.NB1_in, NBall.NB2_in, NBall.NB3_in);
            [UB1 UB2] = deal(UBall.UB1_in, UBall.UB2_in);
            [XB1 XB2] = deal(XBall.XB1_in, XBall.XB2_in);

            dims = [4 7];
            axs_full=prep_tiles(AYfig_,dims);

            axs = axs_full;
            for FB = FBall
                explen = length(FB.exp);
                for i=1:explen
                    expi = FB.exp(i);

                    plot(axs(i), expi.omega, expi.tau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    title(axs(i), expi.label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+1:end);
            end

            plot(axs_full(24), PF1.omega, PF1.tau_comp, PF1.specs, 'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            plot(axs_full(24), PFR.omega, PFR.tau_comp, PFR.specs, 'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            plot(axs_full(25), NB1.omega, NB1.tau_comp, NB1.specs, 'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            plot(axs_full(25), NB2.omega, NB2.tau_comp, NB2.specs, 'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            plot(axs_full(25), NB3.omega, NB3.tau_comp, NB3.specs, 'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            plot(axs_full(26), UB1.omega, UB1.tau_comp, UB1.specs, 'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            plot(axs_full(26), UB2.omega, UB2.tau_comp, UB2.specs, 'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            plot(axs_full(27), XB1.omega, XB1.tau_comp, XB1.specs, 'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            plot(axs_full(27), XB2.omega, XB2.tau_comp, XB2.specs, 'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            ylabel(axs_full(1:dims(2):22), '$$\tau_i$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full(22:end), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full,'XScale',axscale,'YScale',axscale);
            xlim(axs_full,xlim_set);
            % ylim(axs_full,ylim_set);

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Bingham_fluid_params_fitcheck(obj,AYfig_,FB1,FB2)
            [tdim1,tdim2] = deal(2,2);
            axs=prep_tiles(AYfig_,[tdim1 tdim2]);
            axi=0;

            for FB = [FB1 FB2]
                explen=length(FB.exp);
                [tauy_vec mup_vec q_vec] = deal(nan(explen,1));
                [tauy_vec_old mup_vec_old] = deal(nan(explen,1));
                [mup_bounds_mat tauy_bounds_mat] = deal(nan(explen,3));
                color_mat=nan(explen,3);
                for i=1:(explen)
                    expi=FB.exp(i);

                    [tauy_vec(i) mup_vec(i)] = deal(expi.tau_y_Bingham, expi.mu_p_Bingham);
                    [tauy_vec_old(i) mup_vec_old(i)] = deal(expi.tau_y, expi.mu_p);
                    [mup_bounds_mat(i,:) tauy_bounds_mat(i,:)] = expi.determine_Bingham_fluid_bounds;
                    [color_mat(i,:) q_vec(i)] = deal(FB.exp(i).color, FB.exp(i).q);
                end
                scatter(axs(axi+1), q_vec, tauy_vec, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                scatter(axs(axi+1), q_vec, tauy_bounds_mat(:,1), 'v','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                scatter(axs(axi+1), q_vec, tauy_bounds_mat(:,2), 'p','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                scatter(axs(axi+1), q_vec, tauy_bounds_mat(:,3), '^','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);

                scatter(axs(axi+3), q_vec, mup_vec, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                scatter(axs(axi+3), q_vec, mup_bounds_mat(:,1), 'v','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                scatter(axs(axi+3), q_vec, mup_bounds_mat(:,2), 'p','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                scatter(axs(axi+3), q_vec, mup_bounds_mat(:,3), '^','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);

                axi=axi+1;
            end
            ylabel(axs(1:2), '$$\tau_y$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(3:4), '$$\mu_p$$', 'Interpreter', 'LaTeX','FontSize',14)

            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            % set(axs,'YScale','log')
            % set(axs,'YScale','linear')
            set(axs(1:2),'YScale','log')
            set(axs(3:4),'YScale','log')

            % ylim(axs(1:2), [1e-3 5e2]);
            % ylim(axs(3:4), [1e-4 1e2]);

            xlim(axs(1:2:3), [0 2]);
            xlim(axs(2:2:4), [0 16]);
            % ylim(axs(1:2), [1e-5 1]);
            % ylim(axs(3:4), [0 0.6]);

            fig_out=AYfig_;
        end
        function fig_out = FB1_FB2_Bingham_fits_spread(obj,AYfig_,FB1,FB2,log_flag_)
            if (nargin==4)
                log_flag = false;
            else
                log_flag = log_flag_;
            end

            if (log_flag)
                axscale = 'log';
                omega_test = logspace(-2, 2, 100);
            else
                axscale = 'linear';
            end

            AYfig_.init_tiles([4,6]);
            axs_full = AYfig_.ax_tile;
            hold(axs_full, 'on');
            box(axs_full,'on');

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    par = [FB.exp(i).mu_p_Bingham FB.exp(i).tau_y_Bingham];
                    omega_fit = FB.exp(i).omega_fit_Bingham;
                    tau_fit = FB.exp(i).tau_fit_Bingham;

                    if (log_flag)
                        ti = glass_particles.taui_pred_Bingham(par(1),par(2),omega_test);
                        plot(axs(i), omega_test, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                        plot(axs(i), omega_fit, tau_fit, 'o', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', FB.exp(i).label);

                        [omega_not_fit i_nf] = setdiff(FB.exp(i).omega, omega_fit);
                        tau_not_fit = FB.exp(i).tau(i_nf);

                        plot(axs(i), omega_not_fit, tau_not_fit, 'x', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', FB.exp(i).label);

                        [oc tc] = glass_particles.get_plug_crit(par(1),par(2));
                        plot(axs(i), [min(omega_fit) oc oc], [tc tc min(tau_fit)], '--', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 1);

                    else
                        omega_test = linspace(0, max(omega_fit), 100);
                        ti = glass_particles.taui_pred_Bingham(par(1),par(2),omega_test);
                        plot(axs(i), omega_test, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                        plot(axs(i), omega_fit, tau_fit, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    end

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            ylabel(axs_full, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full, 'YScale', axscale, 'XScale', axscale);

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Bingham_fluid_fits_full(obj,AYfig_,FB1,FB2,full_flag_)
            if (nargin==4)
                full_flag=true;
            else
                full_flag=full_flag_;
            end

            if (full_flag)
                omega_linspace=logspace(-2, 2,100);
                % omega_linspace=linspace(0, 15,100);
                xscale_set = 'log';
                yscale_set = 'log';
                % xscale_set = 'linear';
                % yscale_set = 'linear';
            else
                % omega_linspace=linspace(0, 15,100);
                omega_linspace=logspace(-2, 2,100);
                % omega_linspace=2*logspace(-3, 1,100);
                % xscale_set = 'linear';
                yscale_set = 'linear';
                xscale_set = 'log';
                % yscale_set = 'log';
            end

            AYfig_.init_tiles([4,6]);
            axs_full = AYfig_.ax_tile;
            hold(axs_full, 'on');
            box(axs_full,'on');

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    par = [FB.exp(i).mu_p_Bingham FB.exp(i).tau_y_Bingham];
                    ti = glass_particles.taui_pred_Bingham(par(1),par(2),omega_linspace);

                    par_old = [FB.exp(i).mu_p FB.exp(i).tau_y];
                    % [gi_old,ti_old] = FB.exp(i).comp_gammai_taui(omega_linspace,par_old);
                    ti_old = glass_particles.taui_pred_Bingham(par_old(1),par_old(2),omega_linspace);

                    inds_plot = FB.exp(i).omega <= max(omega_linspace);
                    omega_plot = FB.exp(i).omega(inds_plot);
                    tau_plot = FB.exp(i).tau(inds_plot);

                    plot(axs(i), omega_plot, tau_plot, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), omega_linspace, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                    % plot(axs(i), omega_linspace, ti_old, ':', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);

                    inds_pre_crit = FB.exp(i).omega < FB.exp(i).omega_cap_Bingham_use;
                    omega_fit = FB.exp(i).omega(inds_pre_crit);
                    tau_fit = FB.exp(i).tau(inds_pre_crit);

                    plot(axs(i), omega_fit, tau_fit, FB.specs, 'Color', [0 0 0], 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    [oc tc] = glass_particles.get_plug_crit(par(1),par(2));
                    plot(axs(i), tc*ones(2,1), [1e-2 oc], '--', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 1);
                    plot(axs(i), [1e-2 tc], oc*ones(2,1), '--', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 1);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            ylabel(axs_full, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full, 'YScale', yscale_set, 'XScale', xscale_set);

            axis(axs_full, [1e-2 1e2 1e-2 3e2]);

            fig_out = AYfig_;
        end
        function fig_out = FB_appmu_vs_omega_Bingham_spread(obj, AYfig_, FB1, FB2, EC000, EC050, EC075, EC100)
            AYfig_.init_tiles([4,6]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');
            len1 = length(FB1.exp)+1;

            for i = 1:length(FB1.exp)
                % plot(axs(i), FB1.exp(i).omega, FB1.exp(i).appmu, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);

                appmu_plot = glass_particles.compute_appmu_Bingham(FB1.exp(i).tau_y_Bingham, FB1.exp(i).omega, FB1.exp(i).tau);
                appmu_plot_old = glass_particles.compute_appmu_Bingham(FB1.exp(i).tau_y, FB1.exp(i).omega, FB1.exp(i).tau);
                omega_range = [min(FB1.exp(i).omega) max(FB1.exp(i).omega)];

                plot(axs(i), FB1.exp(i).omega, appmu_plot, 's', 'Color', FB1.exp(i).color,'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', FB1.exp(i).label);
                plot(axs(i), FB1.exp(i).omega, appmu_plot_old, FB1.exp(i).specs,'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
                plot(axs(i), omega_range, FB1.exp(i).mu_p_Bingham*ones(size(omega_range)), '--', 'Color', [0 0 0],'LineWidth', 1, 'MarkerSize', 1);

            end

            for i = 1:length(FB2.exp)
                % plot(axs(i+len1), FB2.exp(i).omega, FB2.exp(i).appmu, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);

                appmu_plot = glass_particles.compute_appmu_Bingham(FB2.exp(i).tau_y_Bingham, FB2.exp(i).omega, FB2.exp(i).tau);
                appmu_plot_old = glass_particles.compute_appmu_Bingham(FB2.exp(i).tau_y, FB2.exp(i).omega, FB2.exp(i).tau);
                omega_range = [min(FB2.exp(i).omega) max(FB2.exp(i).omega)];

                plot(axs(i+len1), FB2.exp(i).omega, appmu_plot_old, FB2.exp(i).specs, 'Color', FB2.exp(i).color, 'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
                plot(axs(i+len1), FB2.exp(i).omega, appmu_plot, 's', 'Color', FB2.exp(i).color,'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', FB2.exp(i).label);
                plot(axs(i+len1), omega_range, FB2.exp(i).mu_p_Bingham*ones(size(omega_range)), '--', 'Color', [0 0 0],'LineWidth', 1, 'MarkerSize', 1);
            end

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:6:19), '$$\mu_{app}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(19:end), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            % axis(axs,obj.omega_appmu_range)
            axis(axs, [1e-2 2e2 1e-4 1e4])

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Bingham_fluid_params(obj,AYfig_,FB1,FB2)
            [tdim1,tdim2] = deal(2,2);
            axs=prep_tiles(AYfig_,[tdim1 tdim2]);
            axi=0;

            for FB = [FB1 FB2]
                explen=length(FB.exp);
                [tauy_vec mup_vec q_vec] = deal(nan(explen,1));
                [tauy_vec_old mup_vec_old] = deal(nan(explen,1));
                [mup_bounds_mat tauy_bounds_mat] = deal(nan(explen,3));
                color_mat=nan(explen,3);
                for i=1:(explen)
                    expi=FB.exp(i);

                    [tauy_vec(i) mup_vec(i)] = deal(expi.tau_y_Bingham, expi.mu_p_Bingham);
                    [tauy_vec_old(i) mup_vec_old(i)] = deal(expi.tau_y, expi.mu_p);
                    [mup_bounds_mat(i,:) tauy_bounds_mat(i,:)] = expi.determine_Bingham_fluid_bounds;
                    [color_mat(i,:) q_vec(i)] = deal(FB.exp(i).color, FB.exp(i).q);
                end
                scatter(axs(1), q_vec, tauy_vec/expi.tau_static, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                scatter(axs(2), q_vec, tauy_vec/expi.tau_static, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                % scatter(axs(1+axi), q_vec, tauy_vec/expi.tau_static, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                % scatter(axs(1+axi), q_vec, tauy_vec_old, 'p','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                % scatter(axs(1+axi), q_vec, tauy_bounds_mat(:,1)/expi.tau_static, 'v','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                % scatter(axs(1+axi), q_vec, tauy_bounds_mat(:,3)/expi.tau_static, '^','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);

                scatter(axs(3), q_vec, mup_vec, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                scatter(axs(4), q_vec, mup_vec, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                % scatter(axs(3+axi), q_vec, mup_vec, FB.specs,'CData',color_mat,'LineWidth',FB.LW_L,'SizeData',3*FB.MS_L*FB.MS_L);
                % scatter(axs(3+axi), q_vec, mup_vec_old, 'p','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                % scatter(axs(3+axi), q_vec, mup_bounds_mat(:,1), 'v','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);
                % scatter(axs(3+axi), q_vec, mup_bounds_mat(:,3), '^','CData',[0 0 0],'LineWidth',2*FB.LW_L,'SizeData',15);

                axi=axi+1;
            end
            ylabel(axs(1:2), '$$\tau_y / \tau_{q = 0}$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(3:4), '$$\mu_p$$', 'Interpreter', 'LaTeX','FontSize',14)

            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            % set(axs,'YScale','log')
            % set(axs,'YScale','linear')
            set(axs(1:2),'YScale','log')

            ylim(axs(1:2), [1e-3 5e2]);
            ylim(axs(3:4), [1e-4 1e2]);

            xlim(axs(1:2:3), [0 2]);
            xlim(axs(2:2:4), [0 16]);
            ylim(axs(1:2), [1e-5 1]);
            ylim(axs(3:4), [0 0.6]);

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
        function fig_out = FB1_FB2_Bingham_fluid_omega_fits(obj,AYfig_,FB1,FB2,full_flag_)
            if (nargin==4)
                full_flag=true;
            else
                full_flag=full_flag_;
            end

            if (full_flag)
                omega_max = 100;
                xscale_set = 'log';
                yscale_set = 'log';
            else
                omega_max = 15;
                xscale_set = 'linear';
                yscale_set = 'linear';
            end

            axs_full=prep_tiles(AYfig_,[4 6]);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    par = [FB.exp(i).mu_p_Bingham FB.exp(i).tau_y_Bingham];

                    inds_plot = FB.exp(i).omega <= omega_max;
                    omega_plot = FB.exp(i).omega(inds_plot);
                    tau_plot = FB.exp(i).tau(inds_plot);

                    ti = linspace(min(tau_plot), max(tau_plot), 100);
                    oi = glass_particles.omegai_pred_Bingham(par(1),par(2),ti);

                    plot(axs(i), tau_plot, omega_plot, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), ti, oi, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                    % plot(axs(i), omega_linspace, ti_old, ':', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            xlabel(axs_full, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            ylabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full, 'YScale', yscale_set, 'XScale', xscale_set);

            fig_out = AYfig_;
        end
    end
end

function ax_out = set_UXB_Bingham_summary_axes(AYfig_)
    ax_out = prep_tiles(AYfig_,[4,4]);
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
