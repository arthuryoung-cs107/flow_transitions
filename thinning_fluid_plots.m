classdef thinning_fluid_plots < bingham_plots
    properties

    end
    methods
        function obj = thinning_fluid_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@bingham_plots(write_figs_, write_all_figs_, figs_to_write_);
        end
        function fig_out = FB_Carreau_mueffi_vs_omegai(obj,AYfig_,FBall,muinf_flag)
            axs_full=prep_tiles(AYfig_,[1 2]);

            axi=1;
            for FB = FBall
                explen = length(FB.exp);
                for i=1:explen
                    expi = FB.exp(i);
                    if (muinf_flag)
                        pari = expi.Carreau_fit_muinf_params;
                        [mu0i, lambdai, ni, ki, muinfi] = deal(pari(1),pari(2),pari(3),pari(4),pari(5));
                        [ti,gi,mui] = expi.comp_Carreau_fluid(expi.omega, mu0i, lambdai, ni, ki, muinfi);
                    else
                        mui = expi.mu_effi_Carreau;
                    end

                    plot(axs_full(axi), expi.omega, mui, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                end
                axi = axi+1;
            end
            ylabel(axs_full, '$$\mu_{C}(\omega_i)$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)
            set(axs_full, 'YScale', 'log', 'XScale', 'log');

            fig_out = AYfig_;
        end
        function fig_out = FB_Carreau_fluid_params(obj,AYfig_,FBall,muinf_flag)
            if (muinf_flag)
                [tdim1,tdim2] = deal(4,2);
            else
                [tdim1,tdim2] = deal(3,2);
            end
            axs=prep_tiles(AYfig_,[tdim1 tdim2]);
            axi=0;
            for FB = FBall
                explen=length(FB.exp);
                Cparam_tens = nan(explen,4,4);
                q_vec = nan(explen,1);
                color_mat=nan(explen,3);
                for i=1:(explen)
                    expi=FB.exp(i);
                    if (muinf_flag)
                        Cparam_i = expi.Carreau_fit_muinf_params;
                        [mu0_b lambda_b n_b k_b muinf_b] = expi.determine_Carreau_fluid_bounds_muinf(expi.omega_fit_Bingham, expi.tau_fit_Bingham);

                        Cparam_tens(i,:,1) = [Cparam_i(1) mu0_b];
                        Cparam_tens(i,:,2) = [Cparam_i(2) lambda_b];
                        Cparam_tens(i,:,3) = [Cparam_i(3) n_b];
                        Cparam_tens(i,:,4) = [Cparam_i(5) muinf_b];
                    else
                        Cparam_i = expi.Carreau_fit_params;
                        [mu0_b lambda_b n_b k_b muinf_b] = expi.determine_Carreau_fluid_bounds(expi.omega_fit_Carreau_muinf, expi.tau_fit_Carreau_muinf);

                        Cparam_tens(i,:,1) = [Cparam_i(1) mu0_b];
                        Cparam_tens(i,:,2) = [Cparam_i(2) lambda_b];
                        Cparam_tens(i,:,3) = [Cparam_i(3) n_b];
                    end

                    [color_mat(i,:) q_vec(i)] = deal(FB.exp(i).color, FB.exp(i).q);
                end

                for j = 1:3
                    scatter(axs(axi+((j-1)*tdim2)+1), q_vec, Cparam_tens(:,1,j), 'p','CData',color_mat,'LineWidth',2,'SizeData',5);
                    scatter(axs(axi+((j-1)*tdim2)+1), q_vec, Cparam_tens(:,2,j), 'v','CData',color_mat,'LineWidth',2,'SizeData',5);
                    scatter(axs(axi+((j-1)*tdim2)+1), q_vec, Cparam_tens(:,3,j), 'o','CData',color_mat,'LineWidth',2,'SizeData',5);
                    scatter(axs(axi+((j-1)*tdim2)+1), q_vec, Cparam_tens(:,4,j), '^','CData',color_mat,'LineWidth',2,'SizeData',5);
                end
                if (muinf_flag)
                    scatter(axs(axi+(3*tdim2)+1), q_vec, Cparam_tens(:,1,4), 'p','CData',color_mat,'LineWidth',2,'SizeData',5);
                    scatter(axs(axi+(3*tdim2)+1), q_vec, Cparam_tens(:,2,4), 'v','CData',color_mat,'LineWidth',2,'SizeData',5);
                    scatter(axs(axi+(3*tdim2)+1), q_vec, Cparam_tens(:,3,4), 'o','CData',color_mat,'LineWidth',2,'SizeData',5);
                    scatter(axs(axi+(3*tdim2)+1), q_vec, Cparam_tens(:,4,4), '^','CData',color_mat,'LineWidth',2,'SizeData',5);
                end
                axi=axi+1;
            end
            % set(axs(1:4),'YScale','log')
            set(axs,'XScale', 'log','YScale','log')
            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(1:2), '$$\mu_0$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(3:4), '$$\lambda$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel(axs(5:6), '$$n$$', 'Interpreter', 'LaTeX','FontSize',14)
            if (muinf_flag)
                set(axs(7:8),'YScale','log')
                ylabel(axs(7:8), '$$\mu_{\infty}$$', 'Interpreter', 'LaTeX','FontSize',14)
            end

            fig_out=AYfig_;
        end
        function fig_out = FB_Carreau_fits_spread(obj,AYfig_,FBall,muinf_flag,log_flag_)
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

            axs_full=prep_tiles(AYfig_,[4 6]);

            axs=axs_full;
            for FB = FBall
                explen = length(FB.exp);
                for i=1:explen
                    if (muinf_flag)
                        % par_i = FB.exp(i).fit_Carreau_model_muinf;
                        % FB.exp(i).Carreau_fit_muinf_params = par_i;
                        par_i = FB.exp(i).Carreau_fit_muinf_params;
                        omega_fit = FB.exp(i).omega_fit_Carreau_muinf;
                        tau_fit = FB.exp(i).tau_fit_Carreau_muinf;
                    else
                        par_i = FB.exp(i).Carreau_fit_params;
                        omega_fit = FB.exp(i).omega_fit_Carreau;
                        tau_fit = FB.exp(i).tau_fit_Carreau;
                    end

                    if (log_flag)
                        ti = glass_particles.taui_pred_Carreau(par_i,omega_test);

                        plot(axs(i), omega_test, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                        plot(axs(i), omega_fit, tau_fit, 'o', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', FB.exp(i).label);

                        [omega_not_fit i_nf] = setdiff(FB.exp(i).omega, omega_fit);
                        tau_not_fit = FB.exp(i).tau(i_nf);

                        plot(axs(i), omega_not_fit, tau_not_fit, 'x', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', FB.exp(i).label);

                        omega_crit = 1/(par_i(4)*par_i(2));
                        if ((omega_crit>omega_test(1))&&(omega_crit<omega_test(end)))
                            tau_crit = glass_particles.taui_pred_Carreau(par_i,omega_crit);
                            plot(axs(i), omega_crit*ones(1,2), [min([min(ti) min(tau_fit)]) tau_crit], '- p', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 8);
                        end
                    else
                        omega_test = linspace(0, max(omega_fit), 100);
                        ti = glass_particles.taui_pred_Carreau(par_i,omega_test);

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
        function fig_out = FB_Carreau_err_vs_omegai(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    expi = FB.exp(i);
                    pari = expi.Carreau_fit_params;
                    [mu0i, lambdai, ni, ki, muinfi] = deal(pari(1),pari(2),pari(3),pari(4),pari(5));
                    [ti,gi,mui] = expi.comp_Carreau_fluid(expi.omega_fit_Carreau, mu0i, lambdai, ni, ki, muinfi);

                    tau_err = abs(expi.tau_fit_Carreau - ti)./expi.tau_fit_Carreau;

                    plot(axs(i), expi.omega_fit_Carreau, tau_err, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            ylabel(axs_full, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            fig_out = AYfig_;
        end
        function fig_out = ALL_dtdo_vs_omegai_spread(obj,AYfig_,FBall,PFall,NBall,UBall,XBall,log_flag_)
            if (nargin==7)
                log_flag = true;
            else
                log_flag = log_flag_;
            end
            if (log_flag)
                axscale = 'log';
                xlim_set = [1e-2 1e2];
                ylim_set = [-1 1];
            else
                axscale = 'linear';
                xlim_set = [0 130];
                ylim_set = [-1 1];
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
                    plot(axs(i), expi.dtdo_o, expi.dtdo_d, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    title(axs(i), expi.label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+1:end);
            end

            plot(axs_full(24), PF1.dtdo_o, PF1.dtdo_d, PF1.specs, 'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            plot(axs_full(24), PFR.dtdo_o, PFR.dtdo_d, PFR.specs, 'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            plot(axs_full(25), NB1.dtdo_o, NB1.dtdo_d, NB1.specs, 'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            plot(axs_full(25), NB2.dtdo_o, NB2.dtdo_d, NB2.specs, 'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            plot(axs_full(25), NB3.dtdo_o, NB3.dtdo_d, NB3.specs, 'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            plot(axs_full(26), UB1.dtdo_o, UB1.dtdo_d, UB1.specs, 'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            plot(axs_full(26), UB2.dtdo_o, UB2.dtdo_d, UB2.specs, 'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            plot(axs_full(27), XB1.dtdo_o, XB1.dtdo_d, XB1.specs, 'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            plot(axs_full(27), XB2.dtdo_o, XB2.dtdo_d, XB2.specs, 'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            ylabel(axs_full(1:dims(2):22), '$$\frac{d \tau}{d \omega}$$ ', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full(22:end), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full,'XScale',axscale,'YScale','linear');
            xlim(axs_full,xlim_set);
            ylim(axs_full,ylim_set);

            fig_out = AYfig_;
        end
        function fig_out = ALL_d2tdo2_vs_omegai_spread(obj,AYfig_,FBall,PFall,NBall,UBall,XBall,log_flag_)
            if (nargin==7)
                log_flag = true;
            else
                log_flag = log_flag_;
            end
            if (log_flag)
                axscale = 'log';
                xlim_set = [1e-2 1e2];
                ylim_set = [-1 1];
            else
                axscale = 'linear';
                xlim_set = [0 130];
                ylim_set = [-1 1];
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
                    plot(axs(i), expi.d2tdo2_o, expi.d2tdo2_d, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    title(axs(i), expi.label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+1:end);
            end

            plot(axs_full(24), PF1.d2tdo2_o, PF1.d2tdo2_d, PF1.specs, 'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            plot(axs_full(24), PFR.d2tdo2_o, PFR.d2tdo2_d, PFR.specs, 'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);

            plot(axs_full(25), NB1.d2tdo2_o, NB1.d2tdo2_d, NB1.specs, 'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            plot(axs_full(25), NB2.d2tdo2_o, NB2.d2tdo2_d, NB2.specs, 'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            plot(axs_full(25), NB3.d2tdo2_o, NB3.d2tdo2_d, NB3.specs, 'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);

            plot(axs_full(26), UB1.d2tdo2_o, UB1.d2tdo2_d, UB1.specs, 'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            plot(axs_full(26), UB2.d2tdo2_o, UB2.d2tdo2_d, UB2.specs, 'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);

            plot(axs_full(27), XB1.d2tdo2_o, XB1.d2tdo2_d, XB1.specs, 'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            plot(axs_full(27), XB2.d2tdo2_o, XB2.d2tdo2_d, XB2.specs, 'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            ylabel(axs_full(1:dims(2):22), '$$\frac{d^2 \tau}{d \omega^2}$$ ', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full(22:end), '$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full,'XScale',axscale,'YScale','linear');
            xlim(axs_full,xlim_set);
            ylim(axs_full,ylim_set);

            fig_out = AYfig_;
        end
        function fig_out = FB_power_fluid_fits(obj,AYfig_,FB,tile_dims_)
            axs=prep_tiles(AYfig_,tile_dims_);

            omega_cap = 6;
            omega_linspace=linspace(0,omega_cap,100);

            for i=1:length(FB.exp)
                indi=FB.exp(i).omega<omega_cap;
                [Ki,ni] = FB.exp(i).fit_power_fluid(indi);
                FB.exp(i).power_fit_params=[Ki,ni];
                [ti,gi] = FB.exp(i).comp_power_fluid(omega_linspace, Ki, ni);

                plot(axs(i), FB.exp(i).omega(indi), FB.exp(i).tau(indi), FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 2*FB.MS, 'DisplayName', FB.exp(i).label);

                plot(axs(i), omega_linspace, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
            end

            ylabel(axs, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)
            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_fluid_params(obj,AYfig_,FB1,FB2)
            [tdim1,tdim2] = deal(2,3);
            axs=prep_tiles(AYfig_,[tdim1 tdim2]);
            axi=0;

            for FB = [FB1 FB2]
                explen=length(FB.exp);
                [mu0_vec lambda_vec n_vec q_vec] = deal(nan(explen,1));
                color_mat=nan(explen,3);

                for i=1:(explen)
                    % [mu0_b lambda_b n_b k_b muinf_b] = FB.exp(i).determine_Carreau_fluid_bounds(FB.exp(i).omega(indi), FB.exp(i).tau(indi));
                    expi=FB.exp(i);

                    [mu0_vec(i) lambda_vec(i) n_vec(i)] = deal(expi.mu0_Carreau, expi.lambda_Carreau, expi.n_Carreau);
                    [color_mat(i,:) q_vec(i)] = deal(FB.exp(i).color, FB.exp(i).q);

                    % fprintf('(%s) i:%d, mu0: %e (%.1f%%), l: %e (%.1f%%), n: %e (%.1f%%), k: %e (%.1f%%), muinf: %e (%.1f%%)\n', ...
                    % FB.exp(i).label,i, ...
                    % mu0i, 100*(mu0i-mu0_b(1))/(mu0_b(3)-mu0_b(1)), ...
                    % lambdai, 100*(lambdai-lambda_b(1))/(lambda_b(3)-lambda_b(1)), ...
                    % ni, 100*(ni-n_b(1))/(n_b(3)-n_b(1)), ...
                    % ki, 100*(ki-k_b(1))/(k_b(3)-k_b(1)), ...
                    % muinfi, 100*(muinfi-muinf_b(1))/(muinf_b(3)-muinf_b(1)) );
                end
                scatter(axs(1+axi), q_vec, mu0_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                scatter(axs(2+axi), q_vec, lambda_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                scatter(axs(3+axi), q_vec, n_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                % scatter(axs(4+axi), q_vec, k_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                % scatter(axs(5+axi), q_vec, muinf_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                axi=axi+tdim2;
            end
            ylabel([axs(1) axs(1+tdim2)], '$$\mu_0$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(2) axs(2+tdim2)], '$$\lambda$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(3) axs(3+tdim2)], '$$n$$', 'Interpreter', 'LaTeX','FontSize',14)
            % ylabel([axs(4) axs(9)], '$$k$$', 'Interpreter', 'LaTeX','FontSize',14)
            % ylabel([axs(5) axs(10)], '$$\mu_{\infty}$$', 'Interpreter', 'LaTeX','FontSize',14)

            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            set(axs,'YScale','log')
            set([axs(3) axs(3+tdim2)],'YScale','linear')

            fig_out=AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_fit_results(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[2 5]);

            omega_cap = 15;
            omega_linspace=linspace(0,omega_cap,100);

            axi=0;

            for FB = [FB1 FB2]
                explen=length(FB.exp);
                [mu0_vec lambda_vec n_vec k_vec muinf_vec q_vec] = deal(nan(explen,1));
                color_mat=nan(explen,3);

                for i=1:(explen)
                    indi=FB.exp(i).omega<omega_cap;
                    [mu0i lambdai ni ki muinfi] = FB.exp(i).fit_Carreau_fluid(omega_cap,true);

                    [mu0_b lambda_b n_b k_b muinf_b] = FB.exp(i).determine_Carreau_fluid_bounds(FB.exp(i).omega(indi), FB.exp(i).tau(indi));

                    [mu0_vec(i) lambda_vec(i) n_vec(i) k_vec(i) muinf_vec(i)] = deal(mu0i,lambdai,ni,ki,muinfi);
                    [color_mat(i,:) q_vec(i)] = deal(FB.exp(i).color, FB.exp(i).q);

                    FB.exp(i).Carreau_fit_params=[mu0i lambdai ni ki muinfi];


                    fprintf('(%s) i:%d, mu0: %e (%.1f%%), l: %e (%.1f%%), n: %e (%.1f%%), k: %e (%.1f%%), muinf: %e (%.1f%%)\n', ...
                    FB.exp(i).label,i, ...
                    mu0i, 100*(mu0i-mu0_b(1))/(mu0_b(3)-mu0_b(1)), ...
                    lambdai, 100*(lambdai-lambda_b(1))/(lambda_b(3)-lambda_b(1)), ...
                    ni, 100*(ni-n_b(1))/(n_b(3)-n_b(1)), ...
                    ki, 100*(ki-k_b(1))/(k_b(3)-k_b(1)), ...
                    muinfi, 100*(muinfi-muinf_b(1))/(muinf_b(3)-muinf_b(1)) );
                end
                scatter(axs(1+axi), q_vec, mu0_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                scatter(axs(2+axi), q_vec, lambda_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                scatter(axs(3+axi), q_vec, n_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                scatter(axs(4+axi), q_vec, k_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                scatter(axs(5+axi), q_vec, muinf_vec, FB.specs,'CData',color_mat,'LineWidth',2*FB.LW_L,'SizeData',5*FB.MS_L*FB.MS_L);
                axi=axi+5;
            end
            ylabel([axs(1) axs(6)], '$$\mu_0$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(2) axs(7)], '$$\lambda$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(3) axs(8)], '$$n$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(4) axs(9)], '$$k$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(5) axs(10)], '$$\mu_{\infty}$$', 'Interpreter', 'LaTeX','FontSize',14)

            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            set(axs,'YScale','log')
            set([axs(3) axs(8)],'YScale','linear')

            fig_out=AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_viscosity(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);

            omega_cap = 126;
            omega_linspace=linspace(0,omega_cap,100);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen=length(FB.exp);
                for i=1:(explen)
                    exp=FB.exp(i);
                    [taui, gammai, mueffi] = exp.comp_Carreau_fluid(omega_linspace);
                    plot(axs(i), omega_linspace, mueffi, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    plot(axs(i), omega_linspace, exp.mu_p*ones(size(omega_linspace)), '-', 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 2*FB.MS, 'DisplayName', FB.exp(i).label);
                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            set(axs_full, 'YScale', 'log')
            ylim(axs_full, [1e-2 1e4]);
            fig_out=AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_fluid_fits_full(obj,AYfig_,FB1,FB2,full_flag_)
            if (nargin==4)
                full_flag=true;
            else
                full_flag=full_flag_;
            end

            if (full_flag)
                omega_linspace=linspace(1e-2, 1e2,500);
                xscale_set = 'log';
                yscale_set = 'log';
            else
                omega_linspace=linspace(0, 15,500);
                xscale_set = 'linear';
                yscale_set = 'linear';
            end

            axs_full=prep_tiles(AYfig_,[4 6]);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    pari = FB.exp(i).Carreau_fit_params;
                    [mu0i, lambdai, ni, ki, muinfi] = deal(pari(1),pari(2),pari(3),pari(4),pari(5));
                    [ti,gi,mui] = FB.exp(i).comp_Carreau_fluid(omega_linspace, mu0i, lambdai, ni, ki, muinfi);

                    inds_plot = FB.exp(i).omega <= max(omega_linspace);
                    omega_plot = FB.exp(i).omega(inds_plot);
                    tau_plot = FB.exp(i).tau(inds_plot);

                    % plot(axs(i), omega_linspace, mui, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                    % plot(axs(i), omega_linspace, FB.exp(i).mu_p_Bingham*ones(size(omega_linspace)), ':', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), omega_plot, tau_plot, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), omega_linspace, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                    % plot(axs(i), 1/(lambdai*ki), FB.exp(i).comp_Carreau_fluid(1/(lambdai*ki), mu0i, lambdai, ni, ki, muinfi), 'p', 'Color', [0 0 0], 'LineWidth', 2);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            ylabel(axs_full, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            % ylabel(axs_full, '$$\mu$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            set(axs_full, 'YScale', yscale_set, 'XScale', xscale_set);

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_fluid_fits(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);

            omega_cap = 15;
            omega_linspace=linspace(0,omega_cap,100);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    indi=FB.exp(i).omega<omega_cap;
                    pari = FB.exp(i).Carreau_fit_params;
                    [mu0i, lambdai, ni, ki, muinfi] = deal(pari(1),pari(2),pari(3),pari(4),pari(5));
                    [ti,gi,mui] = FB.exp(i).comp_Carreau_fluid(omega_linspace, mu0i, lambdai, ni, ki, muinfi);

                    plot(axs(i), FB.exp(i).omega(indi), FB.exp(i).tau(indi), FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), omega_linspace, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                    plot(axs(i), 1/(lambdai*ki), FB.exp(i).comp_Carreau_fluid(1/(lambdai*ki), mu0i, lambdai, ni, ki, muinfi), 'p', 'Color', [0 0 0], 'LineWidth', 2);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            ylabel(axs_full, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs_full, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            fig_out = AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_fluid_G_vs_Res(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);

            s2g=0.5*(obj.r_o+obj.r_i)/obj.r_o;
            g2s=1/s2g;


            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    exp = FB.exp(i);
                    [exp.G_Carreau exp.Re_s_Carreau] = FB.exp(i).comp_G_Res_Carreau_fluid(FB.exp(i).Carreau_fit_params);

                    plot(axs(i), g2s*exp.Re_b_Carreau, exp.G_b_Carreau, 's', 'Color', FB.exp(i).color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', FB.exp(i).label);
                    fplot(axs(i), @(Re) obj.G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 0.5, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            set(axs_full, 'XScale', 'log', 'YScale', 'log')
            ylabel(axs_full, '$$ G $$', 'Interpreter', 'Latex')
            xlabel(axs_full, '$$ Re_s = \frac{r_o-r_i}{2 r_o} Re_b$$', 'Interpreter', 'Latex')
            % axis(axs_full,obj.Res_G_range)
            fig_out=AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_fluid_Grat_vs_Res(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                q_vec = nan(explen,1);
                Re_sc_vec = nan(explen,1);
                for i=1:explen
                    [taui, gammai, mueffi] = FB.exp(i).comp_Carreau_fluid(FB.exp(i).omega);
                    tau_rat=FB.exp(i).tau./taui;

                    plot(axs(i), FB.exp(i).Re_s_Carreau, tau_rat, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    inds_full=1:length(tau_rat);
                    inds_trans = inds_full(logical((FB.exp(i).omega>15).*(tau_rat>1.5)));
                    if (~isempty(inds_trans))
                        Re_sc = FB.exp(i).Re_s_Carreau(inds_trans(1));
                        tau_rat_trans = tau_rat(inds_trans(1));
                    else
                        Re_sc = NaN;
                        tau_rat_trans = NaN;
                    end

                    q_vec(i) = FB.exp(i).q;
                    Re_sc_vec(i) = Re_sc;

                    plot(axs(i), Re_sc, tau_rat_trans, 'p', 'Color', [0 0 0], 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    fprintf('%s : Re_{s,c} = %e\n', FB.exp(i).label, Re_sc)

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            % set(axs_full, 'XScale', 'log', 'YScale', 'log')
            set(axs_full, 'XScale', 'log')
            ylim(axs_full, [1e-1 50]);

            % plot(axs_full(11), q_vec, Re_sc_vec, 'p', 'Color', [0 0 0], 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
            % set(axs_full(11), 'XScale', 'linear', 'YScale', 'log')

            fig_out=AYfig_;
        end
        function fig_out = FB1_FB2_Carreau_fluid_Res_vs_omegai(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);

            axs=axs_full;
            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    plot(axs(i), FB.exp(i).omega, FB.exp(i).Re_s_Carreau, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            set(axs_full, 'XScale', 'log', 'YScale', 'log')
            % set(axs_full, 'XScale', 'log')
            % ylim(axs_full, [1e-2 1e4]);
            fig_out=AYfig_;
        end
        function fig_out = ALL_G_vs_Res_CarreauFB(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall)
            axs=prep_tiles(AYfig_, [3,2]);
            s2g=0.5*(obj.r_o+obj.r_i)/obj.r_o;
            g2s=1/s2g;

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
                if (FB1.exp(i).q <1)
                    legend_set_e(i) = plot(axs(5), FB1.exp(i).Re_s, FB1.exp(i).G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
                else
                    legend_set_e(i) = plot(axs(5), g2s*FB1.exp(i).Re_b_Carreau, FB1.exp(i).G_b_Carreau, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
                end
            end
            % fplot(axs(5), @(Re) (FB1.powerfit.b).*(Re).^(FB1.powerfit.m), [71 10000],'-', 'Color', FB1.color,'Linewidth', 2,'DisplayName', 'FB1 $$\beta Re_s^{\alpha}$$');

            for i=1:length(FB2.exp)
                if (FB2.exp(i).q < 1)
                    legend_set_f(i) = plot(axs(6), FB2.exp(i).Re_s, FB2.exp(i).G, FB2.specs, 'Color',     FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
                else
                    legend_set_f(i) = plot(axs(6), g2s*FB2.exp(i).Re_b_Carreau, FB2.exp(i).G_b_Carreau, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
                end
            end
            % fplot(axs(6), @(Re) (FB2.powerfit.b).*(Re).^(FB2.powerfit.m), [71 10000],'-', 'Color', FB2.color,'Linewidth', 2, 'DisplayName', 'FB2 $$\beta Re_s^{\alpha}$$');

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
        function fig_out = FB1_FB2_Carreau_tau_vs_gamma(obj,AYfig_,FB1,FB2)
            axs_full=prep_tiles(AYfig_,[4 6]);
            axs=axs_full;
            % gamma_logspace=logspace(-6,6,1000);
            % gamma_logspace=linspace(0,250,10000);
            gamma_logspace=logspace(-3,4,1000);

            for FB = [FB1 FB2]
                explen = length(FB.exp);
                for i=1:explen
                    exp =FB.exp(i);
                    pari=FB.exp(i).Carreau_fit_params;
                    [mu0, l, n, k, muinf] = deal(pari(1),pari(2),pari(3),pari(4),pari(5));

                    mueff = @(s) muinf+(mu0-muinf)*((1+(l*s).*(l*s)).^(0.5*(n-1)));
                    tau = @(s) mueff(s).*(s);

                    % plot(axs(i), exp.omega, exp.S_Carreau_typ, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    % plot(axs(i), tau(gamma_logspace), mueff(gamma_logspace), '-', 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    yyaxis(axs(i), 'left')
                    plot(axs(i),exp.omega, exp.S_Carreau_typ./abs(exp.gammai_Carreau), exp.specs, 'Color', [0 0 1], 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    set(axs(i), 'XScale', 'log', 'YScale', 'log')
                    ylim(axs(i), [min(exp.S_Carreau_typ./abs(exp.gammai_Carreau)) max(exp.S_Carreau_typ./abs(exp.gammai_Carreau))]);

                    yyaxis(axs(i), 'right')
                    plot(axs(i),exp.omega, exp.mu_Carreau_typ./exp.mu_effi_Carreau, 'o', 'Color', [1 0 0], 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    set(axs(i), 'XScale', 'log', 'YScale', 'log')
                    ylim(axs(i), [1/max(exp.S_Carreau_typ./abs(exp.gammai_Carreau)) 1/min(exp.S_Carreau_typ./abs(exp.gammai_Carreau))]);

                    % plot(axs(i),exp.omega, (exp.S_Carreau_typ./abs(exp.gammai_Carreau)).*(exp.mu_Carreau_typ./exp.mu_effi_Carreau), 'o', 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    % plot(axs(i), abs(exp.gammai_Carreau), exp.S_Carreau_typ, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);
                    % plot(axs(i), gamma_linspace, mueff(gamma_linspace), FB.specs, 'Color', FB.exp(i).color, 'LineWidth', FB.LW, 'MarkerSize', FB.MS, 'DisplayName', FB.exp(i).label);

                    title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
                end
                axs=axs_full(explen+2:end);
            end
            set(axs_full, 'XScale', 'log', 'YScale', 'log')
            % ylim(axs_full, [1e-2 1e4 ]);
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
