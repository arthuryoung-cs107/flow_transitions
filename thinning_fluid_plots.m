classdef thinning_fluid_plots < bingham_plots
    properties

    end
    methods
        function obj = thinning_fluid_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@bingham_plots(write_figs_, write_all_figs_, figs_to_write_);
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
        function fig_out = FB1_FB2_Carreau_fit_results(obj,AYfig_,FB1,FB2)
            axs=prep_tiles(AYfig_,[2 3]);

            omega_cap = 15;
            omega_linspace=linspace(0,omega_cap,100);

            axi=0;

            for FB = [FB1 FB2]
                explen=length(FB.exp);
                [mu0_vec lambda_vec n_vec k_vec muinf_vec q_vec] = deal(nan(explen,1));
                color_mat=nan(explen,3);

                for i=1:(explen)
                    indi=FB.exp(i).omega<omega_cap;
                    [mu0i lambdai ni ki muinfi] = FB.exp(i).fit_Carreau_fluid(indi,true);

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
                axi=axi+3;
            end
            ylabel([axs(1) axs(4)], '$$\mu_0$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(2) axs(5)], '$$\lambda$$', 'Interpreter', 'LaTeX','FontSize',14)
            ylabel([axs(3) axs(6)], '$$n$$', 'Interpreter', 'LaTeX','FontSize',14)

            xlabel(axs, '$$q= Q/Q_{inc}$$', 'Interpreter', 'LaTeX','FontSize',14)
            set(axs,'YScale','log')

            fig_out=AYfig_;
        end

        function fig_out = FB_Carreau_fluid_fits(obj,AYfig_,FB,tile_dims_)
            axs_full=prep_tiles(AYfig_,tile_dims_);
            explen = length(FB.exp);
            axs=axs_full(1:explen);

            omega_cap = 15;
            omega_linspace=linspace(0,omega_cap,100);

            for i=1:explen
                indi=FB.exp(i).omega<omega_cap;
                pari = FB.exp(i).Carreau_fit_params;
                [mu0i, lambdai, ni, ki, muinfi] = deal(pari(1),pari(2),pari(3),pari(4),pari(5));
                [ti,gi,mui] = FB.exp(i).comp_Carreau_fluid(omega_linspace, mu0i, lambdai, ni, ki, muinfi);

                plot(axs(i), FB.exp(i).omega(indi), FB.exp(i).tau(indi), FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 2*FB.MS, 'DisplayName', FB.exp(i).label);
                plot(axs(i), omega_linspace, ti, '-', 'Color', FB.exp(i).color, 'LineWidth', 2, 'DisplayName', FB.exp(i).label);
                plot(axs(i), 1/(lambdai*ki), FB.exp(i).comp_Carreau_fluid(1/(lambdai*ki), mu0i, lambdai, ni, ki, muinfi), 'p', 'Color', [0 0 0], 'LineWidth', 2);

                title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
            end

            ylabel(axs, '$$\tau_w$$ [Pa]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs, '$$\omega_i$$ [rad.s]', 'Interpreter', 'LaTeX','FontSize',12)

            fig_out = AYfig_;
        end
        function fig_out = FB_Carreau_fluid_G_vs_Res(obj,AYfig_,FB,tile_dims_)
            axs_full=prep_tiles(AYfig_,tile_dims_);
            axlen=length(axs_full);
            explen = length(FB.exp);
            axs=axs_full(1:explen);

            for i=1:explen
                [Gi, Re_si] = FB.exp(i).comp_G_Res_Carreau_fluid(FB.exp(i).Carreau_fit_params);

                plot(axs(i), FB.exp(i).Re_s, FB.exp(i).G, FB.specs, 'Color', FB.exp(i).color, 'LineWidth', 3*FB.LW, 'MarkerSize', 2*FB.MS, 'DisplayName', FB.exp(i).label);

                plot(axs(i), Re_si, Gi, 's', 'Color', FB.exp(i).color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', FB.exp(i).label);

                title(axs(i), FB.exp(i).label, 'Interpreter', 'Latex', 'Fontsize', 14)
            end
        end
    end
end

function ax_out = prep_tiles(AYfig_,dims_)
    AYfig_.init_tiles(dims_);
    ax_out=AYfig_.ax_tile;
    hold(ax_out,'on');
    box(ax_out, 'on');
end
