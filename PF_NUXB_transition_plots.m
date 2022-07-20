classdef PF_NUXB_transition_plots < main_plots
    properties

    end
    methods
        function obj = PF_NUXB_transition_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@main_plots(write_figs_, write_all_figs_, figs_to_write_);
        end
        function fig_out = PF_NUXB_alpha_vs_Res_trans(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, NBall, UBall, XBall)
            AYfig_.init_tiles([3,3]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            pLW = 2*0.75;
            pMS = 2*5;

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.Re_s, PF1.alpha, [PF1.specs '-'],'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_b(1) = plot(axs(2), PFR.Re_s, PFR.alpha, [PFR.specs '-'],'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            legend_set_c(1) = plot(axs(3), NB1.Re_s, NB1.alpha, [NB1.specs '-'],'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_d(1) = plot(axs(4), NB2.Re_s, NB2.alpha, [NB2.specs '-'],'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_e(1) = plot(axs(5), NB3.Re_s, NB3.alpha, [NB3.specs '-'],'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            legend_set_f(1) = plot(axs(6), UB1.Re_s, UB1.alpha, [UB1.specs '-'],'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_g(1) = plot(axs(7), UB2.Re_s, UB2.alpha, [UB2.specs '-'],'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
            legend_set_h(1) = plot(axs(8), XB1.Re_s, XB1.alpha, [XB1.specs '-'],'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_i(1) = plot(axs(9), XB2.Re_s, XB2.alpha, [XB2.specs '-'],'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);

            legend_set_a(2) = plot(axs(1), PF1.Re_s_TV, PF1.alpha_TV, PF1.specs,'Color',[0 0 0],'LineWidth',PF1.LW,'MarkerSize', PF1.MS, 'DisplayName', [PF1.label ', transitioned']);
            legend_set_b(2) = plot(axs(2), PFR.Re_s_TV, PFR.alpha_TV, PFR.specs,'Color',[0 0 0],'LineWidth',PFR.LW,'MarkerSize', PFR.MS, 'DisplayName', [PFR.label ', transitioned']);
            legend_set_c(2) = plot(axs(3), NB1.Re_s_TV, NB1.alpha_TV, NB1.specs,'Color',[0 0 0],'LineWidth',NB1.LW,'MarkerSize', NB1.MS, 'DisplayName', [NB1.label ', transitioned']);
            legend_set_d(2) = plot(axs(4), NB2.Re_s_TV, NB2.alpha_TV, NB2.specs,'Color',[0 0 0],'LineWidth',NB2.LW,'MarkerSize', NB2.MS, 'DisplayName', [NB2.label ', transitioned']);
            legend_set_e(2) = plot(axs(5), NB3.Re_s_TV, NB3.alpha_TV, NB3.specs,'Color',[0 0 0],'LineWidth',NB3.LW,'MarkerSize', NB3.MS, 'DisplayName', [NB3.label ', transitioned']);
            legend_set_f(2) = plot(axs(6), UB1.Re_s_TV, UB1.alpha_TV, UB1.specs,'Color',[0 0 0],'LineWidth',UB1.LW,'MarkerSize', UB1.MS, 'DisplayName', [UB1.label ', transitioned']);
            legend_set_g(2) = plot(axs(7), UB2.Re_s_TV, UB2.alpha_TV, UB2.specs,'Color',[0 0 0],'LineWidth',UB2.LW,'MarkerSize', UB2.MS, 'DisplayName', [UB2.label ', transitioned']);
            legend_set_h(2) = plot(axs(8), XB1.Re_s_TV, XB1.alpha_TV, XB1.specs,'Color',[0 0 0],'LineWidth',XB1.LW,'MarkerSize', XB1.MS, 'DisplayName', [XB1.label ', transitioned']);
            legend_set_i(2) = plot(axs(9), XB2.Re_s_TV, XB2.alpha_TV, XB2.specs,'Color',[0 0 0],'LineWidth',XB2.LW,'MarkerSize', XB2.MS, 'DisplayName', [XB2.label ', transitioned']);

            legend_set_a(3) = plot(axs(1), PF1.Re_sc3, PF1.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [PF1.label ', transitioned']);
            legend_set_b(3) = plot(axs(2), PFR.Re_sc3, PFR.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [PFR.label ', transitioned']);
            legend_set_c(3) = plot(axs(3), NB1.Re_sc3, NB1.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [NB1.label ', transitioned']);
            legend_set_d(3) = plot(axs(4), NB2.Re_sc3, NB2.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [NB2.label ', transitioned']);
            legend_set_e(3) = plot(axs(5), NB3.Re_sc3, NB3.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [NB3.label ', transitioned']);
            legend_set_f(3) = plot(axs(6), UB1.Re_sc3, UB1.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [UB1.label ', transitioned']);
            legend_set_g(3) = plot(axs(7), UB2.Re_sc3, UB2.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [UB2.label ', transitioned']);
            legend_set_h(3) = plot(axs(8), XB1.Re_sc3, XB1.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [XB1.label ', transitioned']);
            legend_set_i(3) = plot(axs(9), XB2.Re_sc3, XB2.alpha_tol,'|','Color',[0 0 0],'LineWidth',pLW,'MarkerSize', pMS, 'DisplayName', [XB2.label ', transitioned']);

            set(axs,'XScale', 'log');

            ylabel(axs(1:3:7), '$$\alpha$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(7:9), '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(3), legend_set_c,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(4), legend_set_d,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(5), legend_set_e,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(6), legend_set_f,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(7), legend_set_g,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(8), legend_set_h,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(9), legend_set_i,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);

            axis(axs, [obj.Re_s_range -0.5 obj.alpha_range(2)])

            fig_out = AYfig_;
        end
        function fig_out = PF_NUXB_G_vs_Res_trans(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, NBall, UBall, XBall)
            AYfig_.init_tiles([3,3]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i = 1:9
                fplot(axs(i), @(Re) obj.G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 1, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
            end

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.Re_s, PF1.G, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = plot(axs(1), PF1.Re_s_TV, PF1.G_TV, PF1.specs,'Color', [0 0 0], 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', [PF1.label ', transitioned']);
            legend_set_a(3) = plot(axs(1), PF1.Re_sc2, PF1.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*PF1.LW, 'MarkerSize', 2*PF1.MS, 'DisplayName', [PF1.label ', $$Re_{s,c,2}$$']);
            legend_set_a(4) = plot(axs(1), PF1.Re_sc3, PF1.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*PF1.LW, 'MarkerSize', 4*PF1.MS, 'DisplayName', [PF1.label ', $$Re_{s,c,3}$$']);
            legend_set_a(5) = fplot(axs(1), @(Re) (PF1.powerfit.b).*(Re).^(PF1.powerfit.m), [71 10000],'-', 'Color', PF1.color,'Linewidth', 1, 'DisplayName', [PF1.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_b(1) = plot(axs(2), PFR.Re_s, PFR.G, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            legend_set_b(2) = plot(axs(2), PFR.Re_s_TV, PFR.G_TV, PFR.specs,'Color', [0 0 0], 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', [PFR.label ', transitioned']);
            legend_set_b(3) = plot(axs(2), PFR.Re_sc2, PFR.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*PFR.LW, 'MarkerSize', 2*PFR.MS, 'DisplayName',    [PFR.label ', $$Re_{s,c,2}$$']);
            legend_set_b(4) = plot(axs(2), PFR.Re_sc3, PFR.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*PFR.LW, 'MarkerSize', 4*PFR.MS, 'DisplayName',    [PFR.label ', $$Re_{s,c,3}$$']);
            legend_set_b(5) = fplot(axs(2), @(Re) (PFR.powerfit.b).*(Re).^(PFR.powerfit.m), [71 10000],'-', 'Color', PFR.color,'Linewidth', 1, 'DisplayName', [PFR.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_c(1) = plot(axs(3), NB1.Re_s, NB1.G, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_c(2) = plot(axs(3), NB1.Re_s_TV, NB1.G_TV, NB1.specs,'Color', [0 0 0], 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName',   [NB1.label ', transitioned']);
            legend_set_c(3) = plot(axs(3), NB1.Re_sc2, NB1.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*NB1.LW, 'MarkerSize', 2*NB1.MS, 'DisplayName',      [NB1.label ', $$Re_{s,c,2}$$']);
            legend_set_c(4) = plot(axs(3), NB1.Re_sc3, NB1.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*NB1.LW, 'MarkerSize', 4*NB1.MS, 'DisplayName',      [NB1.label ', $$Re_{s,c,3}$$']);
            legend_set_c(5) = fplot(axs(3), @(Re) (NB1.powerfit.b).*(Re).^(NB1.powerfit.m), [71 10000],'-', 'Color', NB1.color,'Linewidth', 1, 'DisplayName', [NB1.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_d(1) = plot(axs(4), NB2.Re_s, NB2.G, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_d(2) = plot(axs(4), NB2.Re_s_TV, NB2.G_TV, NB2.specs,'Color', [0 0 0], 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName',   [NB2.label ', transitioned']);
            legend_set_d(3) = plot(axs(4), NB2.Re_sc2, NB2.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*NB2.LW, 'MarkerSize', 2*NB2.MS, 'DisplayName',      [NB2.label ', $$Re_{s,c,2}$$']);
            legend_set_d(4) = plot(axs(4), NB2.Re_sc3, NB2.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*NB2.LW, 'MarkerSize', 4*NB2.MS, 'DisplayName',      [NB2.label ', $$Re_{s,c,3}$$']);
            legend_set_d(5) = fplot(axs(4), @(Re) (NB2.powerfit.b).*(Re).^(NB2.powerfit.m), [71 10000],'-', 'Color', NB2.color,'Linewidth', 1, 'DisplayName', [NB2.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_e(1) = plot(axs(5), NB3.Re_s, NB3.G, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            legend_set_e(2) = plot(axs(5), NB3.Re_s_TV, NB3.G_TV, NB3.specs,'Color', [0 0 0], 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName',   [NB3.label ', transitioned']);
            legend_set_e(3) = plot(axs(5), NB3.Re_sc2, NB3.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*NB3.LW, 'MarkerSize', 2*NB3.MS, 'DisplayName',      [NB3.label ', $$Re_{s,c,2}$$']);
            legend_set_e(4) = plot(axs(5), NB3.Re_sc3, NB3.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*NB3.LW, 'MarkerSize', 4*NB3.MS, 'DisplayName',      [NB3.label ', $$Re_{s,c,3}$$']);
            legend_set_e(5) = fplot(axs(5), @(Re) (NB3.powerfit.b).*(Re).^(NB3.powerfit.m), [71 10000],'-', 'Color', NB3.color,'Linewidth', 1, 'DisplayName', [NB3.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_f(1) = plot(axs(6), UB1.Re_s, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_f(2) = plot(axs(6), UB1.Re_s_TV, UB1.G_TV, UB1.specs,'Color', [0 0 0], 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName',   [UB1.label ', transitioned']);
            legend_set_f(3) = plot(axs(6), UB1.Re_sc2, UB1.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*UB1.LW, 'MarkerSize', 2*UB1.MS, 'DisplayName',      [UB1.label ', $$Re_{s,c,2}$$']);
            legend_set_f(4) = plot(axs(6), UB1.Re_sc3, UB1.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*UB1.LW, 'MarkerSize', 4*UB1.MS, 'DisplayName',      [UB1.label ', $$Re_{s,c,3}$$']);
            legend_set_f(5) = fplot(axs(6), @(Re) (UB1.powerfit.b).*(Re).^(UB1.powerfit.m), [71 10000],'-', 'Color', UB1.color,'Linewidth', 1, 'DisplayName', [UB1.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_g(1) = plot(axs(7), UB2.Re_s, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
            legend_set_g(2) = plot(axs(7), UB2.Re_s_TV, UB2.G_TV, UB2.specs,'Color', [0 0 0], 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName',   [UB2.label ', transitioned']);
            legend_set_g(3) = plot(axs(7), UB2.Re_sc2, UB2.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*UB2.LW, 'MarkerSize', 2*UB2.MS, 'DisplayName',      [UB2.label ', $$Re_{s,c,2}$$']);
            legend_set_g(4) = plot(axs(7), UB2.Re_sc3, UB2.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*UB2.LW, 'MarkerSize', 4*UB2.MS, 'DisplayName',      [UB2.label ', $$Re_{s,c,3}$$']);
            legend_set_g(5) = fplot(axs(7), @(Re) (UB2.powerfit.b).*(Re).^(UB2.powerfit.m), [71 10000],'-', 'Color', UB2.color,'Linewidth', 1, 'DisplayName', [UB2.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_h(1) = plot(axs(8), XB1.Re_s, XB1.G, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_h(2) = plot(axs(8), XB1.Re_s_TV, XB1.G_TV, XB1.specs,'Color', [0 0 0], 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName',   [XB1.label ', transitioned']);
            legend_set_h(3) = plot(axs(8), XB1.Re_sc2, XB1.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*XB1.LW, 'MarkerSize', 2*XB1.MS, 'DisplayName',      [XB1.label ', $$Re_{s,c,2}$$']);
            legend_set_h(4) = plot(axs(8), XB1.Re_sc3, XB1.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*XB1.LW, 'MarkerSize', 4*XB1.MS, 'DisplayName',      [XB1.label ', $$Re_{s,c,3}$$']);
            legend_set_h(5) = fplot(axs(8), @(Re) (XB1.powerfit.b).*(Re).^(XB1.powerfit.m), [71 10000],'-', 'Color', XB1.color,'Linewidth', 1, 'DisplayName', 'XB1 $$\beta Re_s^{\alpha}$$');

            legend_set_i(1) = plot(axs(9), XB2.Re_s, XB2.G, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);
            legend_set_i(2) = plot(axs(9), XB2.Re_s_TV, XB2.G_TV, XB2.specs,'Color', [0 0 0], 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName',   [XB2.label ', transitioned']);
            legend_set_i(3) = plot(axs(9), XB2.Re_sc2, XB2.G_c2, 'p','Color', [0 0 0], 'LineWidth', 2*XB2.LW, 'MarkerSize', 2*XB2.MS, 'DisplayName',      [XB2.label ', $$Re_{s,c,2}$$']);
            legend_set_i(4) = plot(axs(9), XB2.Re_sc3, XB2.G_c3, '|','Color', [0 0 0], 'LineWidth', 3*XB2.LW, 'MarkerSize', 4*XB2.MS, 'DisplayName',      [XB2.label ', $$Re_{s,c,3}$$']);
            legend_set_i(5) = fplot(axs(9), @(Re) (XB2.powerfit.b).*(Re).^(XB2.powerfit.m), [71 10000],'-', 'Color', XB2.color,'Linewidth', 1, 'DisplayName', [XB2.label ' $$\beta Re_s^{\alpha}$$']);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:3:7), '$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(7:9), '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(3), legend_set_c,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(4), legend_set_d,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(5), legend_set_e,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(6), legend_set_f,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(7), legend_set_g,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(8), legend_set_h,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);
            legend(axs(9), legend_set_i,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns',1);

            axis(axs,obj.Res_G_range)

            fig_out = AYfig_;
        end
        function fig_out = ALL_G_vs_Res_fits(obj, AYfig_, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall)
            AYfig_.init_tiles([2,3]);
            axs = AYfig_.ax_tile;
            hold(axs, 'on');
            box(axs,'on');

            for i = 1:6
                fplot(axs(i), @(Re) obj.G_obs_Res_slope*(Re), [1e-2 70],'--', 'Color', [0 0 0],'Linewidth', 1, 'DisplayName', '$$ \frac{2 \pi r_i r_o}{(r_o-r_i)^2} Re_s $$')
            end

            %% handling pure fluid plot component
            legend_set_a(1) = plot(axs(1), PF1.Re_s, PF1.G, PF1.specs,'Color', PF1.color, 'LineWidth', PF1.LW, 'MarkerSize', PF1.MS, 'DisplayName', PF1.label);
            legend_set_a(2) = fplot(axs(1), @(Re) (PF1.powerfit.b).*(Re).^(PF1.powerfit.m), [71 10000],'-', 'Color', PF1.color,'Linewidth', 1, 'DisplayName', [PF1.label ' $$\beta Re_s^{\alpha}$$']);
            legend_set_a(3) = plot(axs(1), PFR.Re_s, PFR.G, PFR.specs,'Color', PFR.color, 'LineWidth', PFR.LW, 'MarkerSize', PFR.MS, 'DisplayName', PFR.label);
            legend_set_a(4) = fplot(axs(1), @(Re) (PFR.powerfit.b).*(Re).^(PFR.powerfit.m), [71 10000],'-', 'Color', PFR.color,'Linewidth', 1, 'DisplayName', [PFR.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_b(1) = plot(axs(2), NB1.Re_s, NB1.G, NB1.specs,'Color', NB1.color, 'LineWidth', NB1.LW, 'MarkerSize', NB1.MS, 'DisplayName', NB1.label);
            legend_set_b(2) = plot(axs(2), NB2.Re_s, NB2.G, NB2.specs,'Color', NB2.color, 'LineWidth', NB2.LW, 'MarkerSize', NB2.MS, 'DisplayName', NB2.label);
            legend_set_b(3) = plot(axs(2), NB3.Re_s, NB3.G, NB3.specs,'Color', NB3.color, 'LineWidth', NB3.LW, 'MarkerSize', NB3.MS, 'DisplayName', NB3.label);
            legend_set_b(4) = fplot(axs(2), @(Re) (NBall.powerfit.b).*(Re).^(NBall.powerfit.m), [71 10000],'-', 'Color', NBall.color,'Linewidth', 1, 'DisplayName', [NBall.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_c(1) = plot(axs(3), UB1.Re_s, UB1.G, UB1.specs,'Color', UB1.color, 'LineWidth', UB1.LW, 'MarkerSize', UB1.MS, 'DisplayName', UB1.label);
            legend_set_c(2) = plot(axs(3), UB2.Re_s, UB2.G, UB2.specs,'Color', UB2.color, 'LineWidth', UB2.LW, 'MarkerSize', UB2.MS, 'DisplayName', UB2.label);
            legend_set_c(3) = fplot(axs(3), @(Re) (UBall.powerfit.b).*(Re).^(UBall.powerfit.m), [71 10000],'-', 'Color', UBall.color,'Linewidth', 1, 'DisplayName', [UBall.label ' $$\beta Re_s^{\alpha}$$']);

            legend_set_d(1) = plot(axs(4), XB1.Re_s, XB1.G, XB1.specs,'Color', XB1.color, 'LineWidth', XB1.LW, 'MarkerSize', XB1.MS, 'DisplayName', XB1.label);
            legend_set_d(2) = plot(axs(4), XB2.Re_s, XB2.G, XB2.specs,'Color', XB2.color, 'LineWidth', XB2.LW, 'MarkerSize', XB2.MS, 'DisplayName', XB2.label);
            legend_set_d(3) = fplot(axs(4), @(Re) (XBall.powerfit.b).*(Re).^(XBall.powerfit.m), [71 10000],'-', 'Color', XBall.color,'Linewidth', 1, 'DisplayName', [XBall.label ' $$\beta Re_s^{\alpha}$$']);

            for i=1:length(FB1.exp)
              legend_set_e(i) = plot(axs(5), FB1.exp(i).Re_s, FB1.exp(i).G, FB1.specs, 'Color', FB1.exp(i).color, 'LineWidth', FB1.LW, 'MarkerSize', FB1.MS, 'DisplayName', FB1.exp(i).label);
            end
            legend_set_e(i+1) = fplot(axs(5), @(Re) (FB1.powerfit.b).*(Re).^(FB1.powerfit.m), [71 10000],'-', 'Color', FB1.color,'Linewidth', 1, 'DisplayName', [FB1.label ' $$\beta Re_s^{\alpha}$$']);

            for i=1:length(FB2.exp)
              legend_set_f(i) = plot(axs(6), FB2.exp(i).Re_s, FB2.exp(i).G, FB2.specs, 'Color', FB2.exp(i).color,'LineWidth', FB2.LW, 'MarkerSize', FB2.MS, 'DisplayName', FB2.exp(i).label);
            end
            legend_set_f(i+1) = fplot(axs(6), @(Re) (FB2.powerfit.b).*(Re).^(FB2.powerfit.m), [71 10000],'-', 'Color', FB2.color,'Linewidth', 1, 'DisplayName', [FB2.label ' $$\beta Re_s^{\alpha}$$']);

            set(axs,'YScale', 'log', 'XScale', 'log');

            ylabel(axs(1:3:4), '$$G$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
            xlabel(axs(4:6), '$$Re_s$$', 'Interpreter', 'LaTeX','FontSize',12)

            legend(axs(1), legend_set_a,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(2), legend_set_b,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(3), legend_set_c,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(4), legend_set_d,'Location', 'NorthWest', 'Interpreter', 'Latex');
            legend(axs(5), legend_set_e,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);
            legend(axs(6), legend_set_f,'Location', 'NorthWest', 'Interpreter', 'Latex', 'NumColumns', 2);

            axis(axs,obj.Res_G_range)

            fig_out = AYfig_;
        end
    end
end
