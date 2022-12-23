close all
mp = main_plots(write_figs, write_all_figs, figs_to_write);

fig_specs{1} = {'Name', 'PF_NB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(1, :), mp.dim2];};
fig_specs{2} = {'Name', 'PF_NB_Grat_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(2, :), mp.dim2];};
fig_specs{3} = {'Name', 'PF_NB_cf_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(3, :), mp.dim21];};
fig_specs{4} = {'Name', 'UB_XB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(4, :), mp.dim2];};
fig_specs{5} = {'Name', 'UB_XB_cf_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(5, :), mp.dim21];};
fig_specs{6} = {'Name', 'PF_NB_UB_XB_FB_G_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(6, :), mp.dim222];};
fig_specs{7} = {'Name', 'FB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(7, :), mp.dim2];};
fig_specs{8} = {'Name', 'FB_muapp_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(8, :), mp.dim2];};
fig_specs{9} = {'Name', 'FB_mup_tauy_vs_q_22'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(9, :), mp.dim22];};
fig_specs{10} = {'Name', 'FB_Carreau_par_vs_q'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(10, :), mp.dim32];};
fig_specs{11} = {'Name', 'FB_dimensional_regime_plot'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(1, :), mp.dim2];};

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------
fig_num = fig_num + 1;
PF_NB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PF_NB_T_vs_omega(PF_NB_T_vs_omega_AYfig, PF1, PFR, NB1, NB2, NB3);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
PF_NB_Grat_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PF_NB_Grat_vs_Res(PF_NB_Grat_vs_Res_AYfig, PF1, PFR, NB1, NB2, NB3, RM10, RM20);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
PF_NB_cf_alpha_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PF_NB_cf_alpha_vs_Res(PF_NB_cf_alpha_vs_Res_AYfig, PF1, PFR, NB1, NB2, NB3, RV, LSa, LSb, RK);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
UB_XB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.UB_XB_T_vs_omega(UB_XB_T_vs_omega_AYfig, UB1, UB2, XB1, XB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
UB_XB_cf_alpha_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.UB_XB_cf_alpha_vs_Res(UB_XB_cf_alpha_vs_Res_AYfig, UB1, UB2, XB1, XB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   6  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
ALL_G_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.ALL_G_vs_Res(ALL_G_vs_Res_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   7  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_T_vs_omega(FB_T_vs_omega_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   8  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB_appmu_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_appmu_vs_omega(FB_appmu_vs_omega_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------   9  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB_tauyrat_mup_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_tauyrat_mup_vs_q(FB_tauyrat_mup_vs_q_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  10  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB_Carreau_par_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_Carreau_par_vs_q(FB_Carreau_par_vs_q_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  11  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB_dimensional_regime_plot_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_dimensional_regime_plot(FB_dimensional_regime_plot_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% ----------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


ftt = flow_transitions_tables();
ftt.write_FB_fit_results(FB1,FB2,'FB_Bingham_Carreau_results')

mp.write_figures(figs, save_dir, save_type);
