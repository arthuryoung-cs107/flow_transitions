close all
mp = main_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%%%%%%%%% figure 1
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(1, :), mp.dim2];};
PF_NB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PF_NB_T_vs_omega(PF_NB_T_vs_omega_AYfig, PF1, PFR, NB1, NB2, NB3);

%%%%%%%%%%%%%%% figure 2
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_Grat_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(2, :), mp.dim2];};
PF_NB_Grat_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PF_NB_Grat_vs_Res(PF_NB_Grat_vs_Res_AYfig, PF1, PFR, NB1, NB2, NB3, RM10, RM20);

%%%%%%%%%%%%%%% figure 3
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_cf_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(3, :), mp.dim21];};
PF_NB_cf_alpha_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PF_NB_cf_alpha_vs_Res(PF_NB_cf_alpha_vs_Res_AYfig, PF1, PFR, NB1, NB2, NB3, RV, LSa, LSb, RK);

%%%%%%%%%%%%%%% figure 4
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UB_XB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(4, :), mp.dim2];};
UB_XB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.UB_XB_T_vs_omega(UB_XB_T_vs_omega_AYfig, UB1, UB2, XB1, XB2);

%%%%%%%%%%%%%%% figure 5
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UB_XB_cf_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(5, :), mp.dim21];};
UB_XB_cf_alpha_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.UB_XB_cf_alpha_vs_Res(UB_XB_cf_alpha_vs_Res_AYfig, UB1, UB2, XB1, XB2);

%%%%%%%%%%%%%%% figure 6
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_UB_XB_FB_G_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(6, :), mp.dim222_short];};
ALL_G_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.ALL_G_vs_Res(ALL_G_vs_Res_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall);

%%%%%%%%%%%%%%% figure 7
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(7, :), mp.dim2];};
FB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_T_vs_omega(FB_T_vs_omega_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%%%%%%%%% figure 8
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_muapp_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(8, :), mp.dim2];};
FB_appmu_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_appmu_vs_omega(FB_appmu_vs_omega_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%%%%%%%%% figure 9
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_mup_tauy_vs_q'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(9, :), mp.dim22];};
FB_tauyrat_mup_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_tauyrat_mup_vs_q(FB_tauyrat_mup_vs_q_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);
% FB_tauyrat_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = mp.FB_tauyrat_mup_vs_q(FB_tauyrat_mup_vs_q_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%%%%%%%%% figure 10
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_par_vs_q'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(10, :), mp.dim32];};
FB_Carreau_par_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_Carreau_par_vs_q(FB_Carreau_par_vs_q_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 11
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_dimensional_regime_map'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(11, :), mp.dim2];};
FB_dimensional_regime_plot_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_dimensional_regime_plot(FB_dimensional_regime_plot_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 12
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_phi_vs_omegai_vs_q'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(12, :), mp.dim1];};
FB1_phi_vs_omegai_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_phi_vs_omegai_vs_q(FB1_phi_vs_omegai_vs_q_AYfig, FB1_phi_exp);


%%%%%%%%%%%%%%% figure 13
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UXB_NB_PL_Grat_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(1, :), mp.dim2];};
UXB_NB_PL_Grat_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.UXB_NB_PL_Grat_vs_Res(UXB_NB_PL_Grat_vs_Res_AYfig,UB1,UB2,XB1,XB2,PF1,PFR,NB1,NB2,NB3);

%%%%%%%%%%%%%%% figure 14
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_UB_XB_FB_G_vs_Reb'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(12, :), mp.dim222_short];};
ALL_G_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.ALL_G_vs_Reb(ALL_G_vs_Reb_AYfig, [FB1, FB2], {PF1;PFR}, NBall, UBall, XBall);

%%%%%%%%%%%%%%% figure 15
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_fits_Grat_vs_Reb'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(3, :), mp.dim32_tall];};
FB_Bingham_fits_Grat_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_Bingham_fits_Grat_vs_Reb(FB_Bingham_fits_Grat_vs_Reb_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 16
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_mup_tauy_vs_q_compact'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(4, :), mp.dim2];};
FB_tauyrat_mup_vs_q_compact_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_tauyrat_mup_vs_q_compact(FB_tauyrat_mup_vs_q_compact_AYfig, FB1, FB2, {EC000;EC050;EC075;EC100});

%%%%%%%%%%%%%%% figure 17
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_par_vs_q_compact'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(10, :), mp.dim2];};
FB_Carreau_par_vs_q_compact_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_Carreau_par_vs_q_compact(FB_Carreau_par_vs_q_compact_AYfig, FB1, FB2);


%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% ----------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------


ftt = flow_transitions_tables();
% ftt.write_FB_fit_results(FB1,FB2,'FB_Bingham_Carreau_results')

% mp.write_figures(figs, save_dir, save_type);
