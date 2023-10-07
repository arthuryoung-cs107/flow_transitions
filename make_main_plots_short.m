close all
mp = main_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%%%%%%%%% figure 4
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(1, :), mp.dim2];};
FB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_T_vs_omega(FB_T_vs_omega_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%%%%%%%%% figure 5
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_dimensional_regime_map'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(2, :), mp.dim2];};
FB_dimensional_regime_plot_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_dimensional_regime_plot(FB_dimensional_regime_plot_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 6
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_phi_vs_omegai_vs_q'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(3, :), mp.dim1];};
FB1_phi_vs_omegai_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_phi_vs_omegai_vs_q(FB1_phi_vs_omegai_vs_q_AYfig, FB1_phi_exp);

%%%%%%%%%%%%%%% figure 9
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_mup_tauy_vs_q_compact'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(4, :), mp.dim2];};
FB_tauyrat_mup_vs_q_compact_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_tauyrat_mup_vs_q_compact(FB_tauyrat_mup_vs_q_compact_AYfig, FB1, FB2, {EC000;EC050;EC075;EC100});

%%%%%%%%%%%%%%% figure 10
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_par_vs_q_compact'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(5, :), mp.dim2];};
FB_Carreau_par_vs_q_compact_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_Carreau_par_vs_q_compact(FB_Carreau_par_vs_q_compact_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 11
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_UB_XB_FB_G_vs_Reb'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(7, :) mp.dim222_short];};
ALL_G_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.ALL_G_vs_Reb(ALL_G_vs_Reb_AYfig, [FB1, FB2], {PF1;PFR}, NBall, UBall, XBall,{EC000;EC050;EC075;EC100});

%%%%%%%%%%%%%%% figure 12
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NUXB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(8, :), mp.dim22_med];};
PL_NUXB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PL_NUXB_T_vs_omega(PL_NUXB_T_vs_omega_AYfig, {PF1;PFR}, NBall, UBall, XBall);

%%%%%%%%%%%%%%% figure 13
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NUXB_Grat_vs_Res'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(9, :), mp.dim22_med];};
PL_NUXB_Grat_vs_Res_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.PL_NUXB_Grat_vs_Res(PL_NUXB_Grat_vs_Res_AYfig, {PF1;PFR}, NBall, UBall, XBall, [RM10 RM20]);

%%%%%%%%%%%%%%% figure 16
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UXB_FB_Grat_Rebstar_compact'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(10,:) mp.dim1];};
UXB_FB_Grat_Rebstar_compact_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = mp.UXB_FB_Grat_Rebstar_compact(UXB_FB_Grat_Rebstar_compact_AYfig,UBall,XBall,[FB1 FB2],{EC000;EC050;EC075;EC100});

%%%%%%%%%%%%%%% graphical abstract
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_dimensional_regime_map'; 'Renderer', 'painters'; 'Units', 'centimeters'; 'Position', mp.posdim1_graphical_abstract;};
FB1_dimensional_regime_plot_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = mp.FB1_dimensional_regime_plot(FB1_dimensional_regime_plot_AYfig, FB1);

%%%%%%%%%%%%%%% FB1 solid fraction, sliced
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_phi_vs_omegai_vs_q_slices'; 'Renderer', 'painters'; 'Position', [mp.pos_spread(11, :), mp.dim22_tall];};
FB1_phi_vs_omegai_vs_q_slices_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = mp.FB_phi_vs_omegai_vs_q_slices(FB1_phi_vs_omegai_vs_q_slices_AYfig, FB1_phi_exp);


%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% ----------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% mp.write_figures(figs, save_dir, save_type);
