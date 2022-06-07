close all

tp = transition_plots(write_figs, write_all_figs, figs_to_write);

fig_specs{1} = {'Name', 'ALL_G_vs_Res_fits'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{2} = {'Name', 'FB1_FB2_T_alphaT_vs_omega'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{3} = {'Name', 'FB1_FB2_Re_sc_vs_q'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{4} = {'Name', 'PF_NUXB_G_vs_Res_trans'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{5} = {'Name', 'FB1_FB2_G_vs_Res_trans'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{6} = {'Name', 'FB1_FB2_T_vs_omega_trans'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{7} = {'Name', 'FB1_FB2_alphaT_vs_omega'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{8} = {'Name', 'FB1_FB2_alpha_comp_vs_shear'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
ALL_G_vs_Res_fits_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.ALL_G_vs_Res_fits(ALL_G_vs_Res_fits_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_T_alphaT_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB1_FB2_T_alphaT_vs_omega(FB1_FB2_T_alphaT_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_Re_sc_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB1_FB2_Re_sc_vs_q(FB1_FB2_Re_sc_vs_q_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
PF_NUXB_G_vs_Res_trans_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.PF_NUXB_G_vs_Res_trans(PF_NUXB_G_vs_Res_trans_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, NBall, UBall, XBall);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_G_vs_Res_trans_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB1_FB2_G_vs_Res_trans(FB1_FB2_G_vs_Res_trans_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  6  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_T_vs_omega_trans_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB1_FB2_T_vs_omega_trans(FB1_FB2_T_vs_omega_trans_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  7  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_alphaT_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB1_FB2_alphaT_vs_omega(FB1_FB2_alphaT_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  8  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_alpha_comp_vs_shear_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB1_FB2_alpha_comp_vs_shear(FB1_FB2_alpha_comp_vs_shear_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

tp.write_figures(figs, save_dir, save_type);
