close all

tp = transition_plots(write_figs, write_all_figs, figs_to_write);
tfp = thinning_fluid_plots(write_figs, write_all_figs, figs_to_write);
bp = bingham_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'ALL_G_vs_Res_fits'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% ALL_G_vs_Res_fits_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.ALL_G_vs_Res_fits(ALL_G_vs_Res_fits_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_FB2_T_alphaT_vs_omega'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% FB1_FB2_T_alphaT_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.FB1_FB2_T_alphaT_vs_omega(FB1_FB2_T_alphaT_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_FB2_Re_sc_vs_q'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% FB1_FB2_Re_sc_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.FB1_FB2_Re_sc_vs_q(FB1_FB2_Re_sc_vs_q_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'PF_NUXB_G_vs_Res_trans'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% PF_NUXB_G_vs_Res_trans_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.PF_NUXB_G_vs_Res_trans(PF_NUXB_G_vs_Res_trans_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, NBall, UBall, XBall);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_FB2_G_vs_Res_trans'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% FB1_FB2_G_vs_Res_trans_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.FB1_FB2_G_vs_Res_trans(FB1_FB2_G_vs_Res_trans_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  6  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_FB2_T_vs_omega_trans'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% FB1_FB2_T_vs_omega_trans_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.FB1_FB2_T_vs_omega_trans(FB1_FB2_T_vs_omega_trans_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  7  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_alphaT_vs_omega'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
FB1_FB2_alphaT_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB1_FB2_alphaT_vs_omega(FB1_FB2_alphaT_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  8  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_FB2_alpha_comp_vs_shear'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% FB1_FB2_alpha_comp_vs_shear_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.FB1_FB2_alpha_comp_vs_shear(FB1_FB2_alpha_comp_vs_shear_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  9  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'PF_NUXB_alpha_vs_Res_trans'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% PF_NUXB_alpha_vs_Res_trans_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.PF_NUXB_alpha_vs_Res_trans(PF_NUXB_alpha_vs_Res_trans_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, NBall, UBall, XBall);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  10  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_FB2_alphaT_vs_omega'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
% FB1_FB2_alphaT_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.FB1_FB2_alphaT_vs_omega(FB1_FB2_alphaT_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  11  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB_T_vs_omega'; 'Renderer', 'painters'; 'Position', [tp.pos_spread(12, :), tp.dim2];};
% FB_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tp.FB_T_vs_omega(FB_T_vs_omega_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  12  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_G_vs_Res_spread'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
FB_G_vs_Res_spread_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB_G_vs_Res_spread(FB_G_vs_Res_spread_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  13  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_dimensional_regime_plot'; 'Renderer', 'painters'; 'Position', [tp.pos_spread(1, :), tp.dim2];};
FB_dimensional_regime_plot_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB_dimensional_regime_plot(FB_dimensional_regime_plot_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  14  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_dimensionless_regime_plot'; 'Renderer', 'painters'; 'Position', [tp.pos_spread(7, :), tp.dim2];};
FB_dimensionless_regime_plot_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tp.FB_dimensionless_regime_plot(FB_dimensionless_regime_plot_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  15  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Bingham_fluid_full_fits'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
FB1_FB2_Bingham_fluid_fits_full_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = bp.FB1_FB2_Bingham_fluid_fits_full(FB1_FB2_Bingham_fluid_fits_full_AYfig, FB1,FB2,true);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  16  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Bingham_fluid_params'; 'Renderer', 'painters'; 'Position', tfp.posdim_SW;};
FB1_FB2_Bingham_fluid_params_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fluid_params(FB1_FB2_Bingham_fluid_params_AYfig, FB1,FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  17  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Bingham_G_vs_Res'; 'Renderer', 'painters'; 'Position', [tp.pos_spread(6, :), tp.dim2];};
FB1_FB2_Bingham_G_vs_Res_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tp.FB1_FB2_Bingham_G_vs_Res(FB1_FB2_Bingham_G_vs_Res_AYfig, FB1,FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  18  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_appmu_vs_omega_Bingham_spread'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
FB_appmu_vs_omega_Bingham_spread_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = bp.FB_appmu_vs_omega_Bingham_spread(FB_appmu_vs_omega_Bingham_spread_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

tp.write_figures(figs, save_dir, save_type);
