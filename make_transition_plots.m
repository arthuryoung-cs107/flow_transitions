close all

tp = transition_plots(write_figs, write_all_figs, figs_to_write);

fig_specs{1} = {'Name', 'ALL_G_vs_Res_fits'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};
fig_specs{2} = {'Name', 'FB1_FB2_T_alphaT_vs_omega'; 'Renderer', 'painters'; 'Position', tp.posdimfull;};

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
%%%%%%%% -----------------------------------------  end plots  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

tp.write_figures(figs, save_dir, save_type);
