close all

bp = bingham_plots(write_figs, write_all_figs, figs_to_write);

fig_specs{1} = {'Name', 'FB1_FB2_T_vs_omega'; 'Renderer', 'painters'; 'Position', bp.pos_top_row;};
fig_specs{2} = {'Name', 'FB1_FB2_gammai_vs_omega'; 'Renderer', 'painters'; 'Position', bp.pos_bottom_row;};
fig_specs{3} = {'Name', 'FB1_FB2_bingham_visuals'; 'Renderer', 'painters'; 'Position', bp.posdimfull;};

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_T_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = bp.FB1_FB2_T_vs_omega(FB1_FB2_T_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_gammai_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = bp.FB1_FB2_gammai_vs_omega(FB1_FB2_gammai_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_FB2_bingham_visuals_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = bp.FB1_FB2_bingham_visuals(FB1_FB2_bingham_visuals_AYfig, FB1, FB2);
