close all

bp = bingham_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_T_vs_rpm_tau_vs_omega'; 'Renderer', 'painters'; 'Position', bp.posdimfull;};
FB1_FB2_T_vs_rpm_tau_vs_omega_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = bp.FB1_FB2_T_vs_rpm_tau_vs_omega(FB1_FB2_T_vs_rpm_tau_vs_omega_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_bingham_visuals'; 'Renderer', 'painters'; 'Position', bp.posdimfull;};
FB1_FB2_bingham_visuals_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = bp.FB1_FB2_bingham_visuals(FB1_FB2_bingham_visuals_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_radial_shear'; 'Renderer', 'painters'; 'Position', bp.posdimfull;};
FB1_FB2_radial_shear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = bp.FB1_FB2_radial_shear(FB1_FB2_radial_shear_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_analytical_vs_empirical_shear'; 'Renderer', 'painters'; 'Position', bp.posdimfull;};
FB1_FB2_analytical_vs_empirical_shear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = bp.FB1_FB2_analytical_vs_empirical_shear(FB1_FB2_analytical_vs_empirical_shear_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_tau_vs_omega_tauy_fit_range'; 'Renderer', 'painters'; 'Position', bp.posdimfull;};
FB1_tau_vs_omega_tauy_fit_range_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = bp.FB_tau_vs_omega_tauy_fit_range(FB1_tau_vs_omega_tauy_fit_range_AYfig, FB1, [2,5]);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  6  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB2_tau_vs_omega_tauy_fit_range'; 'Renderer', 'painters'; 'Position', bp.posdimfull;};
FB2_tau_vs_omega_tauy_fit_range_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = bp.FB_tau_vs_omega_tauy_fit_range(FB2_tau_vs_omega_tauy_fit_range_AYfig, FB2,[3,5]);



bp.write_figures(figs, save_dir, save_type);
