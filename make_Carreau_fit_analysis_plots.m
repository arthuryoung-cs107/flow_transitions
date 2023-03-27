close all

tfp = thinning_fluid_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%%%%%%%%% figure 1
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_fits_tau_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Carreau_fits_tau_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_Carreau_fits_spread(FB_Carreau_fits_tau_vs_omegai_linear_AYfig,[FB1 FB2],false,false);

%%%%%%%%%%%%%%% figure 2
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_fits_tau_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Carreau_fits_tau_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_Carreau_fits_spread(FB_Carreau_fits_tau_vs_omegai_log_AYfig,[FB1 FB2],false,true);

%%%%%%%%%%%%%%% figure 3
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_fits_muinf_tau_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Carreau_fits_muinf_tau_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_Carreau_fits_spread(FB_Carreau_fits_muinf_tau_vs_omegai_linear_AYfig,[FB1 FB2],true,false);

%%%%%%%%%%%%%%% figure 4
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_fits_muinf_tau_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Carreau_fits_muinf_tau_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_Carreau_fits_spread(FB_Carreau_fits_muinf_tau_vs_omegai_log_AYfig,[FB1 FB2],true,true);

%%%%%%%%%%%%%%% figure 5
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_UB_XB_FB_G_vs_Reb'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(2, :), tfp.dim222_short];};
ALL_G_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.ALL_G_vs_Reb(ALL_G_vs_Reb_AYfig, [FB1, FB2], {PF1;PFR}, NBall, UBall, XBall);

%%%%%%%%%%%%%%% figure 6
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_par_vs_q_old'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(7, :), tfp.dim222];};
FB_Carreau_par_vs_q_old_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Carreau_fluid_params(FB_Carreau_par_vs_q_old_AYfig, [FB1 FB2], false);

%%%%%%%%%%%%%%% figure 7
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_par_vs_q_new'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(10, :), tfp.dim222];};
FB_Carreau_par_vs_q_new_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Carreau_fluid_params(FB_Carreau_par_vs_q_new_AYfig, [FB1 FB2], true);

%%%%%%%%%%%%%%% figure 8
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_mueffi_vs_omegai_old'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(7, :), tfp.dim2];};
FB_Carreau_mueffi_vs_omegai_old_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Carreau_mueffi_vs_omegai(FB_Carreau_mueffi_vs_omegai_old_AYfig, [FB1 FB2], false);

%%%%%%%%%%%%%%% figure 9
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_mueffi_vs_omegai_new'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(10, :), tfp.dim2];};
FB_Carreau_mueffi_vs_omegai_new_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Carreau_mueffi_vs_omegai(FB_Carreau_mueffi_vs_omegai_new_AYfig, [FB1 FB2], true);

%%%%%%%%%%%%%%% figure 10
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_par_vs_q_compact'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(12, :), tfp.dim2];};
FB_Carreau_par_vs_q_compact_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Carreau_par_vs_q_compact(FB_Carreau_par_vs_q_compact_AYfig, FB1, FB2);

tfp.write_figures(figs, save_dir, save_type);
