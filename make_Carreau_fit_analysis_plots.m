close all

tfp = thinning_fluid_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%%%%%%%%% figure 1
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'ALL_tau_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
ALL_tau_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.ALL_tau_vs_omegai_spread(ALL_tau_vs_omegai_linear_AYfig,[FB1,FB2],{PF1;PFR},NBall,UBall,XBall,false);

%%%%%%%%%%%%%%% figure 2
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'ALL_tau_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
ALL_tau_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.ALL_tau_vs_omegai_spread(ALL_tau_vs_omegai_log_AYfig,[FB1,FB2],{PF1;PFR},NBall,UBall,XBall,true);

%%%%%%%%%%%%%%% figure 3
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_fits_tau_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Bingham_fits_tau_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fits_spread(FB_Bingham_fits_tau_vs_omegai_linear_AYfig,FB1,FB2,false);

%%%%%%%%%%%%%%% figure 4
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_fits_tau_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Bingham_fits_tau_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fits_spread(FB_Bingham_fits_tau_vs_omegai_log_AYfig,FB1,FB2,true);


%%%%%%%%%%%%%%% figure 5
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_params'; 'Renderer', 'painters'; 'Position', [1 1 tfp.dim22_tall];};
FB_Bingham_params_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fluid_params_fitcheck(FB_Bingham_params_AYfig,FB1,FB2);

%%%%%%%%%%%%%%% figure 6
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_fits_Grat_vs_Reb'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(1, :), tfp.dim32_tall];};
FB_Bingham_fits_Grat_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Bingham_fits_Grat_vs_Reb(FB_Bingham_fits_Grat_vs_Reb_AYfig, FB1, FB2, PFR);


%%%%%%%%%%%%%%% figure 7
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_UB_XB_FB_G_vs_Reb'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(2, :), tfp.dim222];};
ALL_G_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.ALL_G_vs_Reb(ALL_G_vs_Reb_AYfig, [FB1, FB2], {PF1;PFR}, NBall, UBall, XBall);

%%%%%%%%%%%%%%% figure 8
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_appmu_vs_omegai'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(3, :), tfp.dim2];};
FB_Bingham_appmu_vs_omegai_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_appmu_vs_omegai(FB_Bingham_appmu_vs_omegai_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 9
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_mup_tauy_vs_q'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(6, :), tfp.dim22];};
FB_tauyrat_mup_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_tauyrat_mup_vs_q(FB_tauyrat_mup_vs_q_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%%%%%%%%% figure 10
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_par_vs_q_compact'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(12, :), tfp.dim2];};
FB_Carreau_par_vs_q_compact_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Carreau_par_vs_q_compact(FB_Carreau_par_vs_q_compact_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 11
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_taui_rat_Bingham_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_taui_rat_Bingham_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_taui_rat_Bingham_vs_omegai_spread (FB_taui_rat_Bingham_vs_omegai_linear_AYfig,FB1,FB2,false);

%%%%%%%%%%%%%%% figure 12
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_taui_rat_Bingham_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_taui_rat_Bingham_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_taui_rat_Bingham_vs_omegai_spread (FB_taui_rat_Bingham_vs_omegai_log_AYfig,FB1,FB2,true);

FB1.fitted_Bingham_pars
FB2.fitted_Bingham_pars
