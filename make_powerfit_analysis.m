close all

tfp = thinning_fluid_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%%%%%%%%% figure 1
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'ALL_G_vs_Reb_powerfits_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
ALL_G_vs_Reb_powerfits_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.ALL_G_vs_Reb_powerfits_spread(ALL_tau_vs_omegai_linear_AYfig,{FB1;FB2;EC000;EC050;EC075;EC100},{PF1;PFR},NBall,UBall,XBall,false);

%%%%%%%%%%%%%%% figure 2
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'ALL_G_vs_Reb_powerfits_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
ALL_G_vs_Reb_powerfits_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.ALL_G_vs_Reb_powerfits_spread(ALL_G_vs_Reb_powerfits_log_AYfig,{FB1;FB2;EC000;EC050;EC075;EC100},{PF1;PFR},NBall,UBall,XBall,true);

tfp.write_figures(figs, save_dir, save_type);
