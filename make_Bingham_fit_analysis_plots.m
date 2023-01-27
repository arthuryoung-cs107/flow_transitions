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

% %%%%%%%%%%%%%%% figure 3
% %%%%%%%%%%%%%%%%%%%%%%%%
% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'ALL_dtdo_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% ALL_dtdo_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.ALL_dtdo_vs_omegai_spread(ALL_dtdo_vs_omegai_linear_AYfig,[FB1,FB2],{PF1;PFR},NBall,UBall,XBall,false);
%
% %%%%%%%%%%%%%%% figure 4
% %%%%%%%%%%%%%%%%%%%%%%%%
% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'ALL_dtdo_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% ALL_dtdo_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.ALL_dtdo_vs_omegai_spread(ALL_dtdo_vs_omegai_log_AYfig,[FB1,FB2],{PF1;PFR},NBall,UBall,XBall,true);
%
% %%%%%%%%%%%%%%% figure 5
% %%%%%%%%%%%%%%%%%%%%%%%%
% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'ALL_d2tdo2_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% ALL_d2tdo2_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.ALL_d2tdo2_vs_omegai_spread(ALL_d2tdo2_vs_omegai_linear_AYfig,[FB1,FB2],{PF1;PFR},NBall,UBall,XBall,false);
%
% %%%%%%%%%%%%%%% figure 6
% %%%%%%%%%%%%%%%%%%%%%%%%
% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'ALL_d2tdo2_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% ALL_d2tdo2_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.ALL_d2tdo2_vs_omegai_spread(ALL_d2tdo2_vs_omegai_log_AYfig,[FB1,FB2],{PF1;PFR},NBall,UBall,XBall,true);

%%%%%%%%%%%%%%% figure 7
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_fits_tau_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Bingham_fits_tau_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fits_spread(FB_Bingham_fits_tau_vs_omegai_linear_AYfig,FB1,FB2,false);

%%%%%%%%%%%%%%% figure 8
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_fits_tau_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Bingham_fits_tau_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fits_spread(FB_Bingham_fits_tau_vs_omegai_log_AYfig,FB1,FB2,true);


%%%%%%%%%%%%%%% figure 9
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_params'; 'Renderer', 'painters'; 'Position', [1 1 tfp.dim22_tall];};
FB_Bingham_params_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fluid_params_fitcheck(FB_Bingham_params_AYfig,FB1,FB2);

%%%%%%%%%%%%%%% figure 10
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_fits_Grat_vs_Reb'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(1, :), tfp.dim32_tall];};
FB_Bingham_fits_Grat_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_Bingham_fits_Grat_vs_Reb(FB_Bingham_fits_Grat_vs_Reb_AYfig, FB1, FB2, PFR);


%%%%%%%%%%%%%%% figure 11
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'PL_NB_UB_XB_FB_G_vs_Reb'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(2, :), tfp.dim222];};
ALL_G_vs_Reb_AYfig = AYfig(fig_specs{fig_num}, false);
% figs(fig_num) = tfp.ALL_G_vs_Reb(ALL_G_vs_Reb_AYfig, [FB1, FB2], {PF1;PFR}, NBall, UBall, XBall, {EC000;EC050;EC075;EC100});
figs(fig_num) = tfp.ALL_G_vs_Reb(ALL_G_vs_Reb_AYfig, [FB1, FB2], {PF1;PFR}, NBall, UBall, XBall);

%%%%%%%%%%%%%%% figure 12
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_appmu_vs_omegai'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(3, :), tfp.dim2];};
FB_Bingham_appmu_vs_omegai_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_appmu_vs_omegai(FB_Bingham_appmu_vs_omegai_AYfig, FB1, FB2);

%%%%%%%%%%%%%%% figure 13
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_mup_tauy_vs_q'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(6, :), tfp.dim22];};
FB_tauyrat_mup_vs_q_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_tauyrat_mup_vs_q(FB_tauyrat_mup_vs_q_AYfig, FB1, FB2, EC000, EC050, EC075, EC100);

%%%%%%%%%%%%%%% figure 14
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_mup_tauy_vs_q_compact'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(12, :), tfp.dim2];};
FB_tauyrat_mup_vs_q_compact_AYfig = AYfig(fig_specs{fig_num}, false);
figs(fig_num) = tfp.FB_tauyrat_mup_vs_q_compact(FB_tauyrat_mup_vs_q_compact_AYfig, FB1, FB2, {EC000;EC050;EC075;EC100});

%%%%%%%%%%%%%%% figure 15
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_taui_rat_Bingham_vs_omegai_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_taui_rat_Bingham_vs_omegai_linear_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_taui_rat_Bingham_vs_omegai_spread (FB_taui_rat_Bingham_vs_omegai_linear_AYfig,FB1,FB2,false);

%%%%%%%%%%%%%%% figure 16
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_taui_rat_Bingham_vs_omegai_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_taui_rat_Bingham_vs_omegai_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_taui_rat_Bingham_vs_omegai_spread (FB_taui_rat_Bingham_vs_omegai_log_AYfig,FB1,FB2,true);

%%%%%%%%%%%%%%% figure 15
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_err_vs_num_linear'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Bingham_err_vs_num_linear = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_Bingham_err_vs_num_spread(FB_Bingham_err_vs_num_linear,FB1,FB2,{EC000;EC050;EC075;EC100},false);

%%%%%%%%%%%%%%% figure 16
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Bingham_err_vs_num_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Bingham_err_vs_num_log_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_Bingham_err_vs_num_spread(FB_Bingham_err_vs_num_log_AYfig,FB1,FB2,{EC000;EC050;EC075;EC100},true);

%%%%%%%%%%%%%%% figure 17
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UXB_Bingham_fit_summary_log'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
UXB_Bingham_fit_summary_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.UXB_Bingham_fit_summary(UXB_Bingham_fit_summary_AYfig,UBall,XBall,[FB1 FB2]);

%%%%%%%%%%%%%%% figure 18
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UXB_FB_taustar_vs_Gamma'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
UXB_FB_taustar_vs_Gamma_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.UXB_FB_taustar_vs_Gamma(UXB_FB_taustar_vs_Gamma_AYfig,UBall,XBall,[FB1 FB2]);

%%%%%%%%%%%%%%% figure 19
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UXB_FB_taustar_vs_Gamma_combined'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(1,:) tfp.dim2];};
UXB_FB_taustar_vs_Gamma_combined_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.UXB_FB_taustar_vs_Gamma_combined(UXB_FB_taustar_vs_Gamma_combined_AYfig,UBall,XBall,[FB1 FB2],{EC000;EC050;EC075;EC100});

%%%%%%%%%%%%%%% figure 20
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FBext_Bingham_fits_summary'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FBext_Bingham_fits_summary_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FBext_Bingham_fits_summary(FBext_Bingham_fits_summary_AYfig,{EC000; EC050; EC075; EC100; FB1; FB2});

%%%%%%%%%%%%%%% figure 21
%%%%%%%%%%%%%%%%%%%%%%%%
fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'UXB_FB_taustar_vs_Gamma_compact'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(2,:) tfp.dim1];};
UXB_FB_taustar_vs_Gamma_compact_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.UXB_FB_taustar_vs_Gamma_compact (UXB_FB_taustar_vs_Gamma_compact_AYfig,UBall,XBall,[FB1 FB2],{EC000;EC050;EC075;EC100});


FB1.fitted_Bingham_pars
FB2.fitted_Bingham_pars

for FB = [FB1 FB2]
    for i = 1:length(FB.exp)
        expi = FB.exp(i);
        tauy = expi.tau_y_Bingham;
        omega_fit = expi.omega_fit_Bingham;
        tau_fit = expi.tau_fit_Bingham;
        iearly = omega_fit<1;
        tau_early = tau_fit(iearly);
        tau_late = tau_fit(~iearly);
        delta_tau = max(tau_early) - min(tau_late);
        tau_mean = mean(tau_fit);
        fprintf('%s: tau_y = %e, mean= %e, del_tau = %e, tau_y ratio: %f, mean ratio: %f\n', expi.label, tauy, tau_mean, delta_tau, delta_tau/tauy, delta_tau/tau_mean);
    end
end

tfp.write_figures(figs, save_dir, save_type);
