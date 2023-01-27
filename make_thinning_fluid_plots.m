close all

tfp = thinning_fluid_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_tau_vs_omega_linscale'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% FB1_tau_vs_omega_linscale_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.FB_tau_vs_omega_linscale(FB1_tau_vs_omega_linscale_AYfig, FB1, [2,5]);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB2_tau_vs_omega_linscale'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% FB2_tau_vs_omega_linscale_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.FB_tau_vs_omega_linscale(FB2_tau_vs_omega_linscale_AYfig, FB2,[3,5]);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_power_fluid_fits'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% FB1_power_fluid_fits_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.FB_power_fluid_fits(FB1_power_fluid_fits_AYfig, FB1, [2,5]);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB2_power_fluid_fits'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% FB2_power_fluid_fits_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.FB_power_fluid_fits(FB2_power_fluid_fits_AYfig, FB2,[3,5]);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_fit_results'; 'Renderer', 'painters'; 'Position', tfp.pos_top_row;};
FB1_FB2_Carreau_fit_results_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_fit_results(FB1_FB2_Carreau_fit_results_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  6  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_viscosity'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB1_FB2_Carreau_viscosity_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_viscosity(FB1_FB2_Carreau_viscosity_AYfig, FB1, FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  7  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_fluid_fits'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB1_FB2_Carreau_fluid_fits_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_fluid_fits(FB1_FB2_Carreau_fluid_fits_AYfig, FB1,FB2);

% %%%%%%%% --------------------------------------------------------------------------------------------
% %%%%%%%% -----------------------------------------  9  ---------------------------------------------
% %%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_fluid_G_vs_Res'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB1_FB2_Carreau_fluid_G_vs_Res_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_fluid_G_vs_Res(FB1_FB2_Carreau_fluid_G_vs_Res_AYfig, FB1, FB2);

% %%%%%%%% --------------------------------------------------------------------------------------------
% %%%%%%%% -----------------------------------------  10  ---------------------------------------------
% %%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_fluid_Grat_vs_Res'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB1_FB2_Carreau_fluid_Grat_vs_Res_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_fluid_Grat_vs_Res(FB1_FB2_Carreau_fluid_Grat_vs_Res_AYfig, FB1, FB2);

% %%%%%%%% --------------------------------------------------------------------------------------------
% %%%%%%%% -----------------------------------------  11  ---------------------------------------------
% %%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_fluid_Res_vs_omegai'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB1_FB2_Carreau_fluid_Res_vs_omegai_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_fluid_Res_vs_omegai(FB1_FB2_Carreau_fluid_Res_vs_omegai_AYfig, FB1, FB2);

% %%%%%%%% --------------------------------------------------------------------------------------------
% %%%%%%%% -----------------------------------------  12  ---------------------------------------------
% %%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'ALL_G_vs_Res_CarreauFB'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(6, :), tfp.dim222];};
ALL_G_vs_Res_CarreauFB_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.ALL_G_vs_Res_CarreauFB(ALL_G_vs_Res_CarreauFB_AYfig, PF1, PFR, NB1, NB2, NB3, UB1, UB2, XB1, XB2, FB1, FB2, NBall, UBall, XBall);

% %%%%%%%% --------------------------------------------------------------------------------------------
% %%%%%%%% -----------------------------------------  13  ---------------------------------------------
% %%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB_Carreau_err_vs_omegai'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB_Carreau_err_vs_omegai_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB_Carreau_err_vs_omegai(FB_Carreau_err_vs_omegai_AYfig, FB1, FB2);


%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  13  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

% fig_num = fig_num + 1;
% fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_tau_vs_gamma'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
% FB1_FB2_Carreau_tau_vs_gamma_AYfig = AYfig(fig_specs{fig_num}, true);
% figs(fig_num) = tfp.FB1_FB2_Carreau_tau_vs_gamma(FB1_FB2_Carreau_tau_vs_gamma_AYfig, FB1,FB2);
