close all

tfp = thinning_fluid_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_fluid_fits'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB1_FB2_Carreau_fluid_fits_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_fluid_fits_full(FB1_FB2_Carreau_fluid_fits_AYfig, FB1,FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  2  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Bingham_fluid_fits'; 'Renderer', 'painters'; 'Position', tfp.posdimfull;};
FB1_FB2_Bingham_fluid_fits_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fluid_fits_full(FB1_FB2_Bingham_fluid_fits_AYfig, FB1,FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  3  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Carreau_fluid_params'; 'Renderer', 'painters'; 'Position', tfp.pos_top_row;};
FB1_FB2_Carreau_fluid_params_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Carreau_fluid_params(FB1_FB2_Carreau_fluid_params_AYfig, FB1,FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  4  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Bingham_fluid_params'; 'Renderer', 'painters'; 'Position', tfp.pos_bottom_row;};
FB1_FB2_Bingham_fluid_params_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_fluid_params(FB1_FB2_Bingham_fluid_params_AYfig, FB1,FB2);

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  5  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
fig_specs{fig_num} = {'Name', 'FB1_FB2_Bingham_and_Carreau_fluid_params'; 'Renderer', 'painters'; 'Position', [tfp.pos_spread(6, :), tfp.dim222];};
FB1_FB2_Bingham_and_Carreau_fluid_params_AYfig = AYfig(fig_specs{fig_num}, true);
figs(fig_num) = tfp.FB1_FB2_Bingham_and_Carreau_fluid_params(FB1_FB2_Bingham_and_Carreau_fluid_params_AYfig, FB1,FB2);
