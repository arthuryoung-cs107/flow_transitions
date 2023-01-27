close all

fbrdp = FB_raw_data_plots(write_figs, write_all_figs, figs_to_write);

fig_num = 0;

%%%%%%%% --------------------------------------------------------------------------------------------
%%%%%%%% -----------------------------------------  1  ---------------------------------------------
%%%%%%%% --------------------------------------------------------------------------------------------

fig_num = fig_num + 1;
FB1_torque_vs_rpm_raw_AYfig = FB_raw_data_plots.make_fig('FB1_torque_vs_rpm_raw', FB_raw_data_plots.posdim_full);
figs(fig_num) = fbrdp.FB1_torque_vs_rpm_raw_plot(FB1_torque_vs_rpm_raw_AYfig, FB1);
