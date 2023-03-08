clear

run load_data.m

% output information
write_figs = true;
write_all_figs = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';

run make_powerfit_analysis.m
