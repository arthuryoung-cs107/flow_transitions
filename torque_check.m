clear
close all
run figure_properties.m

% % % 2018 DATA, re-run using new processing class
load ./raw_data_structures/FB1_rawdata.mat

name_list_113micron = {'FB 113 microns, Q = 1.00 l/min'; 'FB 113 microns, Q = 1.20 l/min' ; 'FB 113 microns, Q = 1.40 l/min' ; 'FB 113 microns, Q = 1.60 l/min' ; 'FB 113 microns, Q = 1.80 l/min' ; 'FB 113 microns, Q = 2.00 l/min' ; 'FB 113 microns, Q = 2.20 l/min' ; 'FB 113 microns, Q = 2.40 l/min' ; 'FB 113 microns, Q = 2.60 l/min' ; 'FB 113 microns, Q = 2.80 l/min'};

phi_vec_113microns = [0.5541; 0.5541; 0.5541; 0.5541; 0.5524; 0.5503; 0.5482; 0.5460; 0.5439; 0.5418];

glass113_all = glass113.empty(size(name_list_113micron, 1), 0);
for i = 1:size(name_list_113micron, 1)
  glass113_all(i) = glass113(name_list_113micron{i}, grey15(4+i, :), phi_vec_113microns(i));
end
FB1 = FB_experiment(glass113_all, 'FB1', grey15(1, :), '*', '113 micron glass beads');
FB1.MS = 4;
FB1.LW = 0.5;
FB1.process_raws(raw_FB1_structure)

% % % 2021 DATA, FROM ABHI
load ./raw_data_structures/FB2_rawdata.mat;
name_list_49micron = {'FB 49 microns, Q = 0.00 l/min'; 'FB 49 microns, Q = 0.05 l/min'; 'FB 49 microns, Q = 0.10 l/min'; 'FB 49 microns, Q = 0.15 l/min' ; 'FB 49 microns, Q = 0.30 l/min' ; 'FB 49 microns, Q = 0.50 l/min' ; 'FB 49 microns, Q = 0.75 l/min' ; 'FB 49 microns, Q = 1.00 l/min' ; 'FB 49 microns, Q = 1.25 l/min' ; 'FB 49 microns, Q = 1.50 l/min' ; 'FB 49 microns, Q = 2.00 l/min' ; 'FB 49 microns, Q = 2.50 l/min' ; 'FB 49 microns, Q = 3.00 l/min'};

glass49_all = glass49.empty(size(name_list_49micron, 1), 0);
for i = 1:size(name_list_49micron, 1)
  glass49_all(i) = glass49(name_list_49micron{i}, grey15(i, :));
end
FB2 = FB_experiment(glass49_all, 'FB1', grey15(1, :), 'x', '49 micron glass beads');
FB2.LW = 0.5;
FB2.LW_L = 2;
FB2.process_raws(raw_FB2_structure)
