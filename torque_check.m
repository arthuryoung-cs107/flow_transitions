clear
close all
run figure_properties.m
fig_pos12 = fig_pos_gen(2, 6);

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
FB1.def_pos = fig_pos12(1, :);
FB1.process_raws(raw_FB1_structure)
FB1.inspect_torques()

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
FB2.def_pos = fig_pos12(2, :);
FB2.process_raws(raw_FB2_structure)
FB2.inspect_torques()

load ./raw_data_structures/PF1_rawdata.mat
name_list_PF1 = {'PF1, smooth, run 1'; 'PF1, smooth, run 2'; 'PF1, smooth, run 3'; 'PF1, smooth, run 4';};
PF1_all = fluid.empty(size(name_list_PF1, 1), 0);
for i = 1:size(name_list_PF1, 1)
  PF1_all(i) = fluid(name_list_PF1{i}, blue12(i, :));
end
PF1 = PF1_experiment(PF1_all, blue5, ' .');
PF1.process_raws(raw_PF1_structure);
PF1.inspect_torques()

load ./raw_data_structures/PF2_rawdata.mat
name_list_PF2 = {'PF2, smooth, run 1'; 'PF2, smooth, run 2'; 'PF2, smooth, run 3'; 'PF2, smooth, run 4'; 'PF2, smooth, run 5';};
PF2_all = fluid.empty(size(name_list_PF2, 1), 0);
for i = 1:size(name_list_PF2, 1)
  PF2_all(i) = fluid(name_list_PF2{i}, red12(i, :));
end
PF2 = PF2_experiment(PF2_all, red2, ' o');
PF2.process_raws(raw_PF2_structure);
PF2.inspect_torques()

load ./raw_data_structures/PF3_rawdata.mat
name_list_PF3 = {'PF3, rough, run 1'; 'PF3, rough, run 2'; 'PF3, rough, run 3'; 'PF3, rough, run 4'; 'PF3, rough, run 5';};
PF3_all = fluid.empty(size(name_list_PF3, 1), 0);
for i = 1:size(name_list_PF3, 1)
PF3_all(i) = fluid(name_list_PF3{i}, red12(i+length(name_list_PF2), :));
end
PF3 = PF3_experiment(PF3_all, red5, ' o');
PF3.process_raws(raw_PF3_structure);
PF3.inspect_torques()

load ./raw_data_structures/PF4_rawdata.mat
name_list_PF4 = {'PF4, smooth, run 1'; 'PF4, smooth, run 2'; 'PF4, smooth, run 3'; 'PF4, smooth, run 4'; 'PF4, smooth, run 5';};
PF4_all = fluid.empty(size(name_list_PF4, 1), 0);
for i = 1:size(name_list_PF4, 1)
  PF4_all(i) = fluid(name_list_PF4{i}, grey(i, :));
end
PF4 = PF4_experiment(PF4_all, grey5, ' o');
PF4.process_raws(raw_PF4_structure);

load ./raw_data_structures/NB1_rawdata.mat
name_list_NB1 = {'NB1, rough, run 1'; 'NB1, rough, run 2'; 'NB1, rough, run 3'; 'NB1, rough, run 4'; 'NB1, rough, run 5';};
NB1_all = fluid.empty(size(name_list_NB1, 1), 0);
for i = 1:size(name_list_NB1, 1)
  NB1_all(i) = fluid(name_list_NB1{i}, green12(i, :));
end
NB1 = NB1_experiment(NB1_all, green5, ' s');
NB1.process_raws(raw_NB1_structure);
NB1.inspect_torques()

load ./raw_data_structures/NB2_rawdata.mat
name_list_NB2 = {'NB2, rough, run 1'; 'NB2, rough, run 2'; 'NB2, rough, run 3'; 'NB2, rough, run 4'; 'NB2, rough, run 5'; 'NB2, rough, run 6';};
NB2_all = fluid.empty(size(name_list_NB2, 1), 0);
for i = 1:size(name_list_NB2, 1)
  NB2_all(i) = fluid(name_list_NB2{i}, green12(i, :));
end
NB2 = NB2_experiment(NB2_all, green2, ' d');
NB2.process_raws(raw_NB2_structure);
NB2.inspect_torques()

load ./raw_data_structures/NB3_rawdata.mat
name_list_NB3 = {'NB3, rough, run 1'; 'NB3, rough, run 2'; 'NB3, rough, run 3'; 'NB3, rough, run 4'; 'NB3, rough, run 5';};
NB3_all = fluid.empty(size(name_list_NB3, 1), 0);
for i = 1:size(name_list_NB3, 1)
  NB3_all(i) = fluid(name_list_NB3{i}, green12(i, :));
end
NB3 = NB3_experiment(NB3_all, green4, ' h');
NB3.process_raws(raw_NB3_structure);
NB3.inspect_torques()

load ./raw_data_structures/UB1_rawdata.mat
name_list_UB1 = {'UB1, rough, run 1'; 'UB1, rough, run 2'; 'UB1, rough, run 3'; 'UB1, rough, run 4'; 'UB1, rough, run 5';};
UB1_all = fluid.empty(size(name_list_UB1, 1), 0);
for i = 1:size(name_list_UB1, 1)
  UB1_all(i) = fluid(name_list_UB1{i}, blue12(i, :));
end
UB1 = UB1_experiment(UB1_all, blue2, ' v');
UB1.process_raws(raw_UB1_structure);
UB1.inspect_torques()

load ./raw_data_structures/UB2_rawdata.mat
name_list_UB2 = {'UB2, rough, run 1'; 'UB2, rough, run 2'; 'UB2, rough, run 3'; 'UB2, rough, run 4'; 'UB2, rough, run 5';};
UB2_all = fluid.empty(size(name_list_UB2, 1), 0);
for i = 1:size(name_list_UB2, 1)
  UB2_all(i) = fluid(name_list_UB2{i}, blue12(i+length(name_list_UB1), :));
end
UB2 = UB2_experiment(UB2_all, blue3, ' ^');
UB2.process_raws(raw_UB2_structure);
UB2.inspect_torques()

load ./raw_data_structures/XB1_rawdata.mat
name_list_XB1 = {'XB1, rough, run 1'; 'XB1, rough, run 2'; 'XB1, rough, run 3'; 'XB1, rough, run 4'; 'XB1, rough, run 5';};
XB1_all = fluid.empty(size(name_list_XB1, 1), 0);
for i = 1:size(name_list_XB1, 1)
  XB1_all(i) = fluid(name_list_XB1{i}, orange12(i, :));
end
XB1 = XB1_experiment(XB1_all, orange3, ' <');
XB1.process_raws(raw_XB1_structure);
XB1.inspect_torques()

load ./raw_data_structures/XB2_rawdata.mat
name_list_XB2 = {'XB2, rough, run 1'; 'XB2, rough, run 2'; 'XB2, rough, run 3'; 'XB2, rough, run 4'; 'XB2, rough, run 5';};
XB2_all = fluid.empty(size(name_list_XB2, 1), 0);
for i = 1:size(name_list_XB2, 1)
  XB2_all(i) = fluid(name_list_XB2{i}, orange12(i+length(name_list_XB1), :));
end
XB2 = XB2_experiment(XB2_all, orange1, ' >');
XB2.process_raws(raw_XB2_structure);
XB2.inspect_torques()
