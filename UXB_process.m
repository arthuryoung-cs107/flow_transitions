run figure_properties.m
load UB1_rawdata.mat
name_list_UB1 = {'UB1, rough, run 1'; 'UB1, rough, run 2'; 'UB1, rough, run 3'; 'UB1, rough, run 4'; 'UB1, rough, run 5';};
UB1_all = fluid.empty(size(name_list_UB1, 1), 0);
for i = 1:size(name_list_UB1, 1)
  UB1_all(i) = fluid(name_list_UB1{i}, blue12(i, :));
end
UB1 = UB1_experiment(UB1_all, green5, '- v');
UB1.process_raws(raw_UB1_structure);

load UB2_rawdata.mat
name_list_UB2 = {'UB2, rough, run 1'; 'UB2, rough, run 2'; 'UB2, rough, run 3'; 'UB2, rough, run 4'; 'UB2, rough, run 5';};
UB2_all = fluid.empty(size(name_list_UB2, 1), 0);
for i = 1:size(name_list_UB2, 1)
  UB2_all(i) = fluid(name_list_UB2{i}, blue12(i+length(name_list_UB1), :));
end
UB2 = UB2_experiment(UB2_all, green5, '- v');
UB2.process_raws(raw_UB2_structure);

load XB1_rawdata.mat
name_list_XB1 = {'XB1, rough, run 1'; 'XB1, rough, run 2'; 'XB1, rough, run 3'; 'XB1, rough, run 4'; 'XB1, rough, run 5';};
XB1_all = fluid.empty(size(name_list_XB1, 1), 0);
for i = 1:size(name_list_XB1, 1)
  XB1_all(i) = fluid(name_list_XB1{i}, orange12(i, :));
end
XB1 = XB1_experiment(XB1_all, green5, '- ^');
XB1.process_raws(raw_XB1_structure);

load XB2_rawdata.mat
name_list_XB2 = {'XB2, rough, run 1'; 'XB2, rough, run 2'; 'XB2, rough, run 3'; 'XB2, rough, run 4'; 'XB2, rough, run 5';};
XB2_all = fluid.empty(size(name_list_XB2, 1), 0);
for i = 1:size(name_list_XB2, 1)
  XB2_all(i) = fluid(name_list_XB2{i}, orange12(i+length(name_list_XB1), :));
end
XB2 = XB2_experiment(XB2_all, green5, '- ^');
XB2.process_raws(raw_XB2_structure);
