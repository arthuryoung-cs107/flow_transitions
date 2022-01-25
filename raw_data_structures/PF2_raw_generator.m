clear
path = './anton_parr_data/';
filename = 'Arthur_9_30_2019 7_20 PM_NB_density_long_laminar_run';

raw_PF2_1 = rheometer_excel_parser2018([path filename], 1);
raw_PF2_2 = rheometer_excel_parser2018([path filename], 2);
raw_PF2_3 = rheometer_excel_parser2018([path filename], 3);
raw_PF2_4 = rheometer_excel_parser2018([path filename], 4);
raw_PF2_5 = rheometer_excel_parser2018([path filename], 5);

raw_PF2_structure = {raw_PF2_1; raw_PF2_2; raw_PF2_3; raw_PF2_4; raw_PF2_5;};

save PF2_rawdata.mat
