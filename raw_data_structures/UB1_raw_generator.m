clear
path = './anton_parr_data/';
filename = 'Arthur_10_21_2019 7_47 PM_UBW_45frac_laminar_turbulent_run';

raw_UB1_1 = rheometer_excel_parser2018([path filename], 1);
raw_UB1_2 = rheometer_excel_parser2018([path filename], 2);
raw_UB1_3 = rheometer_excel_parser2018([path filename], 3);
raw_UB1_4 = rheometer_excel_parser2018([path filename], 4);
raw_UB1_5 = rheometer_excel_parser2018([path filename], 5);

raw_UB1_structure = {raw_UB1_1; raw_UB1_2; raw_UB1_3; raw_UB1_4; raw_UB1_5;};

save UB1_rawdata.mat
