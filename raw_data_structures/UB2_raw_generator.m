clear
path = './anton_parr_data/';
filename = 'Arthur_10_22_2019 6_39 PM_UB102_45frac_laminar_turbulent_run.xlsx';

raw_UB2_1 = rheometer_excel_parser2018([path filename], 1);
raw_UB2_2 = rheometer_excel_parser2018([path filename], 2);
raw_UB2_3 = rheometer_excel_parser2018([path filename], 3);
raw_UB2_4 = rheometer_excel_parser2018([path filename], 4);
raw_UB2_5 = rheometer_excel_parser2018([path filename], 5);

raw_UB2_structure = {raw_UB2_1; raw_UB2_2; raw_UB2_3; raw_UB2_4; raw_UB2_5;};

save UB2_rawdata.mat
