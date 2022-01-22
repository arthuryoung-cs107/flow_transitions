clear
path = '/Users/Arthur/Documents/anton_parr/';
filename = 'Arthur_10_17_2019 5_47 PM_XB109_40frac_laminar_turbulent_run';

raw_XB2_1 = rheometer_excel_parser2018([path filename], 1);
raw_XB2_2 = rheometer_excel_parser2018([path filename], 2);
raw_XB2_3 = rheometer_excel_parser2018([path filename], 3);
raw_XB2_4 = rheometer_excel_parser2018([path filename], 4);
raw_XB2_5 = rheometer_excel_parser2018([path filename], 5);

raw_XB2_structure = {raw_XB2_1; raw_XB2_2; raw_XB2_3; raw_XB2_4; raw_XB2_5;};

save XB2_rawdata.mat
