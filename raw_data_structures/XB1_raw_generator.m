clear
path = '/Users/Arthur/Documents/anton_parr/';
filename = 'Arthur_10_16_2019 8_56 PM_XB107_40frac_laminar_turbulent_run.xlsx';

raw_XB1_1 = rheometer_excel_parser2018([path filename], 1);
raw_XB1_2 = rheometer_excel_parser2018([path filename], 2);
raw_XB1_3 = rheometer_excel_parser2018([path filename], 3);
raw_XB1_4 = rheometer_excel_parser2018([path filename], 4);
raw_XB1_5 = rheometer_excel_parser2018([path filename], 5);

raw_XB1_structure = {raw_XB1_1; raw_XB1_2; raw_XB1_3; raw_XB1_4; raw_XB1_5;};

save XB1_rawdata.mat
