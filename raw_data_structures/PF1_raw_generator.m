clear
path = '/Users/Arthur/Documents/anton_parr/';
filename = 'Arthur_10_8_2019 5_55 PM_60w_40g_long_laminar_turbulent_run';

raw_PF1_1 = rheometer_excel_parser2018([path filename], 1);
raw_PF1_2 = rheometer_excel_parser2018([path filename], 2);
raw_PF1_3 = rheometer_excel_parser2018([path filename], 3);
raw_PF1_4 = rheometer_excel_parser2018([path filename], 4);

raw_PF1_structure = {raw_PF1_1; raw_PF1_2; raw_PF1_3; raw_PF1_4;};

save PF1_rawdata.mat
