clear
path = '/Users/Arthur/Documents/anton_parr/';
filename = 'Arthur_1_13_2020 4_24 PM_NBrig_20frac_laminar_turbulent_run.xlsx';

raw_NB1_1 = rheometer_excel_parser2018([path filename], 1);
raw_NB1_2 = rheometer_excel_parser2018([path filename], 2);
raw_NB1_3 = rheometer_excel_parser2018([path filename], 3);
raw_NB1_4 = rheometer_excel_parser2018([path filename], 4);
raw_NB1_5 = rheometer_excel_parser2018([path filename], 5);

raw_NB1_structure = {raw_NB1_1; raw_NB1_2; raw_NB1_3; raw_NB1_4; raw_NB1_5;};

save NB1_rawdata.mat
