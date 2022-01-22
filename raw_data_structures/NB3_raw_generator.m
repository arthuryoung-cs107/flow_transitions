clear
path = '/Users/Arthur/Documents/anton_parr/';
filename = 'Arthur_1_22_2020 5_34 PM_NB_rigor_40frac_PREMIX2.xlsx';

raw_NB3_1 = rheometer_excel_parser2018([path filename], 1);
raw_NB3_2 = rheometer_excel_parser2018([path filename], 2);
raw_NB3_3 = rheometer_excel_parser2018([path filename], 3);
raw_NB3_4 = rheometer_excel_parser2018([path filename], 4);
raw_NB3_5 = rheometer_excel_parser2018([path filename], 5);

raw_NB3_structure = {raw_NB3_1; raw_NB3_2; raw_NB3_3; raw_NB3_4; raw_NB3_5;};

save NB3_rawdata.mat
