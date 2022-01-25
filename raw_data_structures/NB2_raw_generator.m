clear
path = './anton_parr_data/';
filename = 'Arthur_1_15_2020 3_46 PM_NB_rigor_30frac_PREMIX2.xlsx';

raw_NB2_1 = rheometer_excel_parser2018([path filename], 1);
raw_NB2_2 = rheometer_excel_parser2018([path filename], 2);
raw_NB2_3 = rheometer_excel_parser2018([path filename], 3);
raw_NB2_4 = rheometer_excel_parser2018([path filename], 4);
raw_NB2_5 = rheometer_excel_parser2018([path filename], 5);
raw_NB2_6 = rheometer_excel_parser2018([path filename], 6);

raw_NB2_structure = {raw_NB2_1; raw_NB2_2; raw_NB2_3; raw_NB2_4; raw_NB2_5; raw_NB2_6;};

save NB2_rawdata.mat
