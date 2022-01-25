clear
path = './anton_parr_data/';
filename = 'Arthur_1_16_2020 4_30 PM_PureLiq_NBdens_ProfileCyl.xlsx';

raw_PF3_1 = rheometer_excel_parser2018([path filename], 1);
raw_PF3_2 = rheometer_excel_parser2018([path filename], 2);
raw_PF3_3 = rheometer_excel_parser2018([path filename], 3);
raw_PF3_4 = rheometer_excel_parser2018([path filename], 4);
raw_PF3_5 = rheometer_excel_parser2018([path filename], 5);

raw_PF3_structure = {raw_PF3_1; raw_PF3_2; raw_PF3_3; raw_PF3_4; raw_PF3_5;};

save PF3_rawdata.mat
