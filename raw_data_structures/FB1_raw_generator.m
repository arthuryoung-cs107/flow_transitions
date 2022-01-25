clear
path = './anton_parr_data/';
name = 'Arthur_9_8_2018 4_16 PM_1116g_130M_Q1-3_incfluid_1_2';

raw_FB1_Q100 = rheometer_excel_parser2018([path, name], 1);
raw_FB1_Q120 = rheometer_excel_parser2018([path, name], 2);
raw_FB1_Q140 = rheometer_excel_parser2018([path, name], 3);
raw_FB1_Q160 = rheometer_excel_parser2018([path, name], 4);
raw_FB1_Q180 = rheometer_excel_parser2018([path, name], 5);
raw_FB1_Q200 = rheometer_excel_parser2018([path, name], 6);
raw_FB1_Q220 = rheometer_excel_parser2018([path, name], 7);
raw_FB1_Q240 = rheometer_excel_parser2018([path, name], 8);
raw_FB1_Q260 = rheometer_excel_parser2018([path, name], 9);
raw_FB1_Q280 = rheometer_excel_parser2018([path, name], 10);

raw_FB1_structure = {raw_FB1_Q100; raw_FB1_Q120; raw_FB1_Q140; raw_FB1_Q160; raw_FB1_Q180; raw_FB1_Q200; raw_FB1_Q220; raw_FB1_Q240; raw_FB1_Q260; raw_FB1_Q280};

save FB1_rawdata.mat
