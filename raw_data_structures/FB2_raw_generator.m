clear
path = '/Users/Arthur/Documents/anton_parr/';
name = 'Abhi_49micron_data_sorted2';

raw_FB2_1 = rheometer_excel_parser_Abhi([path name], 1);
raw_FB2_2 = rheometer_excel_parser_Abhi([path name], 2);
raw_FB2_3 = rheometer_excel_parser_Abhi([path name], 3);
raw_FB2_4 = rheometer_excel_parser_Abhi([path name], 4);
raw_FB2_5 = rheometer_excel_parser_Abhi([path name], 5);
raw_FB2_6 = rheometer_excel_parser_Abhi([path name], 6);
raw_FB2_7 = rheometer_excel_parser_Abhi([path name], 7);
raw_FB2_8 = rheometer_excel_parser_Abhi([path name], 8);
raw_FB2_9 = rheometer_excel_parser_Abhi([path name], 9);
raw_FB2_10 = rheometer_excel_parser_Abhi([path name], 10);
raw_FB2_11 = rheometer_excel_parser_Abhi([path name], 11);
raw_FB2_12 = rheometer_excel_parser_Abhi([path name], 12);
raw_FB2_13 = rheometer_excel_parser_Abhi([path name], 13);

raw_FB2_structure = {raw_FB2_1; raw_FB2_2; raw_FB2_3; raw_FB2_4; raw_FB2_5; raw_FB2_6; raw_FB2_7; raw_FB2_8; raw_FB2_9; raw_FB2_10; raw_FB2_11; raw_FB2_12; raw_FB2_13;};

save FB2_rawdata.mat
