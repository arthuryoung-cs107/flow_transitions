clear
path = '/Users/Arthur/Documents/anton_parr/';
filename = 'Arthur_1_16_2020 4_30 PM_PureLiq_NBdens_ProfileCyl.xlsx';

raw_PF4_1 = rheometer_excel_parserCB([path 'Arthur_9_25_2019 6_17 PM_CB_NB_dens2'], 1);
raw_PF4_2 = rheometer_excel_parserCB([path 'Arthur_9_25_2019 6_17 PM_CB_NB_dens2'], 2);
raw_PF4_3 = rheometer_excel_parserCB([path 'Arthur_9_25_2019 6_17 PM_CB_NB_dens2_2'], 1);
raw_PF4_4 = rheometer_excel_parserCB([path 'Arthur_9_25_2019 6_17 PM_CB_NB_dens2_3'], 1);
raw_PF4_5 = rheometer_excel_parserCB([path 'Arthur_9_25_2019 6_17 PM_CB_NB_dens2_3'], 2);

raw_PF4_structure = {raw_PF4_1; raw_PF4_2; raw_PF4_3; raw_PF4_4; raw_PF4_5;};

save PF4_rawdata.mat
