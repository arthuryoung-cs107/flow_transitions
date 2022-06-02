r_i = 0.01208;
r_o = 0.025;
L = 0.036;
eta = r_i/r_o; % gap ratio

G_obs_Res_slope = (2*pi*r_i*r_o)/((r_o-r_i)^2);

omega_tau_range = [1e-2 1e2 1e-9 1e-2];
Res_G_range = [1e-2 1e4 1e0 1e8];
Res_alpha_range = [1e0 5e5 -1 2.1];
Res_cf_range = [1e-1 1e4 1e-4 1e4];
Res_Grat_range = [1e-1 2e4 0.8 1e2];
omega_appmu_range = [1e-2 2e2 6e-2 1e7];

omega_range = [1e-2, 2e2];
torque_range = [1e-7, 1e-2];
Re_s_range = [1e-4, 2e4];
G_range = [1e-1, 1e5];
Grat_range = [5e-1, 3e1];
cf_range = [1e1, 1e10];
alpha_range = [0, 3];

pos11 = [0 500];
pos12 = [450 500];
pos13 = [900 500];
pos21 = [0 100];
pos22 = [450 100];
pos23 = [900 100];

pos_spread = [0 500; 180 500; 360 500; 540 500; 720 500; 900 500; 0 100; 180 100; 360 100; 540 100; 720 100; 900 100];

textbox_pos2_a_NW = [0.09, 0.85, 0.1, 0.1];
textbox_pos2_b_NW = [0.575, 0.85, 0.1, 0.1];
textbox_pos2_a_SW = [0.09, 0.175, 0.1, 0.1];
textbox_pos2_b_SW = [0.575, 0.175, 0.1, 0.1];

textbox_pos21_a_NW = [0.12, 0.85, 0.1, 0.1];
textbox_pos21_b_NW = [0.55, 0.85, 0.1, 0.1];
textbox_pos21_c_NW = [0.12, 0.07, 0.1, 0.1];
textbox_pos21_a_SW = [0.12, 0.675, 0.1, 0.1];
textbox_pos21_b_SW = [0.55, 0.675, 0.1, 0.1];
textbox_pos21_c_SW = [0.12, 0.1, 0.1, 0.1];

textbox_pos222_a = [0.4, 0.625, 0.1, 0.1];
textbox_pos222_b = [0.9, 0.625, 0.1, 0.1];
textbox_pos222_c = [0.4, 0.32, 0.1, 0.1];
textbox_pos222_d = [0.9, 0.32, 0.1, 0.1];
textbox_pos222_e = [0.4, 0.008, 0.1, 0.1];
textbox_pos222_f = [0.9, 0.008, 0.1, 0.1];

textbox_pos22_a_NE = [0.425, 0.85, 0.1, 0.1];
textbox_pos22_b_NE = [0.9, 0.85, 0.1, 0.1];
textbox_pos22_c_NE = [0.425, 0.375, 0.1, 0.1];
textbox_pos22_d_NE = [0.9, 0.375, 0.1, 0.1];

txtbx_pos_mat = [   textbox_pos2_a_NW; ...
                    textbox_pos2_b_NW; ...
                    textbox_pos2_a_SW; ...
                    textbox_pos2_b_SW; ...
                    textbox_pos21_a_NW; ...
                    textbox_pos21_b_NW; ...
                    textbox_pos21_c_NW; ...
                    textbox_pos21_a_SW; ...
                    textbox_pos21_b_SW; ...
                    textbox_pos21_c_SW; ...
                    textbox_pos222_a; ...
                    textbox_pos222_b; ...
                    textbox_pos222_c; ...
                    textbox_pos222_d; ...
                    textbox_pos222_e; ...
                    textbox_pos222_f; ...

                    ];



dim2_short = [580 250];
dim2_tall = [580 330];

dim21_short = [580 330];
dim21_tall = [580 450];

dim222_short = [580 600];
dim222_tall = [580 800];

dim22_short = [580 300];
dim22_tall = [580 500];

dim1 = [580 325];

% posdimfull = [1 1 1438 796];
% posdimfull = [1 1 1728 1000];
posdimfull = [0 0 1728 1000];

ax21_short = [0.1 0.1 0.85 0.475; 0.1 0.675 0.4 0.3; 0.545 0.675 0.4 0.3];
ax21_tall = [0.1 0.1 0.85 0.475; 0.1 0.675 0.4 0.3; 0.545 0.675 0.4 0.3];

dim2 = dim2_short;
dim21 = dim21_short;
dim222 = dim222_tall;
dim22 = dim22_short;

ax21 = ax21_short;
