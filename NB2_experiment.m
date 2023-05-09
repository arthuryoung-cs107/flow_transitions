classdef NB2_experiment < experiment
  properties
    statement = '30% solid fraction glycerin water neutrally buoyant suspension in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by neutral_buoy_30frac_rough1_91ml_REHASH_PREMIX_Jan15_2020, neutral_buoy_30frac_rough2_91ml_REHASH_PREMIX_Jan15_2020, neutral_buoy_30frac_rough3_91ml_REHASH_PREMIX_Jan15_2020, neutral_buoy_30frac_rough4_91ml_REHASH_PREMIX_Jan15_2020, neutral_buoy_30frac_rough5_91ml_REHASH_PREMIX_Jan15_2020, neutral_buoy_30frac_rough6_91ml_REHASH_PREMIX_Jan15_2020';

    Re_s_noKD = [66.9887307401391e+000, 78.1530522567938e+000, 89.3174698881127e+000, 100.482429442538e+000, 111.647280201059e+000, 122.811864127998e+000, 133.976984505688e+000, 145.141184793319e+000, 156.305912687497e+000, 167.470803086601e+000, 178.634330541171e+000, 189.799671760389e+000, 200.963581887847e+000, 212.129295295766e+000, 223.292140324960e+000, 334.939149649498e+000, 446.584106529158e+000, 558.229141933649e+000, 669.876694690431e+000, 781.522642102799e+000, 893.164950033877e+000, 1.00481024721990e+003, 1.11645926385996e+003, 1.39556896483669e+003, 1.67468314875853e+003, 1.95381317937012e+003, 2.23292520001032e+003, 2.51205030459500e+003, 2.79112330797446e+003, 3.07031591923620e+003, 3.34952737276292e+003, 3.62863847665427e+003, 3.90772832011923e+003, 4.18690374603164e+003, 4.46601131697309e+003, 5.58263591183497e+003, 6.69916354116474e+003, 7.81577758529717e+003, 8.93228915508627e+003, 10.0489088977678e+003];
    G_rat_noKD = [5.57772369437179e+000, 6.27028424481601e+000, 6.21510791295519e+000, 6.16628452479742e+000, 6.00234510205130e+000, 5.44059016437777e+000, 5.31894356848767e+000, 5.33614430116489e+000, 5.19637143674892e+000, 5.37223746753718e+000, 5.20927081389096e+000, 5.17117137322011e+000, 5.31695931388234e+000, 5.10804621128099e+000, 4.56754428041829e+000, 5.51805325623943e+000, 5.67353031463335e+000, 6.12754525660147e+000, 7.41024625911531e+000, 8.07053574361605e+000, 8.72145727482458e+000, 9.31028624823766e+000, 10.2670429616709e+000, 12.0190762000193e+000, 13.6875256622675e+000, 15.2714166297824e+000, 17.0046040035521e+000, 18.5809750562921e+000, 20.1784336636625e+000, 21.6705655004507e+000, 23.0349816222358e+000, 23.7432539779732e+000, 24.6229676856304e+000, 25.7098412244316e+000, 27.0407120378487e+000, 32.1058575183084e+000, 36.6089599285651e+000, 40.9322642456737e+000, 44.2422739377172e+000, 47.1817618557063e+000];

    mu_torque_full;
    sigma_torque_full;
    G_full;
    G_rat_full;
    cf_full;
    omega_full;
    Re_s_full;
  end
  methods
    function obj = NB2_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = '\textit{NB2}';
      obj.color = color_;
      obj.specs = specs_;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(7, :);

      mu_f = 0.0020580;
      rho_f = 1042.657;
      phi_m = 0.613;
      phi = 0.3;
      rho_p = 1042.657;
      for i=1:length(exp_list_in)
        % obj.exp(i).range_finder_flag = 2;
        obj.exp(i).steady_state_flag = 3;
        obj.exp(i).mu_f = mu_f;
        obj.exp(i).rho_f = rho_f;
        obj.exp(i).phi_m = phi_m;
        obj.exp(i).phi = phi;
        obj.exp(i).rho_p = rho_p;
        obj.exp(i).mu_eff = obj.Krieger_Dougherty(mu_f, phi, phi_m);
      end
      obj.mu_f = mu_f;
      obj.rho_f = rho_f;
      obj.phi_m = phi_m;
      obj.phi = phi;
      obj.rho_p = rho_p;
      obj.mu_eff = obj.Krieger_Dougherty(mu_f, phi, phi_m);

    end
    function process_raws(obj, raws)
      for i=1:obj.len
        obj.exp(i).process_raw(raws{i});
      end
      obj.dat_num = obj.exp(1).dat_num; %% unless otherwise specified
      T_mat = zeros(obj.dat_num, obj.len);
      G_mat = zeros(obj.dat_num, obj.len);
      G_rat_mat = zeros(obj.dat_num, obj.len);
      cf_mat = zeros(obj.dat_num, obj.len);
      omega_mat = zeros(obj.dat_num, obj.len);
      Re_mat = zeros(obj.dat_num, obj.len);
      for i=1:obj.len
        T_mat(:, i) = obj.exp(i).mu_torque;
        G_mat(:, i) = obj.exp(i).G;
        G_rat_mat(:, i) = obj.exp(i).G_rat;
        cf_mat(:, i) = obj.exp(i).cf;
        omega_mat(:, i) = obj.exp(i).omega;
        Re_mat(:, i) = obj.exp(i).Re_s;
      end
      obj.mu_torque = mean(T_mat', 'omitnan');
      obj.sigma_torque = std(T_mat', 'omitnan');
      obj.G = mean(G_mat', 'omitnan');
      obj.G_rat = mean(G_rat_mat', 'omitnan');
      obj.cf = mean(cf_mat', 'omitnan');
      obj.omega = mean(omega_mat', 'omitnan');
      obj.Re_s = mean(Re_mat', 'omitnan');

      obj.mu_torque_full = obj.mu_torque;
      obj.sigma_torque_full = obj.sigma_torque;
      obj.G_full = obj.G;
      obj.G_rat_full = obj.G_rat;
      obj.cf_full = obj.cf;
      obj.omega_full = obj.omega;
      obj.Re_s_full = obj.Re_s;

      obj.mu_torque(15) = NaN;
      obj.sigma_torque(15) = NaN;
      obj.G(15) = NaN;
      obj.G_rat(15) = NaN;
      obj.cf(15) = NaN;
      obj.omega(15) = NaN;
      obj.Re_s(15) = NaN;

      obj.mu_torque = rmmissing(obj.mu_torque);
      obj.sigma_torque = rmmissing(obj.sigma_torque);
      obj.G = rmmissing(obj.G);
      obj.G_rat = rmmissing(obj.G_rat);
      obj.cf = rmmissing(obj.cf);
      obj.omega = rmmissing(obj.omega);
      obj.Re_s = rmmissing(obj.Re_s);

      obj.len = length(obj.mu_torque);
    end
    function gen_powerfit(obj)
        obj.gen_powerfit_internal(100);
    end
  end
end
