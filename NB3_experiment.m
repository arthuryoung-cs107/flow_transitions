classdef NB3_experiment < experiment
  properties
    statement = '40% solid fraction glycerin water neutrally buoyant suspension in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by neutral_buoy_40frac_rough1_90ml_Jan22_2020, neutral_buoy_40frac_rough2_90ml_Jan22_2020, neutral_buoy_40frac_rough3_90ml_Jan22_2020, neutral_buoy_40frac_rough4_90ml_Jan22_2020, neutral_buoy_40frac_rough5_90ml_Jan22_2020';

    Re_s_noKD = [111.647497432277e+000, 122.811809733215e+000, 133.976539724231e+000, 145.140747602649e+000, 156.305647437523e+000, 167.470062902877e+000, 178.634327395914e+000, 189.799894025190e+000, 200.963296717914e+000, 212.129341426196e+000, 223.293626811468e+000, 334.937938350805e+000, 446.586993732099e+000, 558.228856797488e+000, 669.877359834381e+000, 781.522256706819e+000, 893.167213613106e+000, 1.00481413502779e+003, 1.11646586973615e+003, 1.39557756525362e+003, 1.67469813092188e+003, 1.95382668097343e+003, 2.23294121741318e+003, 2.51206273825547e+003, 2.79116809247875e+003, 3.07035516046892e+003, 3.34933294308384e+003, 3.62850566870384e+003, 3.90778759763355e+003, 4.18690723132186e+003, 4.46629848937153e+003, 5.58250351260949e+003, 6.69898909273756e+003, 7.81579228087873e+003, 8.93233780959791e+003, 10.0488667942673e+003];
    G_rat_noKD = [7.46105178781639e+000, 7.28555719497451e+000, 7.33114551138475e+000, 7.62240298445360e+000, 8.28590705868544e+000, 8.13950264371826e+000, 8.26292799715098e+000, 8.18410256984848e+000, 8.34190760812408e+000, 8.65713287907896e+000, 8.58055177394443e+000, 7.93132480051330e+000, 8.32975803101617e+000, 8.75269183295471e+000, 9.20235657673445e+000, 9.97572566330073e+000, 11.1126224739206e+000, 12.0363585609550e+000, 13.0322445315184e+000, 15.0535091253391e+000, 16.7497762476125e+000, 18.5019050160858e+000, 20.1802973835870e+000, 21.8428420720845e+000, 23.4757520830574e+000, 25.1289679568125e+000, 26.7577225584881e+000, 28.2923000390065e+000, 29.8723703313086e+000, 31.6239775705323e+000, 33.3111403534540e+000, 38.8118302697449e+000, 43.0501670595636e+000, 48.5000215248871e+000, 51.9690279411992e+000, 55.2538061091623e+000];

    mu_torque_full;
    sigma_torque_full;
    G_full;
    G_rat_full;
    cf_full;
    omega_full;
    Re_s_full;
  end
  methods
    function obj = NB3_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'NB3';
      obj.color = color_;
      obj.specs = specs_;
      obj.TV_range = 17:37;

      mu_f = 0.0020580;
      rho_f = 1042.657;
      phi_m = 0.613;
      phi = 0.4;
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
    end
    function process_raws(obj, raws)
      for i=1:obj.len
        obj.exp(i).process_raw(raws{i});
      end
      obj.dat_num = obj.exp(1).dat_num;

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

      % truncate the first point
      obj.mu_torque_full = obj.mu_torque;
      obj.sigma_torque_full = obj.sigma_torque;
      obj.G_full = obj.G;
      obj.G_rat_full = obj.G_rat;
      obj.cf_full = obj.cf;
      obj.omega_full = obj.omega;
      obj.Re_s_full = obj.Re_s;

      plot_range = 2:obj.dat_num;
      obj.mu_torque = obj.mu_torque(plot_range);
      obj.sigma_torque = obj.sigma_torque(plot_range);
      obj.G = obj.G(plot_range);
      obj.G_rat = obj.G_rat(plot_range);
      obj.cf = obj.cf(plot_range);
      obj.omega = obj.omega(plot_range);
      obj.Re_s = obj.Re_s(plot_range);

      obj.dat_num = obj.dat_num - 1;
    end
  end
end
