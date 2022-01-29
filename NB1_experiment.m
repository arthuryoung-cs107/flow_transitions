classdef NB1_experiment < experiment
  properties
    statement = '20% solid fraction glycerin water neutrally buoyant suspension in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by neutral_buoy_20frac_rough1_107ml_Jan13_2020, neutral_buoy_20frac_rough2_107ml_Jan13_2020, neutral_buoy_20frac_rough3_107ml_Jan13_2020, neutral_buoy_20frac_rough4_107ml_Jan13_2020, neutral_buoy_20frac_rough5_107ml_Jan13_2020';
    Re_s_noKD = [11.1678821268074e+000, 22.3323058206277e+000, 33.4972302854920e+000, 44.6620645702111e+000, 55.8277916949124e+000, 66.9934054431742e+000, 78.1583727533730e+000, 89.3232919604387e+000, 100.489009289047e+000, 111.652572487340e+000, 122.818876576580e+000, 133.984403827268e+000, 145.149073977321e+000, 156.314621652777e+000, 167.479861585527e+000, 178.644302887205e+000, 189.810613856199e+000, 200.975113078558e+000, 212.141794974954e+000, 223.298287283587e+000, 334.954069669110e+000, 446.596345111655e+000, 558.245224188175e+000, 669.894125583202e+000, 781.542258632607e+000, 893.192390849013e+000, 1.00484226537050e+003, 1.11648811451876e+003, 1.39561731293474e+003, 1.67474004139492e+003, 1.95388061714343e+003, 2.23300662981298e+003, 2.51213647646476e+003, 2.79113845083146e+003, 3.07041304837697e+003, 3.34960284164194e+003, 3.62878955772379e+003, 3.90779364483919e+003, 4.18700180557012e+003, 4.46614918377786e+003, 5.58274213331372e+003, 6.69928670812742e+003, 7.81590218551632e+003, 8.93247398276357e+003, 10.0489717970832e+003]';
    G_rat_noKD = [43.4121349696673e+000, 31.2050645869660e+000, 24.7226879509009e+000, 21.2315923771158e+000, 18.7520384241210e+000, 17.3170970630237e+000, 16.3362666947061e+000, 15.6788392126087e+000, 14.9998542014054e+000, 14.5869547060254e+000, 14.0795033878704e+000, 13.7065337808678e+000, 13.4861025836951e+000, 13.3296927812730e+000, 13.1952285914906e+000, 13.0320993188093e+000, 12.7181805892966e+000, 12.3655473390117e+000, 11.9089005587169e+000, 11.4228147327654e+000, 7.57594772774348e+000, 4.85578088965472e+000, 5.05367211254479e+000, 5.61722323533307e+000, 6.13614662581531e+000, 6.73955990206216e+000, 7.34963672440027e+000, 7.93273682219142e+000, 9.18608147526410e+000, 10.5453281856474e+000, 11.9395544625060e+000, 13.2101222163982e+000, 14.4413625437022e+000, 15.6567753239092e+000, 16.8006490636002e+000, 17.9848888577046e+000, 19.0881137944178e+000, 20.2568040903291e+000, 21.3718831347910e+000, 22.4916813470732e+000, 27.0378241485276e+000, 29.8071277579937e+000, 33.3241193625762e+000, 36.5124564540953e+000, 39.2844452049903e+000]';

    mu_torque_full;
    sigma_torque_full;
    G_full;
    G_rat_full;
    cf_full;
    omega_full;
    Re_s_full;
  end
  methods
    function obj = NB1_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'NB1';
      obj.color = color_;
      obj.specs = specs_;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(6, :);
      obj.TV_range = 22:45;
      obj.TV_lowRes = 219;

      mu_f = 0.0020580;
      rho_f = 1042.657;
      phi_m = 0.613;
      phi = 0.2;
      rho_p = 1042.657;
      for i=1:length(exp_list_in)
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
      obj.dat_num = obj.exp(2).dat_num; %% in this case, we know run 1 is broken

      T_mat = zeros(obj.dat_num, obj.len);
      G_mat = zeros(obj.dat_num, obj.len);
      G_rat_mat = zeros(obj.dat_num, obj.len);
      cf_mat = zeros(obj.dat_num, obj.len);
      omega_mat = zeros(obj.dat_num, obj.len);
      Re_mat = zeros(obj.dat_num, obj.len);

      T_mat(:, 1) = obj.exp(1).mu_torque(8:obj.exp(1).dat_num);
      G_mat(:, 1) = obj.exp(1).G(8:obj.exp(1).dat_num);
      G_rat_mat(:, 1) = obj.exp(1).G_rat(8:obj.exp(1).dat_num);
      cf_mat(:, 1) = obj.exp(1).cf(8:obj.exp(1).dat_num);
      omega_mat(:, 1) = obj.exp(1).omega(8:obj.exp(1).dat_num);
      Re_mat(:, 1) = obj.exp(1).Re_s(8:obj.exp(1).dat_num);

      for i=2:obj.len
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

      % truncate the dimensionless data so that we have nicer plots
      obj.mu_torque_full = obj.mu_torque;
      obj.sigma_torque_full = obj.sigma_torque;
      obj.G_full = obj.G;
      obj.G_rat_full = obj.G_rat;
      obj.cf_full = obj.cf;
      obj.omega_full = obj.omega;
      obj.Re_s_full = obj.Re_s;

      % plot_range = 22:45;
      plot_range = 1:45;
      obj.mu_torque = obj.mu_torque(plot_range);
      obj.sigma_torque = obj.sigma_torque(plot_range);
      obj.G = obj.G(plot_range);
      obj.G_rat = obj.G_rat(plot_range);
      obj.cf = obj.cf(plot_range);
      obj.omega = obj.omega(plot_range);
      obj.Re_s = obj.Re_s(plot_range);

    end
  end
end
