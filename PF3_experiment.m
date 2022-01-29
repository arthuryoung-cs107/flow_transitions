classdef PF3_experiment < experiment
  properties
    statement = 'Glycerin water composition in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by glycol_017_water_083_rough1_90ml_Jan16_2020, glycol_017_water_083_rough2_90ml_Jan16_2020, glycol_017_water_083_rough3_90ml_Jan16_2020, glycol_017_water_083_rough4_90ml_Jan16_2020, glycol_017_water_083_rough5_90ml_Jan16_2020';

    mu_torque_full;
    sigma_torque_full;
    G_full;
    G_rat_full;
    cf_full;
    omega_full;
    Re_s_full;
  end
  methods
    function obj = PF3_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'PF3';
      obj.color = color_;
      obj.specs = specs_;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(5, :);

      mu_f = 0.0020580;
      rho_f = 1042.657;
      phi_m = 0.613;
      phi = 0;
      rho_p = 0;
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

      obj.mu_torque(5) = NaN;
      obj.sigma_torque(5) = NaN;
      obj.G(5) = NaN;
      obj.G_rat(5) = NaN;
      obj.cf(5) = NaN;
      obj.omega(5) = NaN;
      obj.Re_s(5) = NaN;

      obj.mu_torque = rmmissing(obj.mu_torque);
      obj.sigma_torque = rmmissing(obj.sigma_torque);
      obj.G = rmmissing(obj.G);
      obj.G_rat = rmmissing(obj.G_rat);
      obj.cf = rmmissing(obj.cf);
      obj.omega = rmmissing(obj.omega);
      obj.Re_s = rmmissing(obj.Re_s);

      obj.len = length(obj.mu_torque);
    end
  end
end
