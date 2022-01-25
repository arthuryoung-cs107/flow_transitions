classdef UB2_experiment < experiment
  properties
    statement = '45% solid fraction glycerine water UNDER buoyant suspension in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by under_buoy_102dens_45frac_rough1_105ml_oct22, under_buoy_102dens_45frac_rough2_105ml_oct22, under_buoy_102dens_45frac_rough3_105ml_oct22, under_buoy_102dens_45frac_rough4_105ml_oct22, under_buoy_102dens_45frac_rough5_105ml_oct22';
  end
  methods
    function obj = UB2_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'UB2';
      obj.color = color_;
      obj.specs = specs_;
      obj.TV_range = 36:45;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(10, :);

      mu_f = 0.0013273;
      rho_f = 1020.0;
      phi_m = 0.613;
      phi = 0.45;
      rho_p = 1040.0;
      for i=1:length(exp_list_in)
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
      obj.dat_num = obj.exp(2).dat_num; %% in this case, we know run 1 is broken

      %% consider omitting the first run entirely, is very different from the rest of runs. Reads MUCH lower.
      obj.exp(1).mu_rpm = [obj.exp(1).mu_rpm(1:31); NaN; obj.exp(1).mu_rpm(32:length(obj.exp(1).mu_rpm))];
      obj.exp(1).mu_torque = [obj.exp(1).mu_torque(1:31); NaN; obj.exp(1).mu_torque(32:length(obj.exp(1).mu_torque))];
      obj.exp(1).G = [obj.exp(1).G(1:31); NaN; obj.exp(1).G(32:length(obj.exp(1).G))];
      obj.exp(1).G_rat = [obj.exp(1).G_rat(1:31); NaN; obj.exp(1).G_rat(32:length(obj.exp(1).G_rat))];
      obj.exp(1).cf = [obj.exp(1).cf(1:31); NaN; obj.exp(1).cf(32:length(obj.exp(1).cf))];
      obj.exp(1).omega = [obj.exp(1).omega(1:31); NaN; obj.exp(1).omega(32:length(obj.exp(1).omega))];
      obj.exp(1).Re_s = [obj.exp(1).Re_s(1:31); NaN; obj.exp(1).Re_s(32:length(obj.exp(1).Re_s))];
      % obj.exp(3).mu_torque = [obj.exp(3).mu_torque(1:29); NaN; obj.exp(3).mu_torque(30:length(obj.exp(3).mu_torque))];
      % obj.exp(3).G = [obj.exp(3).G(1:29); NaN; obj.exp(3).G(30:length(obj.exp(3).G))];
      % obj.exp(3).cf = [obj.exp(3).cf(1:29); NaN; obj.exp(3).cf(30:length(obj.exp(3).cf))];
      % obj.exp(3).omega = [obj.exp(3).omega(1:29); NaN; obj.exp(3).omega(30:length(obj.exp(3).omega))];
      % obj.exp(3).Re_s = [obj.exp(3).Re_s(1:29); NaN; obj.exp(3).Re_s(30:length(obj.exp(3).Re_s))];

      T_mat = zeros(obj.dat_num, obj.len);
      G_mat = zeros(obj.dat_num, obj.len);
      G_rat_mat = zeros(obj.dat_num, obj.len);
      cf_mat = zeros(obj.dat_num, obj.len);
      omega_mat = zeros(obj.dat_num, obj.len);
      Re_mat = zeros(obj.dat_num, obj.len);

      T_mat(:,1) = obj.exp(1).mu_torque;
      G_mat(:,1) = obj.exp(1).G;
      G_rat_mat(:,1) = obj.exp(1).G_rat;
      cf_mat(:,1) = obj.exp(1).cf;
      omega_mat(:,1) = obj.exp(1).omega;
      Re_mat(:,1) = obj.exp(1).Re_s;
      % T_mat(:,3) = obj.exp(3).mu_torque;
      % G_mat(:,3) = obj.exp(3).G;
      % cf_mat(:,3) = obj.exp(3).cf;
      % omega_mat(:,3) = obj.exp(3).omega;
      % Re_mat(:,3) = obj.exp(3).Re_s;

      % for i=[2 4 5]
      for i=2:5
        T_mat(:, i) = obj.exp(i).mu_torque;
        G_mat(:, i) = obj.exp(i).G;
        G_rat_mat(:, i) = obj.exp(i).G_rat;
        cf_mat(:, i) = obj.exp(i).cf;
        omega_mat(:, i) = obj.exp(i).omega;
        Re_mat(:, i) = obj.exp(i).Re_s;
      end
      obj.mu_torque = mean(T_mat', 'omitnan');
      obj.sigma_torque = std(T_mat','omitnan');
      obj.G = mean(G_mat','omitnan');
      obj.G_rat = mean(G_rat_mat','omitnan');
      obj.cf = mean(cf_mat','omitnan');
      obj.omega = mean(omega_mat','omitnan');
      obj.Re_s = mean(Re_mat','omitnan');
    end
  end
end
