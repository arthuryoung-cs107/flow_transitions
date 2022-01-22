classdef super_experiment < experiment
  properties
    ind_sort;
  end
  methods
    function obj = super_experiment(exp_full_list_in_)
      obj@experiment(exp_full_list_in_);

      omega_all = obj.exp{1}.omega;
      mu_torque_all = obj.exp{1}.mu_torque;
      sigma_torque_all = obj.exp{1}.sigma_torque;
      G_all = obj.exp{1}.G;
      G_rat_all = obj.exp{1}.G_rat;
      cf_all = obj.exp{1}.cf;
      Re_s_all = obj.exp{1}.Re_s;

      for i=2:obj.len
        omega_all = [omega_all obj.exp{i}.omega];
        mu_torque_all = [mu_torque_all obj.exp{i}.mu_torque];
        sigma_torque_all = [sigma_torque_all obj.exp{i}.sigma_torque];
        G_all = [G_all obj.exp{i}.G];
        G_rat_all = [G_rat_all obj.exp{i}.G_rat];
        cf_all = [cf_all obj.exp{i}.cf];
        Re_s_all = [Re_s_all obj.exp{i}.Re_s];
      end

      [obj.omega, obj.ind_sort] = mink(omega_all, length(omega_all));

      obj.mu_torque = mu_torque_all(obj.ind_sort);
      obj.sigma_torque = sigma_torque_all(obj.ind_sort);
      obj.G = G_all(obj.ind_sort);
      obj.G_rat = G_rat_all(obj.ind_sort);
      obj.cf = cf_all(obj.ind_sort);
      obj.Re_s = Re_s_all(obj.ind_sort);
    end
    function process_raws(obj, raws)
      %% do nothing
    end

  end
end
