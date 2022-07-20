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
      reals = ~isnan(omega_all);
      omega_all = omega_all(reals);
      mu_torque_all = mu_torque_all(reals);
      sigma_torque_all = sigma_torque_all(reals);
      G_all = G_all(reals);
      G_rat_all = G_rat_all(reals);
      cf_all = cf_all(reals);
      Re_s_all = Re_s_all(reals);

      [obj.omega, obj.ind_sort] = mink(omega_all, length(omega_all));

      obj.mu_torque = mu_torque_all(obj.ind_sort);
      obj.sigma_torque = sigma_torque_all(obj.ind_sort);
      obj.G = G_all(obj.ind_sort);
      obj.G_rat = G_rat_all(obj.ind_sort);
      obj.cf = cf_all(obj.ind_sort);
      obj.Re_s = Re_s_all(obj.ind_sort);
    end
    function gen_powerfit(obj,Re_s_TV_, G_TV_, alpha_TV_)
        r_i = 0.01208;
        r_o = 0.025;

        if (nargin>1) % use FB_experiment fitting algorithm
            Re_s_TV = rmmissing(Re_s_TV_);
            G_TV = rmmissing(G_TV_);
            alpha_TV = rmmissing(alpha_TV_);
            [obj.Re_s_TV I_R_TV] = sort(Re_s_TV);
            obj.G_TV = G_TV(I_R_TV);
            obj.alpha_TV = alpha_TV(I_R_TV);

            obj.powerfit = fit(obj.Re_s_TV, obj.G_TV,'b*x^m', 'StartPoint', [70, 1]);
            alpha = obj.powerfit.m;
            beta = obj.powerfit.b;
            m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
            obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
            obj.G_c2 = beta*(obj.Re_sc2)^alpha;
        else % use experminet fitting algorithm
            alpha_tol = obj.alpha_tol;

            alpha_vec = obj.alpha;
            full_indices = 1:length(obj.omega);
            I_transitioned = logical((obj.Re_s>50) .* (alpha_vec>alpha_tol));
            I_transition = min(full_indices(I_transitioned));
            obj.TV_range = I_transition:length(obj.omega);
            Re_s_TV = obj.Re_s(obj.TV_range);
            obj.Re_sc1 = Re_s_TV(1);
            obj.Re_s_TV = Re_s_TV;
            obj.G_TV = obj.G(obj.TV_range);
            obj.alpha_TV = alpha_vec(obj.TV_range);

            obj.G_c1 = obj.G_TV(1);
            [obj.Re_sc3, obj.G_c3] = fluid.interp_trans(alpha_tol, obj.Re_s, alpha_vec, obj.G, I_transition);

            obj.powerfit = fit(reshape(Re_s_TV, [], 1), reshape(obj.G_TV,[],1),'b*x^m', 'StartPoint', [70, 1]);
            alpha = obj.powerfit.m;
            beta = obj.powerfit.b;
            m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
            obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
            obj.G_c2 = beta*(obj.Re_sc2)^alpha;
        end
    end
  end
end
