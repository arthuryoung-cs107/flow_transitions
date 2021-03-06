classdef UB1_experiment < experiment
  properties
    statement = '45% solid fraction WATER UNDER buoyant suspension in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by under_buoy_0994dens_45frac_rough1_105ml_oct21, under_buoy_0994dens_45frac_rough2_105ml_oct21, under_buoy_0994dens_45frac_rough3_105ml_oct21, under_buoy_0994dens_45frac_rough4_105ml_oct21, under_buoy_0994dens_45frac_rough5_105ml_oct21';
  end
  methods
    function obj = UB1_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'UB1';
      obj.color = color_;
      obj.specs = specs_;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(9, :);
      obj.TV_range = 37:45;
      obj.TV_lowRes = 360; 

      mu_f = 0.00098069;
      rho_f = 994.0;
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
      obj.mu_f = mu_f;
      obj.rho_f = rho_f;
      obj.phi_m = phi_m;
      obj.phi = phi;
      obj.rho_p = rho_p;
      obj.mu_eff = obj.Krieger_Dougherty(mu_f, phi, phi_m);

    end
  end
end
