classdef PF1_experiment < experiment
  properties
    statement = 'Glycerin water composition in short aspect ratio geometry. Actuated by smooth inner cylinder, formerly generated and processed by glycol_040_water_060_smooth1_90ml_oct8, glycol_040_water_060_smooth2_90ml_oct8, glycol_040_water_060_smooth3_90ml_oct8, glycol_040_water_060_smooth4_90ml_oct8';
  end
  methods
    function obj = PF1_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'PF1';
      obj.color = color_;
      obj.specs = specs_;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(3, :);

      mu_f = 0.0057866;
      rho_f = 1102.76;
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
  end
end
