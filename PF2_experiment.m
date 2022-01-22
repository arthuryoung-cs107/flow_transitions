classdef PF2_experiment < experiment
  properties
    statement = 'Glycerin water composition in short aspect ratio geometry. Actuated by smooth inner cylinder, formerly generated and processed by glycol_0179_water_0821_smooth1_835ml_sep31, glycol_0179_water_0821_smooth2_835ml_sep31, glycol_0179_water_0821_smooth3_835ml_sep31, glycol_0179_water_0821_smooth4_835ml_sep31, glycol_0179_water_0821_smooth5_835ml_sep31';
  end
  methods
    function obj = PF2_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'PF2';
      obj.color = color_;
      obj.specs = specs_;
      obj.TV_range = 6:45; 

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
    end
  end
end
