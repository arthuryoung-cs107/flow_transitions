classdef PF3_experiment < experiment
  properties
    statement = 'Glycerin water composition in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by glycol_017_water_083_rough1_90ml_Jan16_2020, glycol_017_water_083_rough2_90ml_Jan16_2020, glycol_017_water_083_rough3_90ml_Jan16_2020, glycol_017_water_083_rough4_90ml_Jan16_2020, glycol_017_water_083_rough5_90ml_Jan16_2020';
  end
  methods
    function obj = PF3_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'PF3';
      obj.color = color_;
      obj.specs = specs_;    

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
