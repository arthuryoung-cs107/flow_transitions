classdef XB1_experiment < experiment
  properties
    statement = '45% solid fraction glycerine water OVER buoyant suspension in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by excess_buoy_1075dens_40frac_rough1_90ml_oct16, excess_buoy_1075dens_40frac_rough2_90ml_oct16, excess_buoy_1075dens_40frac_rough3_90ml_oct16, excess_buoy_1075dens_40frac_rough4_90ml_oct16, excess_buoy_1075dens_40frac_rough5_90ml_oct16';
  end
  methods
    function obj = XB1_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'XB1';
      obj.color = color_;
      obj.specs = specs_;
      obj.TV_range = 33:45;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(11, :);

      mu_f = 0.0036981;
      rho_f = 1075.0;
      phi_m = 0.613;
      phi = 0.40;
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
  end
end
