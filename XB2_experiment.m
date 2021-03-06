classdef XB2_experiment < experiment
  properties
    statement = '45% solid fraction glycerine water OVER buoyant suspension in short aspect ratio geometry. Actuated by rough inner cylinder, formerly generated and processed by excess_buoy_109dens_40frac_rough1_90ml_oct17, excess_buoy_109dens_40frac_rough2_90ml_oct17, excess_buoy_109dens_40frac_rough3_90ml_oct17, excess_buoy_109dens_40frac_rough4_90ml_oct17, excess_buoy_109dens_40frac_rough5_90ml_oct17';
  end
  methods
    function obj = XB2_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'XB2';
      obj.color = color_;
      obj.specs = specs_;
      fig_pos = fig_pos_gen(2, 6);
      obj.def_pos = fig_pos(12, :);
      obj.TV_range = 37:45;      
      obj.TV_lowRes = 200;

      mu_f = 0.0037046;
      rho_f = 1090.0;
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
      obj.mu_f = mu_f;
      obj.rho_f = rho_f;
      obj.phi_m = phi_m;
      obj.phi = phi;
      obj.rho_p = rho_p;
      obj.mu_eff = obj.Krieger_Dougherty(mu_f, phi, phi_m);

    end
  end
end
