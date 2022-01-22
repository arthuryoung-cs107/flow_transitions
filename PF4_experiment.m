classdef PF4_experiment < experiment
  properties
    statement = 'Glycerin water composition in CUP AND BOB geometry. Varied temperature from 19 to 23 degrees Celsius. Actuated by SMOOTH inner cylinder, formerly generated and processed by CB_glycol_0179_water_0821_19Celsius_19ml_sep25, CB_glycol_0179_water_0821_20Celsius_19ml_sep25, CB_glycol_0179_water_0821_21Celsius_19ml_sep25, CB_glycol_0179_water_0821_22Celsius_19ml_sep25, CB_glycol_0179_water_0821_23Celsius_19ml_sep25';
    temp = [19, 20, 21, 22, 23];
  end
  methods
    function obj = PF4_experiment(exp_list_in, color_, specs_)
      obj@experiment(exp_list_in);
      obj.label = 'PF4';
      obj.color = color_;
      obj.specs = specs_;
      obj.exp(1).mu_f = 0.0021843;
      obj.exp(2).mu_f = 0.002148;
      obj.exp(3).mu_f = 0.0021351;
      obj.exp(4).mu_f = 0.0021119;
      obj.exp(5).mu_f = 0.0021013;

      rho_f = 1042.657;
      phi_m = 0.613;
      phi = 0;
      rho_p = 0;
      for i=1:length(exp_list_in)
        obj.exp(i).r_i = 0.01333;
        obj.exp(i).r_o = 0.01446;
        obj.exp(i).h = 0.040;
        obj.exp(i).rho_f = rho_f;
        obj.exp(i).phi_m = phi_m;
        obj.exp(i).phi = phi;
        obj.exp(i).rho_p = rho_p;
        obj.exp(i).mu_eff = obj.Krieger_Dougherty(obj.exp(i).mu_f, phi, phi_m);
      end
    end
  end
end
