classdef glass49 < glass_particles
  properties
    %% specific to Abhi's experiment
    RD_flow_lmin;
    RU_flow_lmin;
    RD_torque;
    RU_torque;
    RD_rpm;
    RU_rpm;
    omega_min = 1e-2;

    omega_full;
    tau_full;

    fix_tauy_flag=true;
    % fix_tauy_flag=false;

    glass_id = 2;
  end
  methods
    function obj = glass49(name_, color_)
      obj@glass_particles(name_, color_);
      obj.tag = '\textit{FB2}';
      obj.phi = 0.50;
      obj.q_inc = 0.2;
      obj.tau_static = 159.582557507215e+000;

      % obj.FB_Bingham_tauy_bounds = glass_particles.FB2_Bingham_tauy_bounds_prethin;
      % obj.FB_Bingham_mup_bounds = glass_particles.FB2_Bingham_mup_bounds_prethin;
      % obj.FB_Bingham_omega_cap = glass_particles.FB2_Bingham_omega_cap_prethin;
      % obj.FB_Bingham_wscheme = glass_particles.FB2_Bingham_wscheme_prethin;
      obj.FB_Bingham_tauy_bounds = glass_particles.FB2_Bingham_tauy_bounds_incthin;
      obj.FB_Bingham_mup_bounds = glass_particles.FB2_Bingham_mup_bounds_incthin;
      obj.FB_Bingham_omega_cap = glass_particles.FB2_Bingham_omega_cap_incthin;
      obj.FB_Bingham_wscheme = glass_particles.FB2_Bingham_wscheme_incthin;

      % obj.FB_fitted_Bingham_pars = glass_particles.FB2_fitted_Bingham_pars;
      % obj.mu_p_Bingham_intr = glass_particles.FB2_mu_p_qmax_alphac3;

      % obj.qcrit_sgf = 0.8;
      % obj.qcrit_sgf = 0.85;
      obj.qcrit_sgf = 0.9;
      obj.qcrit_ttv = 1.0;
      obj.qcrit_fgm = 2.0;

      obj.ocrit_dgf_ql = 80;
      obj.ocrit_dgf_qh = 20;
      obj.ocrit_fgm_ql = 15;
      obj.ocrit_fgm_qh = 15;

      % obj.qcrit_Bingham2Carreau = 1;
      obj.qcrit_Bingham2Carreau = obj.qcrit_fgm;
      obj.FB_phi_low = glass_particles.FB2_phi_low;
    end
    function process_raw(obj,raw_,i_)
      obj.exp_id = i_;
      obj.rho_b = obj.rho_p * obj.phi + obj.rho_f*(1.0-obj.phi);

      obj.RD_flow_lmin = raw_(:, 2);
      obj.RU_flow_lmin = raw_(:, 12);
      obj.mu_flow_lmin = 0.5*(obj.RD_flow_lmin + flip(obj.RU_flow_lmin));

      obj.Q = mean([obj.RD_flow_lmin; obj.RU_flow_lmin]);

      obj.q = obj.Q/obj.q_inc;

      obj.RD_torque = raw_(:, 3)/1000^(2); % in Newton meters
      obj.RU_torque = raw_(:, 13)/1000^(2);
      obj.mu_torque = 0.5*(obj.RU_torque + flip(obj.RD_torque));
      obj.sigma_torque = (std([obj.RD_torque, flip(obj.RU_torque)]'))'/1000^(2);

      obj.RD_rpm = raw_(:, 4);
      obj.RU_rpm = raw_(:, 14);
      obj.mu_rpm = 0.5*(obj.RU_rpm + flip(obj.RD_rpm));

      obj.omega = 2*pi/60*obj.mu_rpm;
      obj.tau = obj.mu_torque/(2*pi*(obj.r_i^2)*obj.h);

      ind_use =  obj.omega > obj.omega_min;

      obj.omega_full = obj.omega;
      obj.tau_full = obj.tau;

      obj.mu_torque = obj.mu_torque(ind_use);
      obj.sigma_torque = obj.sigma_torque(ind_use);
      obj.mu_rpm = obj.mu_rpm(ind_use);
      obj.omega = obj.omega(ind_use);
      obj.tau = obj.tau(ind_use);

      %% temporary hack
      if ((i_==11 || i_==12) && obj.fix_tauy_flag)
          % obj.compute_tau_y;
          obj.tau_min=min(obj.tau);
          if (i_==11)
              obj.tau_y=21e-3;
          elseif (i_==12)
              obj.tau_y=19e-3;
          end
          obj.tau_c= (obj.r_o*obj.r_o)*obj.tau_y/(obj.r_i*obj.r_i);
      else
          obj.compute_tau_y;
      end

      obj.compute_appmu;
      obj.compute_mu_plastic;

      obj.compute_dimensionless;
      obj.G = obj.mu_torque/((obj.h)*(obj.mu_p*obj.mu_p)/(obj.rho_b));

      obj.omega_crit = glass_particles.determine_omega_crit(obj.alpha_tol_Bingham, obj.omega_crit_min,obj.omega,obj.alpha_T);
      obj.omega_crit_alpha_peak = glass_particles.determine_omega_crit_alpha_peak(obj.omega,obj.alpha_T);

      % obj.omega_cap_Bingham_use = obj.omega_crit;
      % obj.omega_cap_Bingham_use = obj.omega_crit_alpha_peak;
      obj.omega_cap_Bingham_use = obj.FB_Bingham_omega_cap(obj.exp_id);
      % omega_fit = obj.omega;
      % tau_fit = obj.tau;

      obj.fit_Carreau_model;
      obj.fit_Bingham_model(obj.omega_cap_Bingham_use,obj.omega,obj.tau);
      obj.Carreau_fit_muinf_params = obj.fit_Carreau_model_muinf;
      % obj.mu_p_Bingham = obj.FB_fitted_Bingham_pars(i_,1);
      % obj.tau_y_Bingham = obj.FB_fitted_Bingham_pars(i_,2);

      obj.compute_dimensionless_Bingham_new;
    end
  end
end
