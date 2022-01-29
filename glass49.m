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
  end
  methods
    function obj = glass49(name_, color_)
      obj@glass_particles(name_, color_);
      obj.tag = 'FB2';
      obj.phi = 0.50;
      obj.q_inc = 0.2;
      obj.tau_static = 159.582557507215e+000;
      obj.TV_lowRes = 50;
    end
    function process_raw(obj, raw_)
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

      obj.mu_torque = obj.mu_torque(ind_use);
      obj.sigma_torque = obj.sigma_torque(ind_use);
      obj.mu_rpm = obj.mu_rpm(ind_use);
      obj.omega = obj.omega(ind_use);
      obj.tau = obj.tau(ind_use);

      obj.compute_tau_y;
      obj.compute_appmu;
      obj.compute_mu_plastic;

      obj.compute_gamma_S_bingham;
      obj.compute_Re_s;
      obj.G = obj.mu_torque/((obj.h)*(obj.mu_p*obj.mu_p)/(obj.rho_b));
    end
  end
end
