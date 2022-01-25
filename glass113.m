classdef glass113 < fluid
  properties
    q_inc = 1.4;
    tau_static = 225.046219998818e+000;

    Q;
    q;

    mu_flow_lmin;

    tau_y;
    tau_c;
    tau_min;
    mu_p;
    omega_p;
    appmu_min_sim;
    omega_appmu_min_sim;
    appmu_min;
    omega_appmu_min;

    sim_omega;
    sim_appmu;
  end
  methods
    function obj = glass113(name_, color_, phi_)
      obj@fluid(name_, color_);
      obj.phi = phi_;
      obj.phi_m = 0.5869;
      obj.rho_p = 2500;
      obj.rho_f = 1.225;

    end
    function process_raw(obj, raw)
      obj.rho_b = obj.rho_p * obj.phi + obj.rho_f*(1.0-obj.phi);

      obj.range_array_full = obj.range_finder_old(raw);
      obj.range_array = obj.steady_state(obj.range_array_full, raw);
      raw_ = obj.spike_filter(obj.range_array, raw);

      obj.mu_torque = zeros(length(obj.range_array), 1);
      obj.sigma_torque = zeros(length(obj.range_array), 1);
      obj.mu_flow_lmin = zeros(length(obj.range_array), 1);
      obj.mu_rpm = zeros(length(obj.range_array), 1);
      obj.time_durations = zeros(length(obj.range_array), 1);
      obj.meas_points = zeros(length(obj.range_array), 1);

      for i = 1:length(obj.range_array)
          obj.mu_torque(i) = mean(raw_(obj.range_array(1, i):obj.range_array(2, i),3)/1000^(2), 'omitnan');
          obj.sigma_torque(i) = std(raw_(obj.range_array(1, i):obj.range_array(2, i),3)/1000^(2), 'omitnan');
          obj.mu_flow_lmin(i) = mean(raw_(obj.range_array(1, i):obj.range_array(2, i),2), 'omitnan');
          obj.mu_rpm(i) = mean(raw_(obj.range_array(1, i):obj.range_array(2, i), 6), 'omitnan');
          obj.time_durations(i) = raw_(obj.range_array_full(2, i), 8) - raw_(obj.range_array_full(1, i), 8);
          obj.meas_points(i) = raw_(obj.range_array_full(2, i), 1) - raw_(obj.range_array_full(1, i), 1);
      end
      obj.omega = 2*pi/60*obj.mu_rpm;
      obj.tau = obj.mu_torque/(2*pi*(obj.r_i^2)*obj.h);
      obj.Q = mean(obj.mu_flow_lmin);
      obj.q = obj.Q/obj.q_inc;

      obj.compute_tau_y;
      obj.compute_appmu;
      obj.compute_mu_plastic;

      obj.compute_gamma_S_bingham;
      obj.compute_Re_s;
      obj.G = obj.mu_torque/((obj.h)*(obj.mu_p*obj.mu_p)/(obj.rho_b));
    end

    function steady_array = steady_state(obj, range_array, raw)
      steady_array = range_array;
      steady_value = zeros(1, length(range_array));
      steady_dev = zeros(1, length(range_array));
      for i = 1:length(range_array)
          steady_array(1, i) = range_array(2, i)-round((range_array(2, i)-range_array(1, i))*(0.90));
          steady_value(i) = mean(raw(steady_array(1, i):steady_array(2, i), 3), 'omitnan');
          steady_dev(i) = std(raw(steady_array(1, i):steady_array(2, i), 3), 'omitnan');
      end
      for i = 1:length(steady_array)
          for j =  range_array(1, i):range_array(2, i)
              if abs(raw(j, 3) - steady_value(i)) < 1*steady_dev(i)
                  steady_array(1, i) = j;
                  break
              else
                  continue
              end
          end
      end
    end
    function new_raw = spike_filter(obj, steady_array, raw )
      new_raw = raw;
      steady_value = zeros(1, length(steady_array));
      steady_dev = zeros(1, length(steady_array));
      threshold = 2.0;
      for i = 1:length(steady_array)
          for j = steady_array(1, i):steady_array(2, i)
              if raw(j, 3) < 0
                  new_raw(j, 3) = NaN;
              end
          end
      end
      for i = 1:length(steady_array)
          steady_value(i) = mean(raw(steady_array(1, i):steady_array(2, i), 3), 'omitnan');
          steady_dev(i) = std(raw(steady_array(1, i):steady_array(2, i), 3), 'omitnan');
      end
      for i = 1:length(steady_array)
          for j = steady_array(1, i):steady_array(2, i)
              if abs(raw(j, 3)-steady_value(i)) > threshold*steady_dev(i)
                  new_raw(j, 3) = NaN;
              end
          end
      end

    end
    function fig_out = plot_torques(obj, position)
      run figure_properties.m
      fig_out = figure('Position', fig_pos(position, :));
      set(gca, 'YScale', 'log')
      set(gca, 'XScale', 'log')
      ylabel('Torque [N.m]')
      xlabel('rotational speed [rpm]')
      hold on
      plot(obj.mu_rpm, obj.mu_torque, ' -', 'Color', green4, 'LineWidth', 1.0, 'DisplayName', 'mean')
      legend('Show', 'Location', 'SouthEast')
    end
    function compute_tau_y(obj)
      r_i = obj.r_i;
      r_o = obj.r_o;
      h = obj.h;

      % fit_range = 1:3;
      fit_range = 1:10;

      obj.tau_y = min(obj.tau);
      obj.tau_min = obj.tau_y;
      linfit = fit(obj.omega(fit_range), obj.tau(fit_range), 'poly1');
      if ((obj.tau_y > linfit.p2) && (linfit.p2 > 0))
        obj.tau_y = linfit.p2;
      end
      obj.tau_c = ((r_i/r_o)^(-2))*(obj.tau_y);
    end

    function compute_appmu(obj)
      r_i = obj.r_i;
      r_o = obj.r_o;
      h = obj.h;

      obj.appmu = zeros(size(obj.tau));
      for i = 1:length(obj.tau)
          if obj.tau(i) > obj.tau_c % no plug layer
              r_c = r_o;
          else % plug layer
              r_c = r_i*sqrt(obj.tau(i)/obj.tau_y);
          end
          obj.appmu(i) = (obj.tau_y*log(obj.tau_y / obj.tau(i)) + 0.5*obj.tau(i)*(1 + (r_i/r_c)^2 ) )/obj.omega(i);
      end
    end
    function compute_mu_plastic(obj)
      obj.appmu_min = min(obj.appmu);
      obj.omega_appmu_min = obj.omega(find(obj.appmu == obj.appmu_min));
      obj.mu_p = obj.appmu_min;
      obj.omega_p = obj.omega_appmu_min;

      fit_range = (find(obj.appmu == obj.appmu_min) - 3):(min([length(obj.appmu), find(obj.appmu == obj.appmu_min) + 3]));
      obj.sim_omega = (obj.omega(min(fit_range))):0.5:(obj.omega(max(fit_range)));

      appmu_omega_bingham_fit = fit( obj.omega(fit_range), obj.appmu(fit_range), 'poly3');
      obj.sim_appmu = appmu_omega_bingham_fit.p1*(obj.sim_omega).^(3) + appmu_omega_bingham_fit.p2*(obj.sim_omega).^(2) + appmu_omega_bingham_fit.p3*(obj.sim_omega).^(1) + appmu_omega_bingham_fit.p4;

      obj.appmu_min_sim = min(obj.sim_appmu);
      obj.omega_appmu_min_sim = obj.sim_omega(find(obj.sim_appmu == obj.appmu_min_sim));

      if (obj.appmu_min_sim < obj.appmu_min)
        obj.mu_p = obj.appmu_min_sim;
        obj.omega_p = obj.omega_appmu_min_sim;
      end
    end
    function compute_gamma_S_bingham(obj)
      obj.gamma = zeros(size(obj.omega));
      obj.S = zeros(size(obj.omega));
      for i = 1:length(obj.gamma)
          if obj.tau(i) > obj.tau_c % no plug layer
              r_c = obj.r_o;
          else % plug layer
              r_c = (obj.r_i)*sqrt(obj.tau(i)/(obj.tau_y));
          end
          obj.gamma(i) = (1/obj.mu_p)*(2*(obj.mu_p*obj.omega(i) + obj.tau_y*log(r_c/obj.r_i))/(obj.r_i^(-2) - r_c^(-2))*(obj.r_i)^(-2) - obj.tau_y); %% opposite of what we should find taking outward radial direction as positive, as a result of absolute value

          r_t = sqrt(obj.r_i * r_c);
          obj.S(i) = 1/obj.mu_p*(2*(obj.mu_p*obj.omega(i) + obj.tau_y*log(r_c/obj.r_i))/(obj.r_i^(-2) - r_c^(-2))*(r_t)^(-2) - obj.tau_y); %%
      end
    end
    function compute_Re_s(obj)
      nu = obj.mu_p/obj.rho_b;
      obj.Re_s = zeros(size(obj.omega));
      for i = 1:length(obj.S)
        if obj.tau(i) > obj.tau_c
          d = obj.r_o - obj.r_i;
        else
          d = (obj.r_i)*sqrt(obj.tau(i)/(obj.tau_y)) - obj.r_i;
        end
        obj.Re_s(i) = obj.S(i)*d^(2)/nu;
      end
    end
    function string_out = label(obj)
      string_out = ['FB1 q=', num2str(round( obj.q, 1))];
    end
  end
end
