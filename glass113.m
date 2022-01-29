classdef glass113 < glass_particles
  properties
  end
  methods
    function obj = glass113(name_, color_, phi_)
      obj@glass_particles(name_, color_);
      obj.tag = 'FB1';
      obj.phi = phi_;
      obj.q_inc = 1.4;
      obj.tau_static = 225.046219998818e+000;
      obj.TV_lowRes = 50;
    end
    function process_raw(obj, raw)
      obj.rho_b = obj.rho_p * obj.phi + obj.rho_f*(1.0-obj.phi);

      obj.range_array_full = obj.range_finder(raw);
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
    function plot_raw_torques(obj, raw_, raw)
      run figure_properties.m
      fig_pos2 = fig_pos_gen(2, 2);
      fig_full = figure('Name', 'full measuring period', 'Renderer', 'painters', 'Position', fig_pos2(1, :));
      ylabel('Torque [N.m]')
      xlabel('time')
      xlim([0, 1.5]);
      set(gca, 'YScale', 'log')
      hold on

      fig_trimmed = figure('Name', 'trimmed measuring period', 'Renderer', 'painters', 'Position', fig_pos2(3, :));
      ylabel('Torque [N.m]')
      xlabel('time')
      xlim([0, 1.5]);
      set(gca, 'YScale', 'log')
      hold on

      plot_range = 1:obj.dat_num;

      figure(fig_full.Number)
      for i = plot_range
          plot(raw(obj.range_array_full(1, i):obj.range_array_full(2, i), 1)/obj.meas_points(i), raw(obj.range_array_full(1, i):obj.range_array_full(2, i), 3)/(1e6), ' -', 'Color', colors_big(i, :), 'LineWidth', 1.0, 'MarkerFaceColor', colors_big(i, :), 'MarkerEdgeColor', colors_big(i, :), 'DisplayName', [num2str(obj.mu_rpm(i))])
      end
      legend('Show', 'Location', 'NorthEast')

      figure(fig_trimmed.Number)
      for i = plot_range
        plot(raw_(obj.range_array(1, i):obj.range_array(2, i), 1)/obj.meas_points(i), raw_(obj.range_array(1, i):obj.range_array(2, i), 3)/(1e6), ' -', 'Color', colors_big(i, :), 'LineWidth', 1.0, 'MarkerFaceColor', colors_big(i, :), 'MarkerEdgeColor', colors_big(i, :), 'DisplayName', [num2str(obj.mu_rpm(i))])
      end
      legend('Show', 'Location', 'NorthEast')
    end
  end
end
