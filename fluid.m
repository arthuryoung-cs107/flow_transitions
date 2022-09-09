classdef fluid < handle
  properties
    range_finder_flag = 1;
    steady_state_flag = 1;
    spike_filter_flag = 0;
    noKD_flag = 0;
    r_i = 0.01208;
    r_o = 0.025;
    h = 0.036;

    phi;
    phi_m;
    rho_p;
    rho_f;
    mu_eff;
    mu_f;

    name;
    color;
    specs;

    %% specific to our experiment
    range_array_full;
    range_array;
    time_durations;
    meas_points;
    dat_num;

    mu_torque;
    sigma_torque;
    tau;
    appmu;
    G;
    G_rat;
    cf;

    mu_rpm;
    omega;
    gamma;
    S;
    Re_s;

    rho_b;

    true_raw;
  end
  methods
    function obj = fluid(name_, color_)
      obj.name = name_;
      obj.color = color_;
    end
    function process_raw(obj, raw)
      obj.true_raw = raw; 
      obj.rho_b = obj.rho_p * obj.phi + obj.rho_f*(1.0-obj.phi);

      h = obj.h;
      mu = obj.mu_eff;
      rho = obj.rho_b;
      ri = obj.r_i;
      ro = obj.r_o;
      d = ro-ri;
      nu = mu/rho;

      switch obj.range_finder_flag
        case 1
          obj.range_array_full = range_finder(raw);
        case 2
          obj.range_array_full = range_finder_old(raw);
        case 3
          obj.range_array_full = range_finder_weird(raw);
      end

      switch obj.steady_state_flag
        case 1
          obj.range_array = obj.steady_state_old(obj.range_array_full, raw);
        case 2
          obj.range_array = obj.steady_state(obj.range_array_full, raw);
        case 3
          obj.range_array = obj.steady_state_variable(obj.range_array_full, raw);
      end

      switch obj.spike_filter_flag
        case 0
          raw_ = raw; %% apply no spike filter
        case 1
          raw_ = obj.spike_filter(obj.range_array, raw);
        case 2
          raw_ = obj.spike_filter_old(obj.range_array, raw);
      end

      obj.dat_num = length(obj.range_array);

      obj.mu_torque = zeros(obj.dat_num, 1);
      obj.sigma_torque = zeros(obj.dat_num, 1);
      obj.mu_rpm = zeros(obj.dat_num, 1);
      obj.time_durations = zeros(obj.dat_num, 1);
      obj.meas_points = zeros(obj.dat_num, 1);

      obj.omega = zeros(obj.dat_num, 1);
      obj.gamma = zeros(obj.dat_num, 1);
      obj.Re_s = zeros(obj.dat_num, 1);

      obj.tau = zeros(obj.dat_num, 1);
      obj.G = zeros(obj.dat_num, 1);
      obj.G_rat = zeros(obj.dat_num, 1);
      obj.cf = zeros(obj.dat_num, 1);

      for i = 1:obj.dat_num
          obj.mu_torque(i) = mean(raw_(obj.range_array(1, i):obj.range_array(2, i),3)/1000^(2), 'omitnan');
          obj.sigma_torque(i) = std(raw_(obj.range_array(1, i):obj.range_array(2, i),3)/1000^(2), 'omitnan');
          obj.mu_rpm(i) = mean(raw_(obj.range_array(1, i):obj.range_array(2, i), 6), 'omitnan');
          obj.time_durations(i) = raw_(obj.range_array_full(2, i), 8) - raw_(obj.range_array_full(1, i), 8);
          obj.meas_points(i) = raw_(obj.range_array_full(2, i), 1) - raw_(obj.range_array_full(1, i), 1);

          rpm_it = obj.mu_rpm(i);
          T_it = obj.mu_torque(i);

          obj.omega(i) = 2*pi/60*rpm_it;
          obj.gamma(i) = 2*obj.omega(i)/(ri*ri*(1/(ri*ri) - 1/(ro*ro)));
          obj.Re_s(i) = obj.gamma(i)*(ri/ro)*d*d/nu;

          obj.tau(i) = T_it/(2*pi*(ri*ri)*h);
          obj.G(i) = T_it/(h*mu*mu/rho);
          obj.G_rat(i) = obj.G(i)/(2*pi*ri*ri*obj.gamma(i)*rho/mu);
          U_it = ri/ro*obj.gamma(i)*d;
          obj.cf(i) = T_it/(2*pi*rho*ri*ri*h*U_it*U_it);
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
    function range_array_out = range_finder(obj, raw_)
      range_array_out = range_finder(raw_);
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
    function steady_array = steady_state_old(obj, range_array, raw)
      steady_array = range_array; % initialised to same size and values as unprocessed range array
      steady_value = zeros(1, length(range_array)); % vector of steady state values
      steady_dev = zeros(1, length(range_array));
      for i = 1:length(range_array)
        steady_array(1, i) = range_array(2, i)-round((range_array(2, i)-range_array(1, i))/3); % takes last third of the data range for steady state value
        steady_value(i) = mean(raw(steady_array(1, i):steady_array(2, i), 3), 'omitnan');
        steady_dev(i) = std(raw(steady_array(1, i):steady_array(2, i), 3), 'omitnan');
      end
      for i = 1:length(steady_array)
        for j =  range_array(1, i):range_array(2, i)
          if abs(raw(j, 3) - steady_value(i)) < 2*steady_dev(i) % eschews data outside of one standard deviation of steady state value
            steady_array(1, i) = j;
            break
          else
            continue
          end
        end
      end
    end
    function steady_array = steady_state_variable(obj, range_array, raw)
      steady_array = range_array; % initialised to same size and values as unprocessed range array
  %     time = 150;
      time = 90;
  %     time = 600;
      for i = 1:length(range_array)
          if raw(range_array(2, i) , 8) - raw(range_array(1, i) , 8) >= (time - 1)
              interval_raw = raw( range_array(1, i):range_array(2, i) , :);
              end_time = interval_raw(1 , 8) + (time-1) ; %% should be 150, but playing it safe:
              j = 0;
              over_time = 0;
              while over_time == 0
                  j = j+1;
                  if interval_raw(j, 8) >= end_time
                      over_time = 1;
                  end
              end
              cutoff_index = range_array(1, i) + j;
              steady_array(2, i) = cutoff_index;
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
    function new_raw = spike_filter_old(obj, steady_array, raw )
      steady_section = steady_array;
      steady_value = zeros(1, length(steady_array)); % vector of steady state values
      steady_dev = zeros(1, length(steady_array));

      for i = 1:length(steady_array)
          steady_section(1, i) = steady_array(2, i)-round((steady_array(2, i)-steady_array(1, i))/3); % takes last fifth of the steady data range for steady state value
          raw(steady_section(1, i):steady_section(2, i), 3)
          steady_value(i) = mean(raw(steady_section(1, i):steady_section(2, i), 3), 'omitnan'); % finds mean of last fifth of steady data range, torque
          steady_dev(i) = std(raw(steady_section(1, i):steady_section(2, i), 3), 'omitnan'); % finds standard deviation of last fifth of steady data range, torque
      end

      for i = 1:length(steady_array)
          for j = steady_array(1, i):steady_array(2, i)
              if abs(raw(j, 3)-steady_value(i)) > 2.25*steady_dev(i)
                  raw(j, 3) = NaN;
              end
          end
      end
      new_raw = raw;
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
  end
  methods (Static)
      function [Rc, Gc] = interp_trans(alpha_tol_,R_,alpha_,G_,It_)
          R2=R_(It_);
          alpha2=alpha_(It_);
          G2=G_(It_);

          R1=R_(It_-1);
          alpha1=alpha_(It_-1);
          G1=G_(It_-1);

          Rc=(alpha_tol_-alpha2)*((R2-R1)/(alpha2-alpha1)) + R2;
          Gc=(Rc-R2)*((G2-G1)/(R2-R1)) + G2;
      end
  end
end

function length = interval_length_finder(raw_array,i)
    length_found = 0;
    index = 0;
    while length_found == 0
        index = index + 1;
        if isnan(raw_array(i - index, 2)) == 0
            length = raw_array(i - index, 2);
            length_found = 1;
        end
    end
end

function range_array = range_finder(raw_array)
  one_ind_raw = find(raw_array(:, 1)==1);
  n_raw = length(one_ind_raw);
  n = 0;
  for i=1:n_raw
    if raw_array(one_ind_raw(i) + 1, 1) == 2
      n = n + 1;
      one_ind(n) = one_ind_raw(i);
    end
  end
  range_array = zeros(2, n);
  range_array(1, :) = one_ind;
  for i=1:n-1
    index = one_ind(i+1)-1;
    search_on = true;
    while search_on
      if isnan(raw_array(index, 1))
        index = index - 1;
      else
        if raw_array(index-1, 1) == (raw_array(index, 1)-1)
          range_array(2, i) = index;
          search_on = false;
        else
          index = index - 1;
        end
      end
    end
  end
  range_array(2, n) = size(raw_array, 1);
end

function range_array = range_finder_old(raw_array )
  range_array = [];
  i = 0;
  interval_start = (1:5)';
  data_end = 0;
  current_range = [0; 0];
  long_check = size(raw_array);
  while data_end == 0
    i = i + 1;
    if raw_array(i:i+4, 1) == interval_start %% found the beginning of the interval
        length_add = interval_length_finder(raw_array, i);
        current_range = [i; i+length_add-1];
        range_array = [range_array, current_range];
        i = i+length_add-1; %% skips to the end of the interval
    end
    if current_range(2) == long_check(1)
        data_end = 1;
    end
  end
end

function range_array = range_finder_weird(raw_array )
  range_array = [7; 0];
  index_check = 0;
  for i = 7:length(raw_array)
    if index_check == 0
      if isnan(raw_array(i, 8)) == 1 % if beginning of gap or invalid value
        if (isnan(raw_array(i+1, 8)) == 1 && isnan(raw_array(i+2,8)) == 1 && (raw_array(i+7) == 1 || raw_array(i+4) == 1)) % if gap, begins skip routine
          index_check = 1;
          continue
        else % if invalid value, treats same as if it is not equal to nan
          dimensions = size(range_array);
          range_array(2, dimensions(2)) = i;
          continue
        end
      else % if this is not equal to nan, regular data point
        dimensions = size(range_array);
        range_array(2, dimensions(2)) = i;
        continue
      end
    else
      if isnan(raw_array(i+1, 8)) == 0 %if next value detected to be data point
        index_check = 0;
        new_column = [i + 1; 0];
        range_array = [range_array, new_column];
        continue
      else
        continue
      end
    end
  end
end
