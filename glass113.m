classdef glass113 < glass_particles
  properties
    Tf_clean=0;
    nf_clean=0;
    qf_clean=0;
    if_clean=0;

    T_stats=0;
    n_stats=0;
    q_stats=0;

    clean_flag=true;
  end
  methods
    function obj = glass113(name_, color_, phi_)
      obj@glass_particles(name_, color_);
      obj.tag = 'FB1';
      obj.phi = phi_;
      obj.q_inc = 1.4;
      obj.tau_static = 225.046219998818e+000;
    end

    function [Tf_clean nf_clean qf_clean if_clean] = get_clean_data(obj)
        [Tf_raw nf_raw qf_raw if_raw] = concat_sort_raw(obj.true_raw);
        [Tf_clean nf_clean qf_clean if_clean] = clean_raw_bin(Tf_raw,nf_raw,qf_raw,if_raw);
    end
    function [T_out n_out q_out] = get_clean_stats(obj)
        if (length(obj.Tf_clean) == 1)
            [Tf_clean nf_clean qf_clean if_clean]=obj.get_clean_data;
        else
            [Tf_clean nf_clean qf_clean if_clean]=deal(obj.Tf_clean,obj.nf_clean,obj.qf_clean,obj.if_clean);
        end
        [T_out n_out q_out] = comp_bin_stats(Tf_clean,nf_clean,qf_clean,if_clean);
    end
    function [mu0_out, lambda_out, n_out, k_out, muinf_out] = fit_Carreau_fluid(obj, omega_cap_, full_flag_)
        omega_full = 2*pi/60*obj.nf_clean;
        tau_full = obj.Tf_clean/(2*pi*(obj.r_i^2)*obj.h);
        if (nargin==1)
            omega_cap=15;
            full_flag=true;
        elseif (nargin==2)
            omega_cap=15;
            full_flag=true;
        elseif (nargin==3)
            omega_cap=15;
            full_flag=full_flag_;
        end
        ind = omega_full<omega_cap;
        omega_fit=omega_full(ind);
        tau_fit=tau_full(ind);
        % w_fit = glass_particles.compute_distance_weighting(omega_fit);
        w_fit = ones(size(omega_fit))/norm(ones(size(omega_fit)));
        [mu0_out, lambda_out, n_out, k_out, muinf_out] = obj.fit_internal_Carreau_fluid(omega_fit,tau_fit,w_fit,full_flag);
    end
    function process_raw(obj, raw,i_)
      obj.true_raw=raw;
      obj.rho_b = obj.rho_p * obj.phi + obj.rho_f*(1.0-obj.phi);

      if (obj.clean_flag)
          [obj.Tf_clean obj.nf_clean obj.qf_clean obj.if_clean] = obj.get_clean_data;
          [obj.T_stats obj.n_stats obj.q_stats] = obj.get_clean_stats;

          obj.mu_torque=obj.T_stats(:,1);
          obj.sigma_torque=obj.T_stats(:,2);
          obj.mu_flow_lmin=obj.q_stats(:,1);
          obj.mu_rpm=obj.n_stats(:,1);
      else
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
      end
      obj.omega = 2*pi/60*obj.mu_rpm;
      obj.tau = obj.mu_torque/(2*pi*(obj.r_i^2)*obj.h);
      obj.Q = mean(obj.mu_flow_lmin);
      obj.q = obj.Q/obj.q_inc;

      obj.compute_tau_y;
      obj.compute_appmu;
      obj.compute_mu_plastic;

      obj.compute_dimensionless;
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

function [T_out n_out q_out] = comp_bin_stats(Tf_,nf_,qf_,if_)
    [T_out,n_out,q_out] = deal(nan(max(if_),2));
    for i = 1:max(if_)
        % T_out(i,:) = [mean(Tf_(if_==i)) std(Tf_(if_==i))];
        T_out(i,:) = [median(Tf_(if_==i)) std(Tf_(if_==i))];
        n_out(i,:) = [mean(nf_(if_==i)) std(nf_(if_==i))];
        q_out(i,:) = [mean(qf_(if_==i)) std(qf_(if_==i))];
    end
end
function [Tf_clean nf_clean qf_clean if_clean] = clean_raw_bin(Tf_,nf_,qf_, if_, w_)
    if (nargin==4)
        get_clean = @(t_) ~(isoutlier(t_));
    else
        get_clean = @(t_) ~(isoutlier(t_, 'movmedian', w_));
    end
    clean_ind_full=logical(zeros(size(if_)));
    for i = 1:max(if_)
        fcheck = if_==i;
        clean_ind_full(fcheck)=get_clean(Tf_(fcheck));
    end
    Tf_clean = Tf_(clean_ind_full);
    nf_clean = nf_(clean_ind_full);
    qf_clean = qf_(clean_ind_full);
    if_clean = if_(clean_ind_full);
end
function [Tf_out nf_out qf_out if_out] = concat_sort_raw(raw_)
    [raw_short,irange] = true_range_finder(raw_);
    nsets = size(irange,2);
    if_raw=nan(size(raw_short,1),1);
    for i = 1:nsets
        if_raw(irange(1,i):irange(2,i))=i;
    end
    [nf_out,inds]=sort(raw_short(:,6));
    Tf_out=raw_short(inds,3)*(1e-6);
    qf_out=raw_short(inds,2);
    if_out=if_raw(inds);
end
function [raw_short iout] = true_range_finder(raw_in)
    raw_=raw_in;
    raw_(raw_(:,6)<0, 6) = NaN;
    raw_(raw_(:,3)<0, 3) = NaN;
    raw_short=rmmissing(raw_);
    inds_full = 1:size(raw_short,1);
    diff_raw = diff(raw_short(:,1));
    dec_diff = [false;diff_raw<1];
    set_starts = inds_full(dec_diff);
    nsets = sum(dec_diff) + 1;
    iout = nan(2, nsets);
    iout(1,2:end) = set_starts;
    iout(2,1:(nsets-1)) = set_starts-1;
    iout(1,1)=1;
    iout(2,end)=size(raw_short,1);
end
