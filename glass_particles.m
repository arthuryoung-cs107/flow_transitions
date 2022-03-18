classdef glass_particles < fluid
  properties
    tag = 'FB';
    Q;
    q;
    q_inc;
    tau_static;

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

    alpha;
    Re_s_alpha;

    dimless_ind; 
    Re_s_ns;
    cf_ns;
    G_ns;
    alpha_ns;

    TV_range;
    TV_lowRes = 50;
    % TV_lowomega = 10;
    TV_lowomega = 6;
    TV_lowalpha = 1.0;
    powerfit;
  end
  methods
    function obj = glass_particles(name_, color_)
      obj@fluid(name_, color_);
      obj.phi_m = 0.5869;
      obj.rho_p = 2500;
      obj.rho_f = 1.225;
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
    function sort_dimensionless(obj)
      obj.Re_s_ns = obj.Re_s;
      obj.cf_ns = obj.cf;
      obj.G_ns = obj.G;
      obj.alpha_ns = obj.alpha;

      [obj.Re_s, obj.dimless_ind] = mink(obj.Re_s, length(obj.Re_s));
      obj.cf = obj.cf_ns(obj.dimless_ind);
      obj.G = obj.G_ns(obj.dimless_ind);
    end
    function compute_tau_y(obj)
      r_i = obj.r_i;
      r_o = obj.r_o;
      h = obj.h;

      obj.tau_min = min(obj.tau);

      fit_range = 1:10;
      linfit = fit(obj.omega(fit_range), obj.tau(fit_range), 'poly1');

      if ((obj.tau_min > linfit.p2) && (linfit.p2 > 0))
        obj.tau_y = linfit.p2;
      else
        obj.tau_y = 0.99*obj.tau_min; %% safety factor to avoid singularities in evaluation
      end
      obj.tau_c = ((r_i/r_o)^(-2))*(obj.tau_y);
    end
    function compute_tau_y_old(obj)
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
    function compute_dimensionless(obj)
      nu = obj.mu_p/obj.rho_b;
      obj.gamma = zeros(size(obj.omega));
      obj.S = zeros(size(obj.omega));
      obj.Re_s = zeros(size(obj.omega));
      obj.cf = zeros(size(obj.mu_torque));

      for i = 1:length(obj.omega)
          if obj.tau(i) > obj.tau_c % no plug layer
              r_c = obj.r_o;
          else % plug layer
              r_c = (obj.r_i)*sqrt(obj.tau(i)/(obj.tau_y));
          end
          d = r_c - obj.r_i;
          obj.gamma(i) = abs((1/obj.mu_p)*(2*(obj.mu_p*obj.omega(i) + obj.tau_y*log(r_c/obj.r_i))/(obj.r_i^(-2) - r_c^(-2))*(obj.r_i)^(-2) - obj.tau_y)); %% opposite of what we should find taking outward radial direction as positive, as a result of absolute value

          r_t = sqrt(obj.r_i * r_c);
          obj.S(i) = abs(1/obj.mu_p*(2*(obj.mu_p*obj.omega(i) + obj.tau_y*log(r_c/obj.r_i))/(obj.r_i^(-2) - r_c^(-2))*(r_t)^(-2) - obj.tau_y)); %%
          obj.Re_s(i) = obj.S(i)*d^(2)/nu;
          obj.cf(i) = obj.mu_torque(i)/(2*pi*obj.r_i*obj.r_i*obj.h*obj.S(i)*obj.S(i)*d*d);
      end
    end
    function alpha_out = alpha_comp(obj)
      alpha_out = obj.alpha_G();
    end
    function alpha_out = alpha_cf(obj)
      alpha_out = approx_deriv_weighted_central(log(obj.Re_s), log(obj.cf))+2;
    end
    function alpha_out = alpha_G(obj)
      alpha_out = approx_deriv_weighted_central(log(obj.Re_s), log(obj.G));
    end
    function alpha_out = alpha_T(obj)
      alpha_out = approx_deriv_weighted_central(log(obj.omega), log(obj.mu_torque));
    end
    function Res_out = Re_s_alpha_comp(obj)
      Res_out = obj.Re_s;
    end
    function string_out = label(obj)
      if round( obj.q, 1) >= 10.0
        string_out = [obj.tag, ' q=', num2str(round(obj.q, 1),2)];
      else
        string_out = [obj.tag, ' q=', num2str(round(obj.q, 1),'%.1f')];
      end
    end
    function gen_powerfit(obj)
      % obj.TV_range = obj.Re_s > obj.TV_lowRes;
      obj.TV_range = obj.omega > obj.TV_lowomega;
      % obj.TV_range = (obj.omega>obj.TV_lowomega).*(obj.alpha_T>obj.TV_lowalpha);

      obj.powerfit = fit(obj.Re_s(obj.TV_range), obj.G(obj.TV_range),'b*x^m', 'StartPoint', [70, 1]);
    end
  end
end

function prime_vec = approx_deriv_weighted_central(t_in, x_in)
  n = length(t_in);
  prime_vec = nan(size(x_in));

  k = 5; %% number of points considered. Must be odd
  l = (k-1)/2; %% number of points to left and right
  p = l+1; %% index of central point

  x = reshape(x_in(1:k), [1 k]);
  t = reshape(t_in(1:k), [1 k]);
  for i=1:l
    ind = [1:(i-1), i+1:k];
    t_hat = t(ind);
    x_hat = x(ind);
    w = (abs((0.5*(t_hat + t(i)))-t(i))).^(-1);
    m = ((x(i)-x_hat))./(t(i)-t_hat);
    prime_vec(i) = (sum(w.*m))/(sum(w));
  end
  for i=p:(n-l)
    x = reshape(x_in(i-l:i+l), [1 k]);
    t = reshape(t_in(i-l:i+l), [1 k]);
    ind = [1:p-1, p+1:k];
    t_hat = t(ind);
    x_hat = x(ind);
    w = (abs((0.5*(t_hat + t(p)))-t(p))).^(-1);
    m = ((x(p)-x_hat))./(t(p)-t_hat);
    prime_vec(i) = (sum(w.*m))/(sum(w));
  end
  x = reshape(x_in(n-k+1:n), [1 k]);
  t = reshape(t_in(n-k+1:n), [1 k]);
  for i=p+1:k
    ind = [1:i-1, i+1:k];
    t_hat = t(ind);
    x_hat = x(ind);
    w = (abs((0.5*(t_hat + t(i)))-t(i))).^(-1);
    m = ((x(i)-x_hat))./(t(i)-t_hat);
    prime_vec(n-k+i) = (sum(w.*m))/(sum(w));
  end
end
