classdef experiment < handle
  properties
    label;
    color;
    specs;
    % LW = 1.0;
    % MS = 5.0;
    LW = 0.7;
    MS = 3.5;
    LW_L = 1.0;
    MS_L = 8.0;
    def_pos;

    mu_f;
    rho_f;
    phi_m;
    phi;
    rho_p;

    mu_torque;
    sigma_torque;
    G;
    G_rat;
    cf;

    omega;
    Re_s;

    TV_range;

    powerfit;

    exp;
    len; %% this is the number of data sets
    dat_num; %% this is the length of the measurement sets
  end
  methods
    function obj = experiment(exp_list_in)
      obj.exp = exp_list_in;
      obj.len = length(exp_list_in);
    end
    function process_raws(obj, raws)
      for i=1:obj.len
        obj.exp(i).process_raw(raws{i});
      end
      obj.dat_num = obj.exp(1).dat_num; %% unless otherwise specified
      T_mat = zeros(obj.dat_num, obj.len);
      G_mat = zeros(obj.dat_num, obj.len);
      G_rat_mat = zeros(obj.dat_num, obj.len);
      cf_mat = zeros(obj.dat_num, obj.len);
      omega_mat = zeros(obj.dat_num, obj.len);
      Re_mat = zeros(obj.dat_num, obj.len);
      for i=1:obj.len
        T_mat(:, i) = obj.exp(i).mu_torque;
        G_mat(:, i) = obj.exp(i).G;
        G_rat_mat(:, i) = obj.exp(i).G_rat;
        cf_mat(:, i) = obj.exp(i).cf;
        omega_mat(:, i) = obj.exp(i).omega;
        Re_mat(:, i) = obj.exp(i).Re_s;
      end
      obj.mu_torque = mean(T_mat', 'omitnan');
      obj.sigma_torque = std(T_mat', 'omitnan');
      obj.G = mean(G_mat', 'omitnan');
      obj.G_rat = mean(G_rat_mat', 'omitnan');
      obj.cf = mean(cf_mat', 'omitnan');
      obj.omega = mean(omega_mat', 'omitnan');
      obj.Re_s = mean(Re_mat', 'omitnan');
    end
    function mu_eff = Krieger_Dougherty(obj, mu_f, phi, phi_m)
      mu_eff = mu_f*(1-phi/phi_m)^(-1.82);
    end
    function Latex_label = LTX_label(obj)
      Latex_label = ['$$' obj.label '$$'];
    end
    function alpha_out = alpha(obj)
      alpha_out = approx_deriv_weighted_central(log(obj.Re_s), log(obj.cf)) + 2;
    end
    function Res_out = Re_s_alpha(obj)
      Res_out = obj.Re_s;
    end
    function alpha_out = alpha_old(obj)
      alpha_out = approx_deriv_2ndO_legrangian(log(obj.Re_s), log(obj.cf)) + 2;
    end
    function Res_out = Re_s_alpha_old(obj)
      Res_out = obj.Re_s(3:(length(obj.Re_s)-2));
    end
    function gen_powerfit(obj)
      obj.powerfit = fit(obj.Re_s(obj.TV_range)', obj.G(obj.TV_range)','b*x^m', 'StartPoint', [70, 1]);
    end
    function inspect_torques(obj)
      fig_specs = AYfig.specs_gen(obj.label, obj.def_pos);
      fig_out = AYfig.figure(fig_specs);
      hold on
      set(gca, 'YScale', 'log')
      set(gca, 'XScale', 'log')
      for i=1:length(obj.exp)
        plot(obj.exp(i).omega, obj.exp(i).mu_torque, '- o','Color', obj.exp(i).color, 'LineWidth', 1.5, 'DisplayName', obj.exp(i).name)
      end
      ylabel('$$T$$ [N.m]', 'Interpreter', 'LaTeX','FontSize',12)
      xlabel('$$\omega_i$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
      legend('Show')
      hold off
    end
  end
end

function prime_vec = approx_deriv_2ndO_legrangian(t, x)
    prime_vec = zeros(1, length(t)-4 );
    for i = 1:length(prime_vec)
        c = polyfit(t(i:i+3),x(i:i+3),1);
        prime_vec(i) = c(1);
    end
end

function prime_vec = approx_deriv_1stO_legrangian(t_in, x_in)
  n = length(t_in);
  prime_vec = nan(size(x_in));

  h1 = t_in(2)-t_in(1);
  h2 = t_in(3)-t_in(2);
  prime_vec(1) = -(2*h1+h2)/(h1*(h1+h2))*x_in(1) + (h1+h2)/(h1*h2)*x_in(2) - (h1)/(h2*(h1+h2))*x_in(3);
  for i=2:(n-1)
    h1 = t_in(i)-t_in(i-1);
    h2 = t_in(i+1)-t_in(i);
    prime_vec(i) = -(h2)/(h1*(h1+h2))*x_in(i-1) - (h1-h2)/(h1*h2)*x_in(i) + (h1)/(h2*(h1+h2))*x_in(i+1);
  end
  h1 = t_in(n-1)-t_in(n-2);
  h2 = t_in(n)-t_in(n-1);
  prime_vec(n) = (h2)/(h1*(h1+h2))*x_in(n-2) - (h1+h2)/(h1*h2)*x_in(n-1) + (h1+2*h2)/(h2*(h1+h2))*x_in(n);
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
