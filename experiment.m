classdef experiment < handle
  properties
    label;
    color;
    specs;
    LW = 1.0;
    MS = 5.0; 

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
      n = length(obj.Re_s);
      alpha_out = nan(size(obj.cf));
      logcf = log(obj.cf);
      logRe = log(obj.Re_s);

      h1 = logRe(2)-logRe(1);
      h2 = logRe(3)-logRe(2);
      alpha_out(1) = -(2*h1+h2)/(h1*(h1+h2))*logcf(1) + (h1+h2)/(h1*h2)*logcf(2) - (h1)/(h2*(h1+h2))*logcf(3);
      for i=2:(n-1)
        h1 = logRe(i)-logRe(i-1);
        h2 = logRe(i+1)-logRe(i);
        alpha_out(i) = -(h2)/(h1*(h1+h2))*logcf(i-1) - (h1-h2)/(h1*h2)*logcf(i) + (h1)/(h2*(h1+h2))*logcf(i+1);
      end
      h1 = logRe(n-1)-logRe(n-2);
      h2 = logRe(n)-logRe(n-1);
      alpha_out(n) = (h2)/(h1*(h1+h2))*logcf(n-2) - (h1+h2)/(h1*h2)*logcf(n-1) + (h1+2*h2)/(h2*(h1+h2))*logcf(n);
      alpha_out = alpha_out + 2;
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
      figure
      hold on
      set(gca, 'YScale', 'log')
      set(gca, 'XScale', 'log')
      for i=1:obj.len
        plot(obj.exp(i).mu_rpm, obj.exp(i).mu_torque, ' o','Color', obj.exp(i).color, 'LineWidth', 1.5, 'DisplayName', obj.exp(i).name)
      end
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
