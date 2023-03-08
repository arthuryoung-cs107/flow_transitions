classdef experiment < handle
    properties (Constant)
        UB1_Bingham_tauy_bounds = [1e-1 1 1e1];
        UB2_Bingham_tauy_bounds = [1e-1 1 1e1];
        XB1_Bingham_tauy_bounds = [1e-1 1 1e1];
        XB2_Bingham_tauy_bounds = [1e-1 1 1e1];

        UB1_Bingham_mup_bounds = [1e-2 5e-2 1e-1];
        UB2_Bingham_mup_bounds = [1e-2 5e-2 1e-1];
        XB1_Bingham_mup_bounds = [1e-2 5e-2 1e-1];
        XB2_Bingham_mup_bounds = [1e-2 5e-2 1e-1];

        UB1_Bingham_omega_cap = 7.3;
        UB2_Bingham_omega_cap = 5.2;
        XB1_Bingham_omega_cap = 8.5;
        XB2_Bingham_omega_cap = 9.5;

        UB1_Bingham_omega_floor = 0.62;
        UB2_Bingham_omega_floor = 0.41;
        XB1_Bingham_omega_floor = 0.31;
        XB2_Bingham_omega_floor = 0.41;
    end
    properties
        label;
        color;
        specs;
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
        mu_eff;

        mu_torque;
        sigma_torque;
        G;
        G_rat;
        cf;

        omega;
        Re_s;

        alpha_tol=1.3;
        TV_range=1;
        powerfit=struct('b',NaN,'m',NaN);
        powerfit_Grat_Res=struct('b',NaN,'m',NaN);

        Re_s_TV=NaN;
        G_TV=NaN;
        Grat_TV=NaN;
        alpha_TV=NaN;

        % Re_sc=NaN;
        Re_sc1=NaN;
        Re_sc2=NaN;
        Re_sc3=NaN;

        G_c1=NaN;
        G_c2=NaN;
        G_c3=NaN;

        tau_qs;
        omega_fit_Bingham;
        tau_fit_Bingham;
        i_fit_Bingham;

        Bingham_tauy_bounds;
        Bingham_mup_bounds;
        Bingham_omega_cap=10;
        Bingham_omega_floor=0;
        Bingham_wscheme='tmagnorm_odist1';

        tau_y_Bingham;
        mu_p_Bingham;
        gamma_Bingham;
        Re_b_Bingham;
        Re_s_Bingham;
        S_Bingham;
        cf_Bingham;
        G_b_Bingham;
        G_rat_Bingham;
        rc_Bingham;

        exp;
        len; %% this is the number of data sets
        dat_num; %% this is the length of the measurement sets
    end
  methods
    function obj = experiment(exp_list_in)
      obj.exp = exp_list_in;
      obj.len = length(exp_list_in);
    end
    function rho_rel_out = comp_rho_rel(obj)
        rho_b = obj.comp_rho_b(obj.phi_m);
        rho_rel_out = rho_b - obj.rho_f;
    end
    function rho_b_out = comp_rho_b(obj, phi_)
        if (nargin==1)
            phi = obj.phi;
        else
            phi = phi_;
        end
        rho_b_out = obj.rho_p * phi + obj.rho_f*(1.0-phi);
    end
    function [Res_out Grat_out] = comp_Grat_Res_KD(obj,phi_)
        [r_i r_o h] = deal(fluid.r_i_def, fluid.r_o_def, fluid.h_def);
        rho_b = obj.comp_rho_b(phi_);
        mu = obj.Krieger_Dougherty(obj.mu_f,phi_,obj.phi_m);
        nu = mu/rho_b;
        gammai = 2*r_o*r_o*(obj.omega)/(r_o*r_o - r_i*r_i);
        [b eta] = deal(r_o-r_i,r_i/r_o);
        Res_out = b*b*eta*(gammai)/nu;
        Grat_out = obj.mu_torque./(2*pi*r_i*r_i*h*mu*gammai);
    end
    function tau_out = tau_comp(obj)
        tau_out = obj.mu_torque/(2*pi*(fluid.r_i_def*fluid.r_i_def)*fluid.h_def);
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
    function fit_Bingham_model(obj)
        tau_full = obj.tau_comp;
        ind = logical(double(obj.omega < obj.Bingham_omega_cap).*double(obj.omega > obj.Bingham_omega_floor));
        [omega_fit tau_fit] = deal(obj.omega(ind),tau_full(ind));

        obj.tau_qs = mean(tau_fit);

        [obj.omega_fit_Bingham obj.tau_fit_Bingham obj.i_fit_Bingham] = deal(omega_fit, tau_fit,ind);
        w_fit = glass_particles.compute_equidistant_weighting(omega_fit, tau_fit, obj.Bingham_wscheme);

        [mp_b ty_b] = deal(obj.Bingham_mup_bounds, obj.Bingham_tauy_bounds);

        [obj.mu_p_Bingham obj.tau_y_Bingham] = glass_particles.fit_internal_Bingham_fluid(omega_fit, tau_fit, w_fit, mp_b, ty_b);

        ti = glass_particles.taui_pred_Bingham(obj.mu_p_Bingham,obj.tau_y_Bingham,obj.omega);
        obj.gamma_Bingham = (ti-obj.tau_y_Bingham)/obj.mu_p_Bingham;
        obj.rc_Bingham = glass_particles.determine_rc_Bingham(obj.mu_p_Bingham, obj.tau_y_Bingham, obj.omega);

        [o_ rc] = deal(obj.omega, obj.rc_Bingham);
        [h ri ro] = deal(fluid.h_def, fluid.r_i_def, fluid.r_o_def);
        [mp_ ty_] = deal(obj.mu_p_Bingham, obj.tau_y_Bingham);
        rho_b = obj.comp_rho_b;

        rc2 = rc.*rc;
        rt = sqrt(ri * rc);
        rt2 = rt.*rt;
        etac = ri./rc;
        etac2 = etac.*etac;
        d = rc-ri;

        obj.Re_b_Bingham = rho_b*(obj.omega*ri*(ro-ri))/mp_;
        obj.S_Bingham = (((2./rt2).*((mp_*o_-ty_*log(etac)))./(1/(ri*ri) - 1./(rc2)))-ty_)/mp_;
        obj.Re_s_Bingham = rho_b*(obj.S_Bingham.*(d.*d))/mp_;
        obj.cf_Bingham = obj.mu_torque./(2*pi*ri*ri*h*obj.S_Bingham.*obj.S_Bingham.*d.*d);
        obj.G_b_Bingham = (rho_b/(h*mp_*mp_))*obj.mu_torque;
        obj.G_rat_Bingham = (obj.mu_torque.*(rc.*rc-ri*ri))./(4*pi*ri*ri*h*(rc.*rc).*(mp_*o_ + ty_*log(rc/ri)));

    end
    function gen_powerfit(obj)
      r_i = 0.01208;
      r_o = 0.025;
      alpha_tol = obj.alpha_tol;

      alpha_vec = obj.alpha;
      full_indices = 1:length(obj.omega);
      I_transitioned = logical((obj.Re_s>50) .* (alpha_vec>alpha_tol));
      I_transition = min(full_indices(I_transitioned));
      obj.TV_range = I_transition:length(obj.omega);
      Re_s_TV = obj.Re_s(obj.TV_range);
      obj.Re_sc1 = Re_s_TV(1);
      obj.Re_s_TV = Re_s_TV;
      obj.G_TV = obj.G(obj.TV_range);
      obj.alpha_TV = alpha_vec(obj.TV_range);

      obj.G_c1 = obj.G_TV(1);
      [obj.Re_sc3, obj.G_c3] = fluid.interp_trans(alpha_tol, obj.Re_s, alpha_vec, obj.G, I_transition);

      obj.powerfit = fit(reshape(Re_s_TV, [], 1), reshape(obj.G_TV,[],1),'b*x^m', 'StartPoint', [70, 1]);
      alpha = obj.powerfit.m;
      beta = obj.powerfit.b;
      m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
      obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
      obj.G_c2 = beta*(obj.Re_sc2)^alpha;
    end
    function gen_powerfit_internal(obj, Re_sc_floor_)
        [r_i r_o] = deal(fluid.r_i_def,fluid.r_o_def);
        alpha_full = obj.alpha;
        ifull = 1:length(obj.Re_s);
        iTV = obj.Re_s > Re_sc_floor_;
        obj.TV_range = ifull(iTV);
        obj.Re_s_TV = obj.Re_s(obj.TV_range);
        obj.G_TV = obj.G(obj.TV_range);
        obj.Grat_TV = obj.G_rat(obj.TV_range);
        obj.alpha_TV = alpha_full(obj.TV_range);
        obj.Re_sc1 = obj.Re_s_TV(1);
        obj.G_c1 = obj.G_TV(1);

        [obj.Re_sc3, obj.G_c3] = fluid.interp_trans(obj.alpha_tol, obj.Re_s, alpha_full, obj.G, obj.TV_range(1));

        Res_fit = reshape(obj.Re_s_TV,[],1);
        G_fit = reshape(obj.G_TV,[],1);
        Grat_fit = reshape(obj.Grat_TV,[],1);

        obj.powerfit = fit(Res_fit,G_fit,'b*x^m', 'StartPoint', [70, 1]);
        obj.powerfit_Grat_Res = fit(Res_fit,Grat_fit,'b*x^m', 'StartPoint', [1/(sqrt(70)), 0.5]);
        alpha = obj.powerfit.m;
        beta = obj.powerfit.b;
        m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
        obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
        obj.G_c2 = beta*(obj.Re_sc2)^alpha;
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
    function Re_b_out = Re_b_comp(obj)
        Re_b_out = ((fluid.r_o_def + fluid.r_i_def)/(2*fluid.r_o_def))*obj.Re_s;
    end
    function o_out = dtdo_o(obj)
        os = mink(obj.omega,length(obj.omega));
        o_out = 0.5*(os(1:(length(os)-1))+os(2:end));
    end
    function o_out = d2tdo2_o(obj)
        os = mink(obj.omega,length(obj.omega));
        o_out = os(2:(length(os)-1));
    end
    function d_out = dtdo_d(obj)
        d_out = fluid.comp_deriv_central(obj.omega,obj.tau_comp);
    end
    function d_out = d2tdo2_d(obj)
        d_out = fluid.comp_2deriv_central(obj.omega,obj.tau_comp);
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

  % k = 5; %% number of points considered. Must be odd
  k = 3; %% number of points considered. Must be odd
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
