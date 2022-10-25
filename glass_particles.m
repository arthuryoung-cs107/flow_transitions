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

    dimless_ind;
    Re_s_ns;
    cf_ns;
    G_ns;
    alpha_ns;

    alpha_tol=1.0;
    TV_range=1;
    powerfit=struct('b',NaN,'m',NaN);

    Re_s_TV=NaN;
    G_TV=NaN;
    alpha_TV=NaN;

    Re_sc1=NaN;
    Re_sc2=NaN;
    Re_sc3=NaN;
    G_c1=NaN;
    G_c2=NaN;
    G_c3=NaN;

    omega_c3_Ta=NaN;
    T_c3_Ta=NaN;

    TV_range_Ta=1;
    omega_TV_Ta=NaN;
    Re_s_TV_Ta=NaN;
    T_TV_Ta=NaN;
    G_TV_Ta=NaN;

    Re_sc1_Ta=NaN;
    Re_sc2_Ta=NaN;
    Re_sc3_Ta=NaN;
    G_c1_Ta=NaN;
    G_c2_Ta=NaN;
    G_c3_Ta=NaN;

    tau_y_fit;
    tau_y_fit_range=1:10;
    % tau_y_fit_range=1:3;

    power_fit_params;
    Carreau_fit_params=0;
    mu0_Carreau;
    lambda_Carreau;
    n_Carreau;
    k_Carreau;
    muinf_Carreau;

    Re_s_Carreau=0;
    G_Carreau=0;
    Re_s_Carreau_proper=0;
    G_Carreau_proper=0;

    G_rat_Carreau;

    S_Carreau_typ;
    mu_Carreau_typ;
    mu_effi_Carreau;
    gammai_Carreau;
    tau_Carreau;

    Re_b_Carreau;
    G_b_Carreau;
  end
  methods
    function obj = glass_particles(name_, color_)
      obj@fluid(name_, color_);
      obj.phi_m = 0.5869;
      obj.rho_p = 2500;
      obj.rho_f = 1.225;
    end
    function Re_s_out = get_Re_s(obj)
        if (obj.q<1)
            Re_s_out = obj.Re_s;
        else
            Re_s_out = obj.Re_s_Carreau;
        end
    end
    function G_out = get_G(obj)
        if (obj.q<1)
            G_out = obj.G;
        else
            G_out = obj.G_b_Carreau;  
        end
    end

    function fit_Carreau_model(obj)
        s2g=0.5*(obj.r_o+obj.r_i)/obj.r_o;
        g2s=1/s2g;

        [mu0 lam n k muinf] = obj.fit_Carreau_fluid;

        obj.mu0_Carreau = mu0;
        obj.lambda_Carreau = lam ;
        obj.n_Carreau = n ;
        obj.k_Carreau = k ;
        obj.muinf_Carreau = muinf;

        obj.Carreau_fit_params = [mu0 lam n k muinf];

        [obj.G_Carreau_proper obj.Re_s_Carreau_proper] = obj.comp_G_Res_Carreau_fluid(obj.Carreau_fit_params);
        obj.G_Carreau = obj.G_b_Carreau;
        obj.Re_s_Carreau = g2s*obj.Re_b_Carreau;

        [tau, gamma, mueff] = obj.comp_Carreau_fluid(obj.omega);
        obj.tau_Carreau = tau;
        obj.G_rat_Carreau = obj.tau./tau;

    end
    function [gamma_out, tau_out] = gamma_tau_analytical(obj,omega_)
        if (length(omega_(:))==2)
            omega=logspace(omega_(1),omega_(2),100);
        else
            omega=omega_;
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

      fit_range=obj.tau_y_fit_range;
      linfit = fit(obj.omega(fit_range), obj.tau(fit_range), 'poly1');
      obj.tau_y_fit=linfit;

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
    function string_out = label(obj)
      if round( obj.q, 1) >= 10.0
        string_out = [obj.tag, ' q=', num2str(round(obj.q, 1),2)];
      else
        string_out = [obj.tag, ' q=', num2str(round(obj.q, 1),'%.1f')];
      end
    end
    function gen_powerfit(obj)
        r_i = 0.01208;
        r_o = 0.025;
        alpha_tol = obj.alpha_tol;

        obj.alpha_T_transitions(alpha_tol, r_i, r_o)

        alpha_vec = obj.alpha_G;
        full_indices = 1:length(obj.Re_s);
        I_transitioned = logical((obj.Re_s>10).*(alpha_vec>alpha_tol));
        I_transition = min(full_indices(I_transitioned));
        TV_range = I_transition:length(obj.Re_s);
        if (length(TV_range)>=2)
            obj.TV_range = TV_range;
            obj.Re_s_TV = obj.Re_s(TV_range);
            obj.G_TV = obj.G(TV_range);
            obj.alpha_TV=alpha_vec(TV_range);

            obj.Re_sc1 = obj.Re_s_TV(1);
            obj.G_c1 = obj.G_TV(1);
            [obj.Re_sc3,obj.G_c3] = fluid.interp_trans(alpha_tol,obj.Re_s, alpha_vec, obj.G, I_transition);

            obj.powerfit = fit(reshape(obj.Re_s_TV, [], 1), reshape(obj.G_TV, [], 1),'b*x^m', 'StartPoint', [70, 1]);
            alpha = obj.powerfit.m;
            beta = obj.powerfit.b;
            m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
            obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
            obj.G_c2 = beta*(obj.Re_sc2)^alpha;
        end
    end
    function alpha_T_transitions(obj, alpha_tol_, r_i_, r_o_)
        alpha_vec = obj.alpha_T;
        full_indices = 1:length(obj.omega);
        I_transitioned = logical((obj.omega>1) .* (alpha_vec>alpha_tol_));
        I_transition = min(full_indices(I_transitioned));
        TV_range = I_transition:length(obj.omega);
        if (length(TV_range)>=2)
            obj.TV_range_Ta =TV_range;
            [Re_s_TV, I_TV] = sort(obj.Re_s_ns(obj.TV_range_Ta));
            G_TV = obj.G_ns(obj.TV_range_Ta);
            [obj.Re_sc1_Ta,ic1] = min(Re_s_TV);
            obj.omega_TV_Ta = obj.omega(obj.TV_range_Ta);
            obj.T_TV_Ta = obj.mu_torque(obj.TV_range_Ta);
            obj.Re_s_TV_Ta = Re_s_TV;
            obj.G_TV_Ta = G_TV;
            obj.G_c1_Ta = G_TV(ic1);
            [obj.Re_sc3_Ta, obj.G_c3_Ta, obj.omega_c3_Ta, obj.T_c3_Ta] = interp_trans_Ta(alpha_tol_, obj.omega, obj.Re_s_ns, alpha_vec, obj.G_ns, obj.mu_torque, I_transition);

            powerfit = fit(reshape(Re_s_TV, [], 1), reshape(G_TV, [], 1),'b*x^m', 'StartPoint', [70, 1]);
            alpha = powerfit.m;
            beta = powerfit.b;
            m = (2*pi*r_i_*r_o_)/((r_o_-r_i_)^2);
            obj.Re_sc2_Ta = exp((1/(alpha-1))*log(m/beta));
            obj.G_c2_Ta = beta*(obj.Re_sc2_Ta)^alpha;
        end
    end
    function gammai_out = comp_gammai_ro(obj,omegai_)
        mu_p = obj.mu_p;
        tau_y = obj.tau_y;
        r_o = obj.r_o;
        r_i = obj.r_i;
        gammai_out = (1/mu_p)*(tau_y+2*(omegai_*mu_p+tau_y*log(r_o/r_i))/((r_i*r_i)/(r_o*r_o)-1));
    end
    function taui_out = comp_taui_ro(obj,gammai_)
        taui_out=obj.mu_p*abs(gammai_)+obj.tau_y;
    end
    function [gammai_out,taui_out] = comp_gammai_taui_ro(obj,omegai_)
        gammai_out=obj.comp_gammai_ro(omegai_);
        taui_out=obj.comp_taui_ro(gammai_out);
    end
    function [gammai_out,taui_out] = comp_gammai_taui_rc(obj,omegai_)
        tau_y=obj.tau_y;
        mu_p=obj.mu_p;

        gammai_solve = @(t,o) (tau_y + 2*(o*mu_p + 0.5*tau_y*log(t/tau_y))./((tau_y./t)-1))/mu_p;
        bingham_solve = @(t,o) (mu_p * abs(gammai_solve(t,o))) + tau_y - t;
        % bingham_solve = @(t,o) (mu_p * gammai_solve(t,o)) + tau_y - t;

        % bounds = [tau_y,max(obj.tau)];
        bounds = [tau_y,1e3];

        gammai_out = nan(size(omegai_));
        taui_out = nan(size(omegai_));

        for i = 1:length(omegai_)
            solve_i = @(t) abs(bingham_solve(t,omegai_(i)));
            taui_out(i) = fminbnd(solve_i,bounds(1),bounds(2));
            gammai_out(i) = gammai_solve(taui_out(i),omegai_(i));
        end
    end
    function [rc_out,gammai_out,taui_out,ogt_plug_out] = comp_bingham_shear_rc(obj,omegai_)
        [gammai_out taui_out] = obj.comp_gammai_taui_rc(omegai_);
        rc_out=obj.r_i*sqrt(taui_out/obj.tau_y);

        if (nargout==4)
            full_indices=1:length(omegai_);
            shear_indices=full_indices(rc_out>obj.r_o); %% shear region
            if (isempty(shear_indices)) %% if we never reach shear flow
                index_shear=length(omegai_);
            else
                index_shear=shear_indices(1);
            end
            ogt_plug_out = [omegai_(index_shear) gammai_out(index_shear) taui_out(index_shear)];
        end
    end
    function [gammai_out,taui_out,o_g_t_plug_out] = comp_gammai_taui(obj,omegai_)
        [gammai_ro taui_ro] = obj.comp_gammai_taui_ro(omegai_);
        [rc_vec gammai_out taui_out o_g_t_plug_out] = obj.comp_bingham_shear_rc(omegai_);
        shear_inds=rc_vec>obj.r_o;
        gammai_out(shear_inds)=gammai_ro(shear_inds);
        taui_out(shear_inds)=taui_ro(shear_inds);
    end
    function [r_mat,omega_mat,gamma_mat,tau_mat] = comp_radial_shear(obj,omegai_,point_count_)
        if (nargin==2)
            point_count=30;
        else
            point_count=point_count_;
        end

        [rc_vec, g_vec, t_vec] = obj.comp_bingham_shear_rc(omegai_);

        tau_y=obj.tau_y;
        mu_p=obj.mu_p;
        r_i=obj.r_i;

        gamma_r = @(o,r,r_c) (1/mu_p)*(tau_y + (2./(r.*r))*((o*mu_p + tau_y*log(r_c/r_i))/(1/(r_c*r_c)-1/(r_i*r_i))));

        r_mat=nan(point_count,length(omegai_));
        omega_mat=nan(point_count,length(omegai_));
        gamma_mat=nan(point_count,length(omegai_));
        tau_mat=nan(point_count,length(omegai_));

        for i = 1:length(omegai_)
            rc_it = rc_vec(i);
            r_vec = linspace(r_i,rc_it,point_count);
            omega_it = omegai_(i)*ones(size(r_vec));
            gamma_it = gamma_r(omegai_(i),r_vec,rc_it);
            tau_it = mu_p*abs(gamma_it)+tau_y;

            r_mat(:,i)=reshape(r_vec,[],1);
            omega_mat(:,i)=reshape(omega_it,[],1);
            gamma_mat(:,i)=reshape(gamma_it,[],1);
            tau_mat(:,i)=reshape(tau_it,[],1);
        end
    end
    function [K_fitted,n_fitted] = fit_power_fluid(obj, ind_in)
        if (nargin==2)
            ind_=ind_in;
        else
            ind_=obj.omega<20;
        end

        omega_fit=obj.omega(ind_);
        tau_fit=obj.tau(ind_);

        eta = obj.r_i/obj.r_o;

        [K_bounds n_bounds] = determine_power_fluid_bounds(obj.omega,obj.tau,eta);
        [K_min K_0 K_max] = deal(K_bounds(1), K_bounds(2), K_bounds(3));
        [n_min n_0 n_max] = deal(n_bounds(1), n_bounds(2), n_bounds(3));

        penalty_func = @(K,n,o,t) t-K*(2*o/(n.*(1-(eta*eta).^(1.0./n)))).^(n);
        % fit_func = @(x) (penalty_func(x(1),x(2),omega_fit,tau_fit))./abs(tau_fit);
        fit_func = @(x) (penalty_func(x(1),x(2),omega_fit,tau_fit));

        x_out = lsqnonlin(fit_func,[K_0 n_0],[K_min n_min],[K_max n_max]);
        [K_fitted n_fitted] = deal(x_out(1),x_out(2));
    end
    function [taui_out gammai_out] = comp_power_fluid(obj,omegai_,K_, n_)
        if (nargin < 4)
            [K n] = obj.fit_power_fluid;
        else
            K = K_;
            n = n_;
        end
        eta = obj.r_i/obj.r_o;
        gammai_out = -1.0*((2*omegai_./(n*(1.0-(eta*eta).^(1.0./n)))).^n);
        taui_out = K*abs(gammai_out);
    end
    function [mu0_out, lambda_out, n_out, k_out, muinf_out] = fit_Carreau_fluid(obj, omega_cap_, full_flag_)
        if (nargin==1)
            omega_cap=15;
            full_flag=true;
        elseif (nargin==2)
            omega_cap=omega_cap_;
            ind=ind_;
            full_flag=true;
        elseif (nargin==3)
            omega_cap=omega_cap_;
            full_flag=full_flag_;
        end
        ind=obj.omega<omega_cap;
        omega_fit=obj.omega(ind);
        tau_fit=obj.tau(ind);
        w_fit = glass_particles.compute_distance_weighting(omega_fit);

        [mu0_out, lambda_out, n_out, k_out, muinf_out] = obj.fit_internal_Carreau_fluid(omega_fit,tau_fit,w_fit,full_flag);
    end
    function [mu0_out, lambda_out, n_out, k_out, muinf_out] = fit_internal_Carreau_fluid(obj, omega_fit, tau_fit, w_fit, full_flag)
        trust_region_reflect = 'trust-region-reflective';
        Lev_Marq = 'levenberg-marquardt';
        functol_try=1e-12;
        steptol_try=1e-12;
        optimtol_try=1e-12;

        % alg_use=trust_region_reflect;
        alg_use=Lev_Marq;
        optimtol_use=optimtol_try;
        functol_use=functol_try;
        steptol_use=steptol_try;

        penalty_func = @(m0,l,n,k,minf,o,t,w) w.*(t-((minf+(m0-minf)*((1+(l*k*o).*(l*k*o)).^(0.5*(n-1)))).*(k*o)));
        [mu0_b lam_b n_b k_b muinf_b] = obj.determine_Carreau_fluid_bounds(omega_fit,tau_fit);
        [m0_min m0_0 m0_max] = deal(mu0_b(1), mu0_b(2), mu0_b(3));
        [l_min l_0 l_max] = deal(lam_b(1), lam_b(2), lam_b(3));
        [n_min n_0 n_max] = deal(n_b(1), n_b(2), n_b(3));

        [k_out muinf_out] = deal(2*(obj.r_o*obj.r_o)/(obj.r_o*obj.r_o-obj.r_i*obj.r_i),0);
        if (full_flag)
            [k_min k_0 k_max] = deal(k_b(1), k_b(2), k_b(3));
            [minf_min minf_0 minf_max] = deal(muinf_b(1), muinf_b(2), muinf_b(3));

            fit_func = @(x) (penalty_func(x(1),x(2),x(3),x(4),x(5),omega_fit,tau_fit,w_fit));

            upper_bound=[m0_max l_max n_max k_max minf_max];
            lower_bound=[m0_min l_min n_min k_min minf_min];
            init_guess=[m0_0 l_0 n_0 k_0 minf_0];
        else
            fit_func = @(x) (penalty_func(x(1),x(2),x(3),k_out,muinf_out,omega_fit,tau_fit,w_fit));

            upper_bound=[m0_max l_max n_max];
            lower_bound=[m0_min l_min n_min];
            init_guess=[m0_0 l_0 n_0];
        end

        % solve_opts=optimoptions(@lsqnonlin,'Display','final', ...
        solve_opts=optimoptions(@lsqnonlin,'Display','off', ...
                                'Algorithm',alg_use, ...
                                'MaxIterations', 1e5, ...
                                'MaxFunctionEvaluations', 1e5, ...
                                'FiniteDifferenceType', 'central', ...
                                'OptimalityTolerance',optimtol_use, ...
                                'FunctionTolerance',functol_use, ...
                                'StepTolerance', steptol_use);

        x_out = lsqnonlin(fit_func, init_guess, ...
                                    lower_bound, ...
                                    upper_bound, ...
                                    solve_opts);
        [mu0_out lambda_out n_out] = deal(x_out(1),x_out(2),x_out(3));
        if (full_flag)
            [k_out,muinf_out] = deal(x_out(4),x_out(5));
        end
    end
    function [taui_out gammai_out mueffi_out] = comp_Carreau_fluid(obj,omegai_,mu0_,l_,n_,k_,muinf_)
        if (nargin == 2)
            if (length(obj.Carreau_fit_params)==1)
                [mu0,l,n,k,muinf] = obj.fit_Carreau_fluid;
            else
                par=obj.Carreau_fit_params;
                [mu0,l,n,k,muinf] = deal(par(1),par(2),par(3),par(4),par(5));
            end
        else
            [mu0,l,n,k,muinf] = deal(mu0_,l_,n_,k_,muinf_);
        end
        gammai_out = -1.0*(k*omegai_);
        mueffi_out = muinf+((mu0-muinf)*((1+((l*gammai_out).*(l*gammai_out))).^(0.5*(n-1))));
        taui_out = mueffi_out.*abs(gammai_out);
    end
    function [G_out Re_s_out] = comp_G_Res_Carreau_fluid(obj,Carreau_par_)
        if (nargin == 2)
            [mu0,l,n,k,muinf] = obj.fit_Carreau_fluid;
        else
            [mu0,l,n,k,muinf] = deal(Carreau_par_(1),Carreau_par_(2),Carreau_par_(3),Carreau_par_(4),Carreau_par_(5));
        end
        [taui gi mueffi] = obj.comp_Carreau_fluid(obj.omega,mu0,l,n,k,muinf);
        [r_i r_o rho_b h]=deal(obj.r_i,obj.r_o,obj.rho_b,obj.h);
        tau_typ = (r_i/r_o)*taui;

        %% solving for the local shear rate at the specified shear stress. k does not come into play
        mu_solve = @(s) muinf+(mu0-muinf)*((1+(l*s).*(l*s)).^(0.5*(n-1)));
        Carreau_solve = @(t,s) (mu_solve(s).*(s)) - t;

        solve_opts=optimset('Display', 'off', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'TolX', 1e-14);

        S=nan(size(tau_typ));
        for i = 1:length(tau_typ)
            solve_i = @(s) abs(Carreau_solve(tau_typ(i),s));
            S(i)=fminbnd(solve_i,0,abs(gi(i)),solve_opts);
        end
        mu_use = mu_solve(S);
        obj.S_Carreau_typ = S;
        obj.mu_Carreau_typ = mu_use;
        obj.mu_effi_Carreau = mueffi;
        obj.gammai_Carreau = gi;

        Re_s_out = ((rho_b*(r_o-r_i)*(r_o-r_i))./mu_use).*(S);
        G_out = (rho_b./(h.*mu_use.*mu_use)).*obj.mu_torque;

        obj.Re_b_Carreau = (rho_b*r_i*(r_o-r_i))./mueffi.*(obj.omega);
        obj.G_b_Carreau = (rho_b./(h.*mueffi.*mueffi)).*obj.mu_torque;

    end
    function [mu0_bout lambda_bout n_bout k_bout muinf_bout] = determine_Carreau_fluid_bounds(obj,o_,t_)
        k_Newt = 2*(obj.r_o*obj.r_o)/(obj.r_o*obj.r_o-obj.r_i*obj.r_i);
        [t_max i_tmax] = max(t_);
        tau_full_max = max(obj.tau);
        o_full_max = max(obj.omega);
        t_min=min(t_);
        o_min=min(o_);
        o_max=max(o_);

        mu0_min = t_min/(k_Newt*o_full_max);
        lambda_min=1/(k_Newt*o_max);
        n_min=0;
        k_min=k_Newt;
        % k_min=floor(k_Newt);
        muinf_min=0;

        mu0_max = (tau_full_max)/(k_Newt*o_min);
        lambda_max=1e5;
        n_max=1;
        k_max=k_Newt;
        % k_max=ceil(k_Newt);
        % muinf_max=mu0_min;
        muinf_max=0;

        mu0_0=t_(1)/(k_Newt*o_(1));;
        lambda_0=1.0/(k_Newt*o_min);
        n_0=0.5;
        k_0=k_Newt;
        muinf_0=0;

        mu0_bout = [mu0_min mu0_0 mu0_max];
        lambda_bout = [lambda_min lambda_0 lambda_max];
        n_bout = [n_min n_0 n_max];
        k_bout = [k_min k_0 k_max];
        muinf_bout = [muinf_min muinf_0 muinf_max];
    end
  end
  methods (Static)
    function w_out = compute_distance_weighting(o_)
        % len_o=length(o_);
        % w_first=o_(2)-o_(1);
        % w_last=o_(len_o)-o_(len_o-1);
        % w_out = [w_first; 0.5*(o_(3:end)-o_(1:(len_o-2))); w_last];
        % w_out = w_out/norm(w_out);

        % w_out = ones(size(o_))/norm(ones(size(o_)));

        del=1/(length(o_)-1)*log(o_(end)/o_(1));
        w_out=del.^(0:(length(o_)-1));
        w_out=w_out/norm(w_out);
    end
  end
end

function [K_bounds_out, n_bounds_out] = determine_power_fluid_bounds(o_,t_,eta_)
    K_min = 0.0;
    n_max = 1.0;
    K_max = 10*max(t_);
    n_min = 1e-10;

    K_bounds_out = [K_min 0.5*(K_min+K_max) K_max];
    n_bounds_out = [n_min 0.5*(n_min+n_max) n_max];
end

function [Rc, Gc, omegac, Tc] = interp_trans_Ta(alpha_tol_,omega_,R_,alpha_,G_,T_,It_)
    R2=R_(It_);
    omega2=omega_(It_);
    alpha2=alpha_(It_);
    G2=G_(It_);
    T2=T_(It_);

    R1=R_(It_-1);
    omega1=omega_(It_-1);
    alpha1=alpha_(It_-1);
    G1=G_(It_-1);
    T1=T_(It_-1);

    omegac=(alpha_tol_-alpha2)*((omega2-omega1)/(alpha2-alpha1)) + omega2;
    Rc=(omegac-omega2)*((R2-R1)/(omega2-omega1)) + R2;
    Gc=(omegac-omega2)*((G2-G1)/(omega2-omega1)) + G2;
    Tc=(omegac-omega2)*((T2-T1)/(omega2-omega1)) + T2;
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
