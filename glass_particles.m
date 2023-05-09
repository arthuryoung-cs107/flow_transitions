classdef glass_particles < fluid
  properties (Constant)

    FB1_EC_Bingham_tauy_bounds =        [   1.0e-8 9.20000000000000e+000 5.000e+2; ...  % 1.0e-8 - 5.000e+2
                                            1.0e-8 9.20000000000000e+000 5.000e+2; ...  % 1.0e-8 - 5.000e+2
                                            1.0e-8 9.20000000000000e+000 5.000e+2; ...  % 1.0e-8 - 5.000e+2
                                            1.0e-8 9.20000000000000e+000 5.000e+2];     % 1.0e-8 - 5.000e+2

    FB1_EC_Bingham_mup_bounds =        [    1.1e-2 25.7394567553242e-003 1.000e-0; ...  % 1.0e-3 - 1.000e-0
                                            1.2e-2 25.7394567553242e-003 1.000e-0; ...  % 1.0e-3 - 1.000e-0
                                            1.0e-2 25.7394567553242e-003 1.000e-0; ...  % 1.0e-3 - 1.000e-0
                                            1.0e-2 25.7394567553242e-003 3.500e-2];     % 1.0e-3 - 3.500e-2

    FB1_EC_Bingham_omega_cap =         [    1.00e3; ...
                                            5.20e1; ...
                                            3.20e1; ... % 5.20e1
                                            3.20e1];    % 5.20e1

    FB1_EC_Bingham_omega_floor =       [    0.0e+0; ...
                                            0.0e+0; ...
                                            0.0e+0; ...
                                            0.0e+0];


    FB1_Bingham_tauy_bounds_incthin =   [   1.0e-8 9.20000000000000e+000 1.000e+2; ...  % 9.200e-0
                                            1.0e-8 7.07546425350865e+000 1.000e+2; ...  % 7.200e-0
                                            1.0e-8 5.90392839762222e+000 1.000e+2; ...  % 6.300e-0
                                            1.0e-8 3.11100000000000e+000 3.111e-0; ...  %
                                            1.0e-8 1.93300000000000e+000 1.933e-0; ...  %
                                            1.0e-8 722.067580588442e-003 8.389e-1; ...  %
                                            1.0e-8 355.600000000000e-003 3.556e-1; ...  %
                                            1.0e-8 121.800000000000e-003 1.218e-1; ...  %
                                            1.0e-8 57.2400000000000e-003 5.724e-2; ...  %
                                            1.0e-8 41.5959091010154e-003 4.315e-2];     %
    FB2_Bingham_tauy_bounds_incthin =   [   1.0e-8 157.451268396902e+000 160.0e-0; ...  %
                                            1.0e-8 120.706722323418e+000 130.0e-0; ...  %
                                            1.0e-8 80.6644129599834e+000 87.00e-0; ...  %
                                            1.0e-8 58.6178522358182e+000 65.00e-0; ...  %
                                            1.0e-8 6.37020357237319e+000 6.600e-0; ...  %
                                            1.0e-8 1.93842075148303e+000 3.939e-0; ...  %
                                            1.0e-8 909.800000000000e-003 9.098e-1; ...  %
                                            1.0e-8 263.600000000000e-003 2.636e-1; ...  %
                                            1.0e-8 122.300000000000e-003 1.223e-1; ...  %
                                            1.0e-8 77.9200000000000e-003 7.792e-2; ...  %
                                            1.0e-8 30.8700000000000e-003 3.087e-2; ...  %
                                            1.0e-8 24.2600000000000e-003 2.426e-2; ...  %
                                            1.0e-8 17.4700000000000e-003 1.747e-2];     %
    FB1_Bingham_mup_bounds_incthin  =   [   1.0e-3 25.7394567553242e-003 2.580e-2; ...
                                            1.0e-3 44.8079134410635e-003 1.000e-0; ...
                                            1.0e-3 51.4140906361686e-003 1.000e-0; ...
                                            1.0e-3 119.865050767594e-003 1.000e-0; ...
                                            1.0e-3 168.020382361456e-003 1.000e-0; ...
                                            1.0e-3 232.140859326886e-003 1.000e-0; ...
                                            1.0e-3 201.153788881823e-003 1.000e-0; ...
                                            1.0e-3 207.576695783095e-003 1.000e-0; ...
                                            1.0e-3 206.794429294975e-003 1.000e-0; ...
                                            1.0e-3 200.229074066884e-003 1.000e-0];
    FB2_Bingham_mup_bounds_incthin  =   [   1.0e-3 18.3260850877242e-003 1.000e-0; ...
                                            1.0e-3 20.9762887171971e-003 1.000e-0; ...
                                            1.0e-3 39.2248738646570e-003 1.000e-0; ...
                                            1.0e-3 61.3418531437856e-003 1.000e-0; ...
                                            1.0e-3 92.0000000000000e-003 0.092e-0; ...
                                            1.0e-3 92.7469606772215e-003 1.000e-0; ...
                                            1.0e-3 82.1356176378636e-003 1.000e-0; ...
                                            1.0e-3 87.1812671974613e-003 1.000e-0; ...
                                            1.0e-3 89.5329959077820e-003 1.000e-0; ...
                                            1.0e-3 89.5696228621380e-003 1.000e-0; ...
                                            1.0e-3 88.3735319414805e-003 1.000e-0; ...
                                            1.0e-3 88.8392010324877e-003 1.000e-0; ...
                                            1.0e-3 88.9701870884362e-003 1.000e-0];
    FB1_Bingham_omega_cap_incthin = [       2.1e1; ... % 1.1e1, 2.1e1
                                            2.1e1; ... % 3.0e1, 3.2e1
                                            2.1e1; ... % 3.0e1, 3.2e1
                                            3.0e1; ...
                                            3.0e1; ...
                                            3.0e1; ...
                                            3.0e1; ...
                                            3.0e1; ...
                                            3.0e1; ...
                                            3.0e1];
    FB2_Bingham_omega_cap_incthin = [       8.0e1; ...
                                            7.1e1; ...
                                            7.1e1; ...
                                            5.0e1; ...
                                            3.1e1; ...
                                            3.1e1; ...
                                            3.1e1; ...
                                            2.8e1; ...
                                            2.8e1; ...
                                            2.8e1; ...
                                            2.8e1; ...
                                            2.8e1; ...
                                            2.8e1];
    FB1_Bingham_wscheme_incthin = {         'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'};
    FB2_Bingham_wscheme_incthin = {         'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'};

    FB1_Bingham_tauy_bounds_prethin =   [   1.0e-8 9.20000000000000e+000 1.000e+3; ...
                                            1.0e-8 7.07546425350865e+000 1.000e+3; ...
                                            1.0e-8 5.90392839762222e+000 1.000e+3; ...
                                            1.0e-8 3.11100000000000e+000 1.000e+3; ...
                                            1.0e-8 1.93300000000000e+000 1.000e+3; ...
                                            1.0e-8 722.067580588442e-003 1.000e+3; ...
                                            1.0e-8 355.600000000000e-003 1.000e+3; ...
                                            1.0e-8 121.800000000000e-003 1.000e+3; ...
                                            1.0e-8 57.2400000000000e-003 1.000e+3; ...
                                            1.0e-8 41.5959091010154e-003 1.000e+3];
    FB2_Bingham_tauy_bounds_prethin =   [   1.0e-8 157.451268396902e+000 1.000e+3; ...
                                            1.0e-8 120.706722323418e+000 1.000e+3; ...
                                            1.0e-8 80.6644129599834e+000 1.000e+3; ...
                                            1.0e-8 58.6178522358182e+000 1.000e+3; ...
                                            1.0e-8 6.37020357237319e+000 1.000e+3; ...
                                            1.0e-8 1.93842075148303e+000 1.000e+3; ...
                                            1.0e-8 909.800000000000e-003 1.000e+3; ...
                                            1.0e-8 263.600000000000e-003 1.000e+3; ...
                                            1.0e-8 122.300000000000e-003 1.000e+3; ...
                                            1.0e-8 77.9200000000000e-003 1.000e+3; ...
                                            1.0e-8 30.8700000000000e-003 1.000e+3; ...
                                            1.0e-8 24.2600000000000e-003 1.000e+3; ...
                                            1.0e-8 17.4700000000000e-003 1.000e+3];
    FB1_Bingham_mup_bounds_prethin  =   [   1.0e-3 25.7394567553242e-003 1.000e+2; ...
                                            1.0e-3 44.8079134410635e-003 1.000e+2; ...
                                            1.0e-3 51.4140906361686e-003 1.000e+2; ...
                                            1.0e-3 119.865050767594e-003 1.000e+2; ...
                                            1.0e-3 168.020382361456e-003 1.000e+2; ...
                                            1.0e-3 232.140859326886e-003 1.000e+2; ...
                                            1.0e-3 201.153788881823e-003 1.000e+2; ...
                                            1.0e-3 207.576695783095e-003 1.000e+2; ...
                                            1.0e-3 206.794429294975e-003 1.000e+2; ...
                                            1.0e-3 200.229074066884e-003 1.000e+2];
    FB2_Bingham_mup_bounds_prethin  =   [   1.0e-3 18.3260850877242e-003 1.000e+2; ...
                                            1.0e-3 20.9762887171971e-003 1.000e+2; ...
                                            1.0e-3 39.2248738646570e-003 1.000e+2; ...
                                            1.0e-3 61.3418531437856e-003 1.000e+2; ...
                                            1.0e-3 92.0000000000000e-003 1.000e+2; ...  % 455.588462950813e-003
                                            1.0e-3 92.7469606772215e-003 1.000e+2; ...
                                            1.0e-3 82.1356176378636e-003 1.000e+2; ...
                                            1.0e-3 87.1812671974613e-003 1.000e+2; ...
                                            1.0e-3 89.5329959077820e-003 1.000e+2; ...
                                            1.0e-3 89.5696228621380e-003 1.000e+2; ...
                                            1.0e-3 88.3735319414805e-003 1.000e+2; ...
                                            1.0e-3 88.8392010324877e-003 1.000e+2; ...
                                            1.0e-3 88.9701870884362e-003 1.000e+2];
    FB1_Bingham_omega_cap_prethin = [       2.10e1; ... % 1
                                            2.10e1; ... % 2
                                            2.10e1; ... % 3
                                            1.10e0; ... % 4
                                            1.10e0; ... % 5
                                            1.10e0; ... % 6
                                            1.10e0; ... % 7
                                            1.10e0; ... % 8
                                            1.10e0; ... % 9
                                            1.10e0];    % 10
    FB2_Bingham_omega_cap_prethin = [       9.50e0; ... % 1
                                            7.01e1; ... % 2
                                            3.91e1; ... % 3
                                            9.50e0; ... % 4
                                            9.50e0; ... % 5
                                            0.06e0; ... % 6
                                            0.06e0; ... % 7
                                            0.09e0; ... % 8
                                            0.18e0; ... % 9
                                            0.57e0; ... % 10
                                            1.00e0; ... % 11
                                            1.00e0; ... % 12
                                            1.00e0];    % 13
    FB1_Bingham_wscheme_prethin = {         'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'};
    FB2_Bingham_wscheme_prethin = {         'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'; ...
                                            'tmagnorm_odist1'};
    FB1_EC_phi_low = [                      0.58; ... % 0.58; ...
                                            0.57; ... % 0.57; ...
                                            0.57; ... % 0.57; ...
                                            0.56];

    FB1_phi_low = [                         0.55; ... % 0.55; ...
                                            0.55; ... % 0.55; ...
                                            0.55; ... % 0.55; ...
                                            0.54; ... % 0.54; ...
                                            0.54; ... % 0.54; ...
                                            0.53; ... % 0.53; ...
                                            0.53; ... % 0.53; ...
                                            0.52; ... % 0.52; ...
                                            0.51; ... % 0.51; ...
                                            0.51];

    FB2_phi_low = [                         0.58; ... % 0.58; ...
                                            0.57; ... % 0.57; ...
                                            0.56; ... % 0.56; ...
                                            0.56; ... % 0.56; ...
                                            0.58; ... % 0.55; ...
                                            0.58; ... % 0.54; ...
                                            0.58; ... % 0.53; ...
                                            0.58; ... % 0.52; ...
                                            0.51; ... % 0.51; ...
                                            0.50; ... % 0.50; ...
                                            0.49; ... % 0.49; ...
                                            0.48; ... % 0.48; ...
                                            0.47];

  end
  properties
    exp_id=0;
    tag = 'FB';
    Q;
    q;
    q_inc;
    tau_static;

    FB_fitted_Bingham_pars;

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

    omega_fit_Carreau;
    tau_fit_Carreau;

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

    FB_Bingham_tauy_bounds;
    FB_Bingham_mup_bounds;
    FB_Bingham_omega_cap;
    FB_Bingham_wscheme;
    mu_p_Bingham_intr=NaN;
    alpha_tol_Bingham =3.0;
    omega_crit_min = 1;
    omega_crit=NaN;
    omega_crit_alpha_peak=NaN;
    omega_cap_Bingham_use=NaN;
    omega_fit_Bingham;
    tau_fit_Bingham;
    i_fit_Bingham;

    tau_y_Bingham;
    mu_p_Bingham;
    gamma_Bingham;
    S_Bingham;
    Re_s_Bingham;
    Re_b_Bingham;
    cf_Bingham;
    G_b_Bingham;
    G_rat_Bingham;
    rc_Bingham;

    Carreau_fit_muinf_params;
    omega_fit_Carreau_muinf;
    tau_fit_Carreau_muinf;

    qcrit_sgf;
    qcrit_ttv;
    qcrit_fgm;

    ocrit_dgf_ql;
    ocrit_dgf_qh;
    ocrit_fgm_ql;
    ocrit_fgm_qh;

    qcrit_Bingham2Carreau;
    FB_phi_low;
  end
  methods
    function obj = glass_particles(name_, color_)
      obj@fluid(name_, color_);
      obj.phi_m = 0.5869;
      obj.rho_p = 2500;
      obj.rho_f = 1.225;
    end
    function rho_rel_out = comp_rho_rel(obj)
        rho_b = obj.rho_p * obj.FB_phi_low(obj.exp_id) + obj.rho_f*(1 - obj.FB_phi_low(obj.exp_id));
        rho_rel_out = rho_b - obj.rho_f;
    end
    function Re_b_out = get_Re_b(obj)
        if (obj.q<obj.qcrit_Bingham2Carreau)
            Re_b_out = obj.Re_b_Bingham;
        else
            Re_b_out = obj.Re_b_Carreau;
        end
    end
    function Re_s_out = get_Re_s(obj)
        if (obj.q<obj.qcrit_Bingham2Carreau)
            Re_s_out = obj.Re_s_Bingham;
        else
            Re_s_out = obj.Re_s_Carreau;
        end
    end
    function G_out = get_G(obj)
        if (obj.q<obj.qcrit_Bingham2Carreau)
            G_out = obj.G_b_Bingham;
        else
            G_out = obj.G_b_Carreau;
        end
    end
    function [mu0_out, lambda_out, n_out, k_out, muinf_out] = fit_Carreau_model_muinf(obj)
        if (isempty(strfind(obj.name,'49'))) % this is FB1
            omega_cap_vec = [   2.0e1; ... % 1.1e1, 2.1e1
                                2.0e1; ... % 3.0e1, 3.2e1
                                2.0e1; ... % 3.0e1, 3.2e1
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1];

        else
            omega_cap_vec = [   2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1; ...
                                2.0e1];
        end
        omega_cap = omega_cap_vec(obj.exp_id);
        i_fit_Carreau_muinf = obj.omega_fit_Bingham<omega_cap;
        omega_fit = obj.omega(i_fit_Carreau_muinf);
        tau_fit = obj.tau(i_fit_Carreau_muinf);
        w_fit = ones(size(omega_fit));

        obj.omega_fit_Carreau_muinf = omega_fit;
        obj.tau_fit_Carreau_muinf = tau_fit;

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
        [mu0_b lam_b n_b k_b muinf_b] = obj.determine_Carreau_fluid_bounds_muinf(omega_fit,tau_fit);
        [m0_min m0_0 m0_max] = deal(mu0_b(1), mu0_b(2), mu0_b(3));
        [l_min l_0 l_max] = deal(lam_b(1), lam_b(2), lam_b(3));
        [n_min n_0 n_max] = deal(n_b(1), n_b(2), n_b(3));
        [k_min k_0 k_max] = deal(k_b(1), k_b(2), k_b(3));
        [minf_min minf_0 minf_max] = deal(muinf_b(1), muinf_b(2), muinf_b(3));

        fit_func = @(x) (penalty_func(x(1),x(2),x(3),x(4),x(5),omega_fit,tau_fit,w_fit));

        upper_bound=[m0_max l_max n_max k_max minf_max];
        lower_bound=[m0_min l_min n_min k_min minf_min];
        init_guess=[m0_0 l_0 n_0 k_0 minf_0];

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
        if (nargout == 1)
            mu0_out = x_out;
        else
            [mu0_out lambda_out n_out k_out muinf_out] = deal(x_out(1),x_out(2),x_out(3),x_out(4),x_out(5));
        end
    end
    function fit_Bingham_model(obj,omega_cap_,omega_fit_,tau_fit_)
        ind = omega_fit_<omega_cap_;
        omega_fit = omega_fit_(ind);
        tau_fit = tau_fit_(ind);
        w_fit = glass_particles.compute_equidistant_weighting(omega_fit,tau_fit,obj.FB_Bingham_wscheme{obj.exp_id});

        [obj.omega_fit_Bingham obj.tau_fit_Bingham, obj.i_fit_Bingham] = deal(omega_fit,tau_fit,ind);

        [mp_b ty_b] = obj.determine_Bingham_fluid_bounds(omega_fit,tau_fit);

        [obj.mu_p_Bingham obj.tau_y_Bingham] = glass_particles.fit_internal_Bingham_fluid(omega_fit, tau_fit, w_fit, mp_b, ty_b);
        obj.rc_Bingham = glass_particles.determine_rc_Bingham(obj.mu_p_Bingham, obj.tau_y_Bingham, obj.omega);
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
    function t_out =  tau_wall_pred(obj, omega_, mp_, ty_)
        if (nargin==2)
            [mp,ty] = deal(obj.mu_p,obj.tau_y);
        elseif (nargin==3)
            [mp,ty] = deal(mp(1),mp(2));
        else
            [mp,ty] = deal(mp(1),mp(2));
        end

        r_o = obj.r_o;
        r_i = obj.r_i;
        ie2 = (r_o/r_i)*(r_o/r_i);
        ir_i2 = 1.0/(r_i*r_i);
        cr = 2.0*ir_i2;
        rp = ie2*ones(size(omega_));

        sm = @(t1,t2) double(t1>t2).*exp(-1.0./((t1-t2).*(t1-t2)));
        rc_f = @(tw,ty) max((r_i*sqrt(tw/ty)).*sm(ty*rp,tw) + r_o*sm(tw,ty*rp), r_i);

        % tw_pred = @(o,t,rc) cr*(ir_i2 - rc.^(-2)).*(mp*o - ty*log(r_i./rc));
        tw_pred = @(o,t,rc) ty - cr*(ir_i2 - rc.^(-2)).*(mp*o - ty*log(rc/r_i));


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
    function compute_dimensionless_Bingham_new(obj)
        [ri,ro,h] = deal(obj.r_i, obj.r_o, obj.h);
        [o_,mp_,ty_] = deal(obj.omega, obj.mu_p_Bingham, obj.tau_y_Bingham);
        nu = mp_/obj.rho_b;
        eta = ri/ro;
        imp_ = 1/mp_;
        ocrit = 0.5*ty_*imp_*(log(eta*eta) + ((1/(eta*eta)) - 1));
        i_full = reshape(1:length(o_), size(o_));
        i_plug = o_<ocrit;
        i_shear = ~i_plug;

        o_shear = o_(i_shear);
        o_plug = o_(i_plug);

        rc = ro*ones(size(o_));
        ti = glass_particles.taui_pred_Bingham(mp_,ty_,o_);
        rc(i_plug) = ri*sqrt(ti(i_plug)/ty_);
        rc2 = rc.*rc;
        rt = sqrt(ri * rc);
        rt2 = rt.*rt;
        etac = ri./rc;
        etac2 = etac.*etac;
        % obj.gamma_Bingham(i_shear) = imp_*(ty_+2*(o_(i_shear)*mp_+ty_*log(1/eta))/(eta*eta-1));
        obj.gamma_Bingham = imp_*(ti-ty_);
        obj.S_Bingham = imp_*(((2./rt2).*((mp_*o_-ty_*log(etac)))./(1/(ri*ri) - 1./(rc2)))-ty_);

        d = rc - ri;

        obj.Re_s_Bingham = (obj.S_Bingham.*(d.*d))/nu;
        obj.Re_b_Bingham = (obj.omega*ri*(ro-ri))/nu;
        % obj.Re_b_Bingham = (obj.omega*ri.*(rc-ri))/nu;
        obj.cf_Bingham = obj.mu_torque./(2*pi*ri*ri*h*obj.S_Bingham.*obj.S_Bingham.*d.*d);
        obj.G_b_Bingham = (obj.rho_b/(h*mp_*mp_))*obj.mu_torque;
        obj.G_rat_Bingham = (obj.mu_torque.*(rc.*rc-ri*ri))./(4*pi*ri*ri*h*(rc.*rc).*(mp_*o_ + ty_*log(rc/ri)));
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
        string_out = [obj.tag, ' \textit{q}=', num2str(round(obj.q, 1),2)];
      else
        string_out = [obj.tag, ' \textit{q}=', num2str(round(obj.q, 1),'%.1f')];
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
    function gammai_out = comp_gammai_ro(obj,omegai_,par_)
        if (nargin==2)
            [mu_p tau_y] = deal(obj.mu_p, obj.tau_y);
        else
            [mu_p tau_y] = deal(par_(1), par_(2));
        end
        r_o = obj.r_o;
        r_i = obj.r_i;
        gammai_out = (1/mu_p)*(tau_y+2*(omegai_*mu_p+tau_y*log(r_o/r_i))/((r_i*r_i)/(r_o*r_o)-1));
    end
    function taui_out = comp_taui_ro(obj,gammai_,par_)
        if (nargin==2)
            [mu_p tau_y] = deal(obj.mu_p, obj.tau_y);
        else
            [mu_p tau_y] = deal(par_(1), par_(2));
        end
        taui_out=mu_p*abs(gammai_)+tau_y;
    end
    function [gammai_out,taui_out] = comp_gammai_taui_ro(obj,omegai_,par_)
        if (nargin==2)
            par = [obj.mu_p, obj.tau_y];
        else
            par = par_;
        end
        gammai_out=obj.comp_gammai_ro(omegai_,par);
        taui_out=obj.comp_taui_ro(gammai_out,par);
    end
    function [gammai_out,taui_out] = comp_gammai_taui_rc(obj,omegai_, par_)
        if (nargin==2)
            [mu_p tau_y] = deal(obj.mu_p, obj.tau_y);
        else
            [mu_p tau_y] = deal(par_(1), par_(2));
        end

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
    function [rc_out,gammai_out,taui_out,ogt_plug_out] = comp_bingham_shear_rc(obj,omegai_,par_)
        if (nargin==2)
            par = [obj.mu_p, obj.tau_y];
        else
            par = par_;
        end

        [gammai_out taui_out] = obj.comp_gammai_taui_rc(omegai_,par);
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
    function [gammai_out,taui_out,o_g_t_plug_out] = comp_gammai_taui(obj,omegai_,par_)
        if (nargin==2)
            par = [obj.mu_p, obj.tau_y];
        else
            par = par_;
        end
        [gammai_ro taui_ro] = obj.comp_gammai_taui_ro(omegai_,par);
        [rc_vec gammai_out taui_out o_g_t_plug_out] = obj.comp_bingham_shear_rc(omegai_,par);
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

        obj.omega_fit_Carreau = omega_fit;
        obj.tau_fit_Carreau = tau_fit;

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

        mueffi_out(mueffi_out<(1.81e-5)) = 1.81e-5;

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
    function [mp_bout ty_bout] = determine_Bingham_fluid_bounds(obj,o_,t_)
        o_max=max(obj.omega);
        o_min=min(obj.omega);
        t_max=max(obj.tau);
        t_min=min(obj.tau);
        ofit_min = min(obj.omega_fit_Bingham);
        ofit_max = max(obj.omega_fit_Bingham);
        tfit_min = min(obj.tau_fit_Bingham);
        tfit_max = max(obj.tau_fit_Bingham);

        ty_min = obj.FB_Bingham_tauy_bounds(obj.exp_id,1);
        ty_max = obj.FB_Bingham_tauy_bounds(obj.exp_id,3);
        ty_0 = obj.FB_Bingham_tauy_bounds(obj.exp_id,2);

        mp_min = obj.FB_Bingham_mup_bounds(obj.exp_id,1);
        mp_max = obj.FB_Bingham_mup_bounds(obj.exp_id,3);
        mp_0 = obj.FB_Bingham_mup_bounds(obj.exp_id,2);

        mp_bout = [mp_min mp_0 mp_max];
        ty_bout = [ty_min ty_0 ty_max];
    end
    function [mu0_bout lambda_bout n_bout k_bout muinf_bout] = determine_Carreau_fluid_bounds_muinf(obj,o_,t_)
        if (isempty(strfind(obj.name,'49'))) % this is FB1
            mu0_bmat = [    27.9811205731978e-003 502.713381750567e+000 8.21789399425747e+003; ...
                            21.6298507103194e-003 383.998202908695e+000 8.25450147238607e+003; ...
                            19.0084821751595e-003 228.292494505235e+000 8.19918363792479e+003; ...
                            9.47669022855846e-003 113.824173138701e+000 8.19972854742755e+003; ...
                            5.88830957203466e-003 70.6643178782695e+000 8.16801672363917e+003; ...
                            2.61559591993504e-003 31.3478091162148e+000 8.11777836043818e+003; ...
                            1.08304546567496e-003 12.9682312248578e+000 8.03392889127985e+003; ...
                            371.029993630266e-006 4.44886857313514e+000 8.06614604898461e+003; ...
                            174.331007841861e-006 2.09173269801252e+000 8.12000780803754e+003; ...
                            131.420318182700e-006 1.57682289156162e+000 8.11892744334501e+003];
            l_bmat = [      18.3841160293422e-003 36.5941838924198e+000 100.000000000000e+003; ...
                            18.2808753228252e-003 36.5751227574102e+000 100.000000000000e+003; ...
                            18.2801413760622e-003 36.5735291343548e+000 100.000000000000e+003; ...
                            18.2805581433744e-003 36.5759597743373e+000 100.000000000000e+003; ...
                            18.2810557857380e-003 36.5490944912188e+000 100.000000000000e+003; ...
                            18.2807328058862e-003 36.4974656591081e+000 100.000000000000e+003; ...
                            18.2806651005088e-003 36.4656852286121e+000 100.000000000000e+003; ...
                            18.2803284109446e-003 36.5110016473680e+000 100.000000000000e+003; ...
                            18.2805690329284e-003 36.5405875957595e+000 100.000000000000e+003; ...
                            18.2803082831570e-003 36.5424515315612e+000 100.000000000000e+003];
            muinf_bmat = [  181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006; ...
                            181.000000000000e-007 181.000000000000e-007 500.000000000000e-006];
        else % this is FB2
            mu0_bmat = [    481.015626732804e-003 5.80066454635963e+003 8.15145521795111e+003; ...
                            382.350672459460e-003 4.46176139554654e+003 7.34389696454814e+003; ...
                            263.719909570974e-003 3.05143083102082e+003 6.82863424569648e+003; ...
                            192.525450055062e-003 2.37451867214429e+003 6.75804501726294e+003; ...
                            19.3130687125380e-003 466.975921160043e+000 6.08379281463824e+003; ...
                            10.4280686256002e-003 144.183253064103e+000 5.58531304196371e+003; ...
                            2.77508126936174e-003 33.3039172009423e+000 5.25992858614400e+003; ...
                            804.243392569149e-006 9.65434800437134e+000 4.92059458923057e+003; ...
                            373.047801570539e-006 4.47813356703683e+000 4.73924340595501e+003; ...
                            237.677983552893e-006 2.85277294847429e+000 4.58110617460558e+003; ...
                            94.1732680347510e-006 1.13020543320845e+000 4.45718723664347e+003; ...
                            74.0178838662886e-006 888.263466993723e-003 4.45434790507383e+003; ...
                            53.3142786109067e-006 639.707364597078e-003 4.13747949599738e+003];
            l_bmat = [      4.87821007024614e-003 36.5985271415181e+000 100.000000000000e+003; ...
                            5.47357728247160e-003 35.2892943221658e+000 100.000000000000e+003; ...
                            5.47341356466936e-003 35.2892943221658e+000 100.000000000000e+003; ...
                            7.80235937952078e-003 36.5985271415181e+000 100.000000000000e+003; ...
                            12.4798905890739e-003 36.5985271415181e+000 100.000000000000e+003; ...
                            12.4794650463798e-003 36.6005401712275e+000 100.000000000000e+003; ...
                            12.4794650463798e-003 36.6032855686421e+000 100.000000000000e+003; ...
                            15.7827103978257e-003 36.6115242326207e+000 100.000000000000e+003; ...
                            15.7830507111362e-003 36.6128061359111e+000 100.000000000000e+003; ...
                            15.7830507111362e-003 36.6082283220235e+000 100.000000000000e+003; ...
                            15.7830507111362e-003 36.6056652462411e+000 100.000000000000e+003; ...
                            15.7830507111362e-003 36.6051160624094e+000 100.000000000000e+003; ...
                            15.7830507111362e-003 36.5994421275713e+000 100.000000000000e+003];
            muinf_bmat = [  18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006; ...
                            18.1000000000000e-006 18.1000000000000e-006 500.000000000000e-006];
        end

        k_Newt = 2*(obj.r_o*obj.r_o)/(obj.r_o*obj.r_o-obj.r_i*obj.r_i);
        mu0_bout = mu0_bmat(obj.exp_id,:);
        lambda_bout = l_bmat(obj.exp_id,:);
        n_bout = [0 0.5 1];
        k_bout = [k_Newt k_Newt k_Newt];
        muinf_bout = muinf_bmat(obj.exp_id,:);
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
    function rc_out = determine_rc_Bingham(mp_,ty_,o_)
        eta = fluid.r_i_def/fluid.r_o_def;
        ocrit = 0.5*(ty_/mp_)*(log(eta*eta) + ((1/(eta*eta)) - 1));

        i_full = reshape(1:length(o_), size(o_));
        i_plug = o_<ocrit;
        i_shear = ~i_plug;

        o_shear = o_(i_shear);
        o_plug = o_(i_plug);
        rc = fluid.r_o_def*ones(size(o_));
        ti = glass_particles.taui_pred_Bingham(mp_,ty_,o_);
        rc(i_plug) = fluid.r_i_def*sqrt(ti(i_plug)/ty_);
        rc_out = rc;
    end
    function [oc_out tc_out] = get_plug_crit(mp_,ty_)
        eta = fluid.r_i_def/fluid.r_o_def;
        oc_out = 0.5*(ty_/mp_)*(log(eta*eta) + ((1/(eta*eta)) - 1));
        tc_out = ty_/(eta*eta);
    end
    function appmu_out = compute_appmu_Bingham(ty_, o_, t_)
      [ri ro h] = deal(glass_particles.r_i_def, glass_particles.r_o_def, glass_particles.h_def);

      rc = ro*ones(size(t_));
      tc = ty_*((ro*ro)/(ri*ri));

      i_shear = t_>tc;
      i_plug = ~i_shear;

      rc(i_plug) = ri*sqrt(t_(i_plug)/ty_);
      etac = ri./rc;
      appmu_out = (ty_*log(etac) + 0.5*t_.*(1 - etac.*etac) )./o_;
    end

    function ocrit_out = determine_omega_crit_alpha_peak(o_, a_)
        [alpha_max imax] = max(a_);
        ocrit_out = 1.0001*o_(imax);
    end

    function ocrit_out = determine_omega_crit(alpha_tol_, omega_min_, o_, a_)
      ifull = 1:length(o_);
      i_transitioned = logical((o_>omega_min_).*(a_>alpha_tol_));

      if (sum(i_transitioned))
          icrit = min(ifull(i_transitioned));
          % ocrit_out=(alpha_tol_-a_(icrit))*((o_(icrit)-o_(icrit-1))/(a_(icrit)-a_(icrit-1))) + o_(icrit);
          ocrit_out=1.000001*o_(icrit);
      else
          ocrit_out = 1.1*max(o_);
      end
    end

    function w_out = compute_equidistant_weighting(o_,t_,wscheme_)
        if (nargin==3)
            wscheme = wscheme_;
        else
            wscheme = 'none';
        end

        w_out = nan(size(o_));
        if (strcmp(wscheme,'tmagnorm_odist1'))
            for i = 1:length(o_)
                diff = o_-o_(i);
                w_out(i) = sum(abs(diff))/abs(t_(i));
            end
        elseif (strcmp(wscheme, 'tmagnorm'))
            for i = 1:length(o_)
                w_out(i) = 1/abs(t_(i));
            end
        elseif (strcmp(wscheme, 'tmagnorm_odist2'))
            for i = 1:length(o_)
                diff = o_-o_(i);
                w_out(i) = sqrt(sum(diff.*diff))/abs(t_(i));
                % w_out(i) = sum(sqrt(abs(diff)))/abs(t_(i));
            end
        else
            for i = 1:length(o_)
                diff = o_-o_(i);
                % w_out(i) = sqrt(sum(diff.*diff));
                % w_out(i) = sqrt(sum(diff.*diff))/abs(t_(i));
                % w_out(i) = sum(abs(diff));
                % w_out(i) = sum(abs(diff))/abs(t_(i));
                % w_out(i) = sum(sqrt(abs(diff)));
                % w_out(i) = sum(sqrt(abs(diff)))/abs(t_(i));
                % w_out(i) = 1/abs(t_(i));
                w_out(i) = 1;
            end
        end
        w_out = w_out/length(w_out);
        w_out = w_out/norm(w_out);
    end
    function o_out = omegai_pred_Bingham(mp_,ty_,t_)
        o_out = zeros(size(t_));
        eta = fluid.r_i_def/fluid.r_o_def;
        tplug = ty_/(eta*eta);
        i_full = reshape(1:length(t_),size(t_));
        i_noflow = t_<ty_;
        i_plug = logical(double(~i_noflow).*double(t_<tplug));
        i_shear = logical(double(~i_noflow) + double(~i_plug));

        tplug = t_(i_plug);
        tshear = t_(i_shear);

        o_out(i_plug) = 0.5*(ty_/mp_)*(log(ty_./tplug) + (tplug/ty_) - 1);
        o_out(i_shear) = (1/mp_)*(ty_*log(eta) + 0.5*t_*(1-(eta*eta)));
    end
    function t_out = taui_pred_Bingham(mp_,ty_,o_,t0_)
        if (nargin==3)
            t0vec = ones(size(o_));
        else
            t0vec = t0_;
        end
        eta = glass_particles.r_i_def/glass_particles.r_o_def;
        imp_ = 1/mp_;
        ocrit = 0.5*ty_*imp_*(log(eta*eta) + ((1/(eta*eta)) - 1));
        t_out = nan(size(o_));
        i_full = reshape(1:length(o_), size(o_));
        i_plug = o_<ocrit;
        i_shear = ~i_plug;

        gi_shear = imp_*(ty_+2*(o_(i_shear)*mp_+ty_*log(1/eta))/(eta*eta-1));
        t_out(i_shear)=mp_*abs(gi_shear)+ty_;

        gammai_solve = @(t,o) (ty_ + 2*(o*mp_ + 0.5*ty_*log(t/ty_))./((ty_./t)-1))/mp_;
        bingham_solve = @(t,o) (mp_ * abs(gammai_solve(t,o))) + ty_ - t;

        bounds = [ty_,1e3];

        if (isempty(gcp('nocreate')))
            for i = reshape(i_full(i_plug),1,[])
                solve_i = @(t) abs(bingham_solve(t,o_(i)));
                t_out(i) = fminbnd(solve_i,bounds(1),bounds(2));
            end
        else
            parfor i = reshape(i_full(i_plug),1,[])
                solve_i = @(t) abs(bingham_solve(t,o_(i)));
                t_out(i) = fminbnd(solve_i,bounds(1),bounds(2));
            end
        end
    end
    function t_out = taui_pred_Carreau(mu0_,l_,n_,k_,muinf_,o_)
        if (nargin==2)
            [mu0 l n k muinf] = deal(mu0_(1),mu0_(2),mu0_(3),mu0_(4),mu0_(5));
            omega = l_;
        else
            [mu0 l n k muinf] = deal(mu0_,l_,n_,k_,muinf_);
            omega = o_;
        end
        gi = -1.0*(k*omega);
        mueffi = muinf+((mu0-muinf)*((1+((l*gi).*(l*gi))).^(0.5*(n-1))));
        t_out = mueffi.*abs(gi);
    end
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
    function [mu_p_out tau_y_out] = fit_internal_Bingham_fluid(omega_fit, tau_fit, w_fit, mp_b, ty_b)
          trust_region_reflect = 'trust-region-reflective';
          Lev_Marq = 'levenberg-marquardt';
          % functol_try=1e-16;
          % steptol_try=1e-16;
          % optimtol_try=1e-16;
          % max_it_try=1e5;
          % max_eval_try=1e5;

          functol_try=1e-8;
          steptol_try=1e-8;
          optimtol_try=1e-8;
          max_it_try=1000;
          max_eval_try=1000;

          alg_use=Lev_Marq;
          optimtol_use=optimtol_try;
          functol_use=functol_try;
          steptol_use=steptol_try;

          r_o = glass_particles.r_o_def;
          r_i = glass_particles.r_i_def;
          ie2 = (r_o/r_i)*(r_o/r_i);
          ir_i2 = 1.0/(r_i*r_i);
          cr = 2.0*ir_i2;
          rp = ie2*ones(size(tau_fit));

          sm = @(t1,t2) double(t1>t2).*exp(-1.0./((t1-t2).*(t1-t2)));
          rc_f = @(tw,ty) max((r_i*sqrt(tw/ty)).*sm(ty*rp,tw) + r_o*sm(tw,ty*rp), r_i);

          penalty_func = @(mp,ty,o,t,w) w.*(t-glass_particles.taui_pred_Bingham(mp,ty,o,t));

          [mp_min mp_0 mp_max] = deal(mp_b(1), mp_b(2), mp_b(3));
          [ty_min ty_0 ty_max] = deal(ty_b(1), ty_b(2), ty_b(3));

          upper_bound=[mp_max ty_max];
          lower_bound=[mp_min ty_min];
          init_guess=[mp_0 ty_0];

          fit_func = @(x) (penalty_func(x(1),x(2),omega_fit,tau_fit,w_fit));

          % solve_opts=optimoptions(@lsqnonlin,'Display','final-detailed', ...
          % solve_opts=optimoptions(@lsqnonlin,'Display','iter', ...
          solve_opts=optimoptions(@lsqnonlin,'Display','off', ...
                                  'Algorithm',alg_use, ...
                                  'MaxIterations', max_it_try, ...
                                  'MaxFunctionEvaluations', max_eval_try, ...
                                  'FiniteDifferenceType', 'central', ...
                                  'OptimalityTolerance',optimtol_use, ...
                                  'FunctionTolerance',functol_use, ...
                                  'StepTolerance', steptol_use);

          x_out = lsqnonlin(fit_func, init_guess, ...
                                      lower_bound, ...
                                      upper_bound, ...
                                      solve_opts);
          [mu_p_out tau_y_out] = deal(x_out(1),x_out(2));
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
