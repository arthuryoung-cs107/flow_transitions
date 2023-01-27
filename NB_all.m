classdef NB_all < super_experiment
  properties
    NB1_in;
    NB2_in;
    NB3_in;
  end
  methods
    function obj = NB_all(NB1_in_, NB2_in_, NB3_in_)
      obj@super_experiment({NB1_in_; NB2_in_; NB3_in_});

      run figure_properties.m
      obj.NB1_in = NB1_in_;
      obj.NB2_in = NB2_in_;
      obj.NB3_in = NB3_in_;

      obj.color = green3;
      obj.specs = ' s';
      obj.label = 'NB';
    end
    function gen_powerfit(obj,exp_list_)
        Res_TV_all = reshape(exp_list_{1}.Re_s_TV,[],1);
        G_TV_all = reshape(exp_list_{1}.G_TV,[],1);
        Grat_TV_all = reshape(exp_list_{1}.Grat_TV,[],1);
        alpha_TV_all = reshape(exp_list_{1}.alpha_TV,[],1);
        for i = 2:length(exp_list_)
            Res_TV_all = [ Res_TV_all; reshape(exp_list_{i}.Re_s_TV,[],1)];
            G_TV_all = [ G_TV_all; reshape(exp_list_{i}.G_TV,[],1)];
            Grat_TV_all = [ Grat_TV_all; reshape(exp_list_{i}.Grat_TV,[],1)];
            alpha_TV_all = [ alpha_TV_all; reshape(exp_list_{i}.alpha_TV,[],1)];
        end
        [obj.Re_s_TV, isort] = mink(Res_TV_all, length(Res_TV_all));
        obj.G_TV = G_TV_all(isort);
        obj.Grat_TV = Grat_TV_all(isort);
        obj.alpha_TV = alpha_TV_all(isort);

        Res_fit = reshape(obj.Re_s_TV,[],1);
        G_fit = reshape(obj.G_TV,[],1);
        Grat_fit = reshape(obj.Grat_TV,[],1);

        obj.powerfit = fit(Res_fit,G_fit,'b*x^m', 'StartPoint', [70, 1]);
        obj.powerfit_Grat_Res = fit(Res_fit,Grat_fit,'b*x^m', 'StartPoint', [1/(sqrt(70)), 0.5]);

        alpha = obj.powerfit.m;
        beta = obj.powerfit.b;
        [r_i r_o] = deal(fluid.r_i_def,fluid.r_o_def);
        m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
        obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
        obj.G_c2 = beta*(obj.Re_sc2)^alpha;
    end
  end
end
