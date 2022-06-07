classdef FB_experiment < experiment
  properties
    statement;
    ind_sort;
    TV_lowomega = 11;
  end
  methods
    function obj = FB_experiment(exp_list_in, label_, color_, specs_, statement_)
      obj@experiment(exp_list_in);
      obj.label = label_;
      obj.color = color_;
      obj.specs = specs_;
      obj.TV_lowRes = 50;

      obj.statement = statement_;
    end
    function process_raws(obj, raws)
      for i=1:obj.len
        obj.exp(i).specs = obj.specs;
        obj.exp(i).process_raw(raws{i});
        obj.exp(i).sort_dimensionless;
        obj.exp(i).alpha = obj.exp(i).alpha_comp();
      end
    end
    function inspect_mu_plastic_fit(obj)
      fig_specs = AYfig.specs_gen(obj.label, obj.def_pos);
      fig_out = AYfig.figure(fig_specs);
      hold on
      set(gca, 'YScale', 'log')
      set(gca, 'XScale', 'log')
      for i=1:length(obj.exp)
        plot(obj.exp(i).omega, obj.exp(i).appmu, ' o','Color', obj.exp(i).color, 'LineWidth', 1.5, 'DisplayName', obj.exp(i).name)
        plot(obj.exp(i).sim_omega, obj.exp(i).sim_appmu, '-','Color', obj.exp(i).color, 'LineWidth', 1.5, 'DisplayName', obj.exp(i).name)
      end
      ylabel('$$\mu_{app}$$ [Pa.s]', 'Interpreter', 'LaTeX','FontSize',12)
      xlabel('$$\omega_{i}$$ [rad/s]', 'Interpreter', 'LaTeX','FontSize',12)
      % legend('Show')
      hold off
    end
    function gen_powerfit(obj)
        r_i = 0.01208;
        r_o = 0.025;

        for i=1:obj.len
          obj.exp(i).gen_powerfit;
        end

        G = obj.exp(1).G;
        Re_s = obj.exp(1).Re_s;
        G_TV = obj.exp(1).G_TV;
        Re_s_TV = obj.exp(1).Re_s_TV;
        for i=2:obj.len
            G = [G; obj.exp(i).G];
            Re_s = [Re_s; obj.exp(i).Re_s];
            G_TV = [G_TV; obj.exp(i).G_TV];
            Re_s_TV = [Re_s_TV; obj.exp(i).Re_s_TV];
        end
        [obj.Re_s I_R] = sort(Re_s);
        obj.G = G(I_R);
        [obj.Re_s_TV I_R_TV] = sort(Re_s_TV);
        obj.G_TV = G_TV(I_R_TV);



        obj.powerfit = fit(rmmissing(obj.Re_s_TV), rmmissing(obj.G_TV),'b*x^m', 'StartPoint', [70, 1]);
        alpha = obj.powerfit.m;
        beta = obj.powerfit.b;
        m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
        obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
        obj.G_c2 = beta*(obj.Re_sc2)^alpha;
    end
  end
end
