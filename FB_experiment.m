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
        obj.exp(i).alpha = obj.exp(i).alpha_comp();
        obj.exp(i).Re_s_alpha = obj.exp(i).Re_s_alpha_comp();
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
      omega_all = obj.exp(1).omega;
      G_all = obj.exp(1).G;
      Re_s_all = obj.exp(1).Re_s;
      for i=2:obj.len
        omega_all = [omega_all; obj.exp(i).omega];
        G_all = [G_all; obj.exp(i).G];
        Re_s_all = [Re_s_all; obj.exp(i).Re_s];
      end
      [obj.omega, obj.ind_sort] = mink(omega_all, length(omega_all));
      obj.G = G_all(obj.ind_sort);
      obj.Re_s = Re_s_all(obj.ind_sort);

      % obj.TV_range = obj.Re_s > obj.TV_lowRes;
      obj.TV_range = obj.omega > obj.TV_lowomega;

      obj.powerfit = fit(obj.Re_s(obj.TV_range), obj.G(obj.TV_range),'b*x^m', 'StartPoint', [70, 1]);
      for i=1:obj.len
        obj.exp(i).gen_powerfit; 
      end
    end
  end
end
