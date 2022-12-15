classdef FB_experiment < experiment
  properties
    statement;
    ind_sort;
  end
  methods
    function obj = FB_experiment(exp_list_in, label_, color_, specs_, statement_)
      obj@experiment(exp_list_in);
      obj.label = label_;
      obj.color = color_;
      obj.specs = specs_;

      obj.statement = statement_;
    end
    function par_out = fitted_Bingham_pars(obj)
        par_out = nan(obj.len, 2);
        for i = 1:obj.len
            par_out(i,:) = [obj.exp(i).mu_p_Bingham obj.exp(i).tau_y_Bingham];
        end
    end
    function process_raws(obj, raws)
      for i=1:obj.len
        obj.exp(i).specs = obj.specs;
        obj.exp(i).process_raw(raws{i},i);
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

        Re_s_TV = reshape(obj.exp(1).Re_s_TV,[],1);
        G_TV = reshape(obj.exp(1).G_TV,[],1);
        alpha_TV = reshape(obj.exp(1).alpha_TV,[],1);
        for i=2:obj.len
            Re_s_TV = [Re_s_TV; reshape(obj.exp(i).Re_s_TV,[],1)];
            G_TV = [G_TV; reshape(obj.exp(i).G_TV,[],1)];
            alpha_TV = [alpha_TV; reshape(obj.exp(i).alpha_TV,[],1)];
        end
        Re_s_TV = rmmissing(Re_s_TV);
        G_TV = rmmissing(G_TV);
        alpha_TV = rmmissing(alpha_TV);
        [obj.Re_s_TV I_R_TV] = sort(Re_s_TV);
        obj.G_TV = G_TV(I_R_TV);
        obj.alpha_TV = alpha_TV(I_R_TV);

        obj.powerfit = fit(obj.Re_s_TV, obj.G_TV,'b*x^m', 'StartPoint', [70, 1]);
        alpha = obj.powerfit.m;
        beta = obj.powerfit.b;
        m = (2*pi*r_i*r_o)/((r_o-r_i)^2);
        obj.Re_sc2 = exp((1/(alpha-1))*log(m/beta));
        obj.G_c2 = beta*(obj.Re_sc2)^alpha;
    end
  end
end
