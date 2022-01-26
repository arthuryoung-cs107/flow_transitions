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
    function process_raws(obj, raws)
      for i=1:obj.len
        obj.exp(i).specs = obj.specs;
        obj.exp(i).process_raw(raws{i});
      end
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

      obj.TV_range = obj.Re_s > 50.0;

      obj.powerfit = fit(obj.Re_s(obj.TV_range), obj.G(obj.TV_range),'b*x^m', 'StartPoint', [70, 1]);
    end
  end
end
