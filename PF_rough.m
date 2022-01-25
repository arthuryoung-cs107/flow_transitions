classdef PF_rough < super_experiment
  properties
    PF2_in;
    PF3_in;
  end
  methods
    function obj = PF_rough(PF2_in_, PF3_in_)
      obj@super_experiment({PF2_in_; PF3_in_});
      obj.mu_torque(15) = NaN;
      obj.sigma_torque(15) = NaN;
      obj.G(15) = NaN;
      obj.G_rat(15) = NaN;
      obj.cf(15) = NaN;
      obj.omega(15) = NaN;
      obj.Re_s(15) = NaN;

      obj.mu_torque = rmmissing(obj.mu_torque);
      obj.sigma_torque = rmmissing(obj.sigma_torque);
      obj.G = rmmissing(obj.G);
      obj.G_rat = rmmissing(obj.G_rat);
      obj.cf = rmmissing(obj.cf);
      obj.omega = rmmissing(obj.omega);
      obj.Re_s = rmmissing(obj.Re_s);

      run figure_properties.m
      obj.PF2_in = PF2_in_;
      obj.PF3_in = PF3_in_;

      obj.TV_range = obj.Re_s > 70.0;
      obj.color = red5;
      obj.specs = ' o';
      obj.label = 'PF2';
    end
  end
end
