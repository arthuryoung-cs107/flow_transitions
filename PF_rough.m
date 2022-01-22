classdef PF_rough < super_experiment
  properties
    PF2_in;
    PF3_in;
  end
  methods
    function obj = PF_rough(PF2_in_, PF3_in_)
      obj@super_experiment({PF2_in_; PF3_in_});

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
