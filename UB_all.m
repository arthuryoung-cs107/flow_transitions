classdef UB_all < super_experiment
  properties
    UB1_in;
    UB2_in;
  end
  methods
    function obj = UB_all(UB1_in_, UB2_in_)
      obj@super_experiment({UB1_in_; UB2_in_});

      run figure_properties.m
      obj.UB1_in = UB1_in_;
      obj.UB2_in = UB2_in_;

      obj.color = blue1;
      obj.specs = ' v';
      obj.label = 'UB';
    end
  end
end
