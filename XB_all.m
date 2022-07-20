classdef XB_all < super_experiment
  properties
    XB1_in;
    XB2_in;
  end
  methods
    function obj = XB_all(XB1_in_, XB2_in_)
      obj@super_experiment({XB1_in_; XB2_in_});

      run figure_properties.m
      obj.XB1_in = XB1_in_;
      obj.XB2_in = XB2_in_;

      obj.color = orange2;
      obj.specs = ' ^';
      obj.label = 'XB';
    end
  end
end
