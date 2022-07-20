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
  end
end
