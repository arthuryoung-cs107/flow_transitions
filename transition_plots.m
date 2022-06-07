classdef transition_plots < PF_NUXB_transition_plots & FB_transition_plots
    properties

    end
    methods
        function obj = transition_plots(write_figs_, write_all_figs_, figs_to_write_)
            obj@PF_NUXB_transition_plots(write_figs_, write_all_figs_, figs_to_write_);
            obj@FB_transition_plots(write_figs_, write_all_figs_, figs_to_write_);
        end
    end
end
