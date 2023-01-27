classdef PF_rough < super_experiment
    properties
        PF2_in;
        PF3_in;
    end
    methods
        function obj = PF_rough(PF2_in_, PF3_in_)
            obj@super_experiment({PF2_in_; PF3_in_});
            obj.alpha_tol = 1.4;

            % bad data point
            obj.mu_torque(15) = NaN;
            obj.sigma_torque(15) = NaN;
            obj.G(15) = NaN;
            obj.G_rat(15) = NaN;
            obj.cf(15) = NaN;
            obj.omega(15) = NaN;
            obj.Re_s(15) = NaN;

            %% fixing duplicate data points
            obj.mu_torque(17) = 0.5*(obj.mu_torque(17)+obj.mu_torque(18));
            obj.sigma_torque(17) = 0.5*(obj.sigma_torque(17)+obj.sigma_torque(18));
            obj.G(17) = 0.5*(obj.G(17)+obj.G(18));
            obj.G_rat(17) = 0.5*(obj.G_rat(17)+obj.G_rat(18));
            obj.cf(17) = 0.5*(obj.cf(17)+obj.cf(18));
            obj.omega(17) = 0.5*(obj.omega(17)+obj.omega(18));
            obj.Re_s(17) = 0.5*(obj.Re_s(17)+obj.Re_s(18));
            obj.mu_torque(18)= NaN;
            obj.sigma_torque(18)= NaN;
            obj.G(18)= NaN;
            obj.G_rat(18)= NaN;
            obj.cf(18)= NaN;
            obj.omega(18)= NaN;
            obj.Re_s(18)= NaN;
            obj.mu_torque(20) = 0.5*(obj.mu_torque(20)+obj.mu_torque(21));
            obj.sigma_torque(20) = 0.5*(obj.sigma_torque(20)+obj.sigma_torque(21));
            obj.G(20) = 0.5*(obj.G(20)+obj.G(21));
            obj.G_rat(20) = 0.5*(obj.G_rat(20)+obj.G_rat(21));
            obj.cf(20) = 0.5*(obj.cf(20)+obj.cf(21));
            obj.omega(20) = 0.5*(obj.omega(20)+obj.omega(21));
            obj.Re_s(20) = 0.5*(obj.Re_s(20)+obj.Re_s(21));
            obj.mu_torque(21)= NaN;
            obj.sigma_torque(21)= NaN;
            obj.G(21)= NaN;
            obj.G_rat(21)= NaN;
            obj.cf(21)= NaN;
            obj.omega(21)= NaN;
            obj.Re_s(21)= NaN;
            obj.mu_torque(23) = 0.5*(obj.mu_torque(23)+obj.mu_torque(24));
            obj.sigma_torque(23) = 0.5*(obj.sigma_torque(23)+obj.sigma_torque(24));
            obj.G(23) = 0.5*(obj.G(23)+obj.G(24));
            obj.G_rat(23) = 0.5*(obj.G_rat(23)+obj.G_rat(24));
            obj.cf(23) = 0.5*(obj.cf(23)+obj.cf(24));
            obj.omega(23) = 0.5*(obj.omega(23)+obj.omega(24));
            obj.Re_s(23) = 0.5*(obj.Re_s(23)+obj.Re_s(24));
            obj.mu_torque(24)= NaN;
            obj.sigma_torque(24)= NaN;
            obj.G(24)= NaN;
            obj.G_rat(24)= NaN;
            obj.cf(24)= NaN;
            obj.omega(24)= NaN;
            obj.Re_s(24)= NaN;

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

            obj.color = red5;
            obj.specs = ' o';
            obj.label = 'PL2';
        end
        function gen_powerfit(obj)
            obj.gen_powerfit_internal(70);
        end
    end
end
