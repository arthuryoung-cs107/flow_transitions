classdef phi_experiment < handle
  properties
    name;
    col_scheme;
    filename = './raw_data_structures/anton_parr_data/Arthur_9_8_2018_height_check_Q1_1-2_9_logrpm_0_1-1000';

  end
  methods
    function obj = phi_experiment(name_, col_scheme_) %% fuck it, just do this. 
      obj.name = name_;
      obj.col_scheme = col_scheme_;
      raw = rheometer_height_change_parser(obj.filename, 1);
      q_raw = zeros(1, length(range_array));
      rpm_raw = zeros(1, length(range_array));
      for i = 1:length(range_array)
          q_raw(i) = nanmean(raw(range_array(1,i):range_array(2,i), 2));
          rpm_raw(i) = nanmean(raw(range_array(1,i):range_array(2,i), 3));
      end

      r_i = 0.01208;
      r_o = 0.025;
      h = 0.036;
      Q_inc = 1.2; %% according to notes.txt in the ~/SURF_2018/ directory

      Vol_theo = pi*(r_i)^2*(h)*(10^2)^3;
      Vol_meas = 19.94;

      rpm_0_h = [96 96 96 96 97 100 101.5 104 105 106 107] - Vol_meas; % yields net powder volume (mL). assume that PC is same vol as SC
      rpm_01_h = [96 96 96 96 97 100 102 104 105 106 107] - Vol_meas;
      rpm_1_h = [96 96 96 96 97 100 102 104.5 105 106 107.5] - Vol_meas;
      rpm_10_h = [96 96 96 96 98 100 102 104.5 105.5 106.5 107.5] - Vol_meas;
      rpm_100_h = [96 96 96 96 99 100.5 102.5 104.5 105.5 106.5 108] - Vol_meas;
      rpm_1000_h = [96 97 98 100 101 102 104 105 106 107 107] - Vol_meas;
      rpm0_vec = [0, 0.1, 1, 10, 100, 1000];

      phi_0 = rpm_0_h.^(-1)*111.6/2.5; % yields fraction between solid space (mL) and net powder volume
      phi_01 = rpm_01_h.^(-1)*111.6/2.5;
      phi_1 = rpm_1_h.^(-1)*111.6/2.5;
      phi_10 = rpm_10_h.^(-1)*111.6/2.5;
      phi_100 = rpm_100_h.^(-1)*111.6/2.5;
      phi_1000 = rpm_1000_h.^(-1)*111.6/2.5;
      phi0_mat = [phi_0', phi_01', phi_1', phi_10', phi_100', phi_1000'];

      q_vector = [0 0.5374 1.0385 1.2739 1.5243 1.7879 2.067 2.2810 2.5 2.75 3.0449]/Q_inc; %true flowrate vector, don't think this data was actually truly recorded in .csv

      % raw data, to be matched to rpm and q values
      rpm1_01_h = [99 99 99.5 99.5 99.5 99.5 102 103 104 105 106] - Vol_meas;
      rpm1_1_h = [99 99 99.5 100 101 101 101.5 103 104 105 106] - Vol_meas;
      rpm1_10_h = [99.5 99.5 100 100 100 102 103 104 105 105.5 106] - Vol_meas;
      rpm1_100_h = [100 100.5 101 101 101.5 101.5 103 104 105 105.5 106] - Vol_meas;
      rpm1_1000_h = [100 100 100 100 100 101 101 101 101 101 102] - Vol_meas;

      rpm2_05_h = [100 100 100 101 101 101.5 102 103 104 105 106] - Vol_meas;
      rpm2_5_h = [100.5 100.5 100.5 101 101 101.5 102 103.5 104.5 105.5 106.5] - Vol_meas;
      rpm2_50_h = [100.5 100.5 100.5 100.5 101 101.5 102.5 103.5 105 106.5 107.5] - Vol_meas;
      rpm2_250_h = [100.5 100.5 100.5 100.5 100.5 101 101 102 103.5 104 105] - Vol_meas;
      rpm2_500_h = [100.5 100.5 100.5 100.5 100.5 100.5 101 102 102 102 102] - Vol_meas;
      rpm2_750_h = [100.5 100.5 100.5 100.5 100.5 100.5 102 102.5 103 103 103] - Vol_meas;

      rpm3_01_h = [100 100 100 100.5 100.5 101 102 103.5 105 105.5 106.5] - Vol_meas;
      rpm3_1_h = [100.5 100.5 100.5 101 101.5 102 102.5 103 105 105.5 106.5] - Vol_meas;
      rpm3_10_h = [100.5 100.5 100.5 100.5 101 101 101.5 102.5 103 105.5 106.5] - Vol_meas;
      rpm3_100_h = [100.5 100.5 100.5 101 101.5 102 103 105 105.5 106 107] - Vol_meas;
      rpm3_500_h = [100.5 100.5 100.5 100.5 101 101 101.5 101.5 102 102.5 103] - Vol_meas;
      rpm3_1000_h = [100.5 101 101 101.5 101.5 102 102 102 102 102.5 103] - Vol_meas;

      phi1_01 = rpm1_01_h.^(-1)*111.6/2.5;
      phi1_1 = rpm1_1_h.^(-1)*111.6/2.5;
      phi1_10 = rpm1_10_h.^(-1)*111.6/2.5;
      phi1_100 = rpm1_100_h.^(-1)*111.6/2.5;
      phi1_1000 = rpm1_1000_h.^(-1)*111.6/2.5;
      phi1_mat = [phi1_01', phi1_1', phi1_10', phi1_100', phi1_1000'];

      phi2_05 = rpm2_05_h.^(-1)*111.6/2.5;
      phi2_5 = rpm2_5_h.^(-1)*111.6/2.5;
      phi2_50 = rpm2_50_h.^(-1)*111.6/2.5;
      phi2_250 = rpm2_250_h.^(-1)*111.6/2.5;
      phi2_500 = rpm2_500_h.^(-1)*111.6/2.5;
      phi2_750 = rpm2_750_h.^(-1)*111.6/2.5;
      phi2_mat = [phi2_05', phi2_5', phi2_50', phi2_250', phi2_500', phi2_750'];

      phi3_01 = rpm3_01_h.^(-1)*111.6/2.5;
      phi3_1 = rpm3_1_h.^(-1)*111.6/2.5;
      phi3_10 = rpm3_10_h.^(-1)*111.6/2.5;
      phi3_100 = rpm3_100_h.^(-1)*111.6/2.5;
      phi3_500 = rpm3_500_h.^(-1)*111.6/2.5;
      phi3_1000 = rpm3_1000_h.^(-1)*111.6/2.5;
      phi3_mat = [phi3_01', phi3_1', phi3_10', phi3_100', phi3_500', phi3_1000'];

      trial1_num_rpm = 5;
      trial2_num_rpm = 6;
      trial3_num_rpm = 6;

      trial1_num_q = 11;
      trial2_num_q = 11;
      trial3_num_q = 11;

      trial1_rpm_vec = zeros(1, trial1_num_rpm);
      trial2_rpm_vec = zeros(1, trial2_num_rpm);
      trial3_rpm_vec = zeros(1, trial3_num_rpm);

      trial1_q_vec = zeros(trial1_num_q, trial1_num_rpm);
      trial2_q_vec = zeros(trial1_num_q, trial1_num_rpm);
      trial3_q_vec = zeros(trial1_num_q, trial1_num_rpm);

      trial1_index_range = [1, trial1_num_rpm*trial1_num_q];
      trial2_index_range = [max(trial1_index_range)+1, (max(trial1_index_range)+1) + trial2_num_rpm*trial2_num_q-1];
      trial3_index_range = [max(trial2_index_range)+1, (max(trial2_index_range)+1) + trial3_num_rpm*trial3_num_q-1];

      %
      for i = 1:trial1_num_rpm %average rpm value across q values for each rpm trials 1
          trial1_rpm_vec(i) = mean(rpm_raw( trial1_index_range(1)+(i-1)*(trial1_num_q):(trial1_index_range(1)+(i)*(trial1_num_q)-1) ));
      end
      for i = 1:trial2_num_rpm %average rpm value across q values for each rpm trials 2
          trial2_rpm_vec(i) = mean(rpm_raw( trial2_index_range(1)+(i-1)*(trial2_num_q):(trial2_index_range(1)+(i)*(trial2_num_q)-1) ));
      end
      for i = 1:trial3_num_rpm %average rpm value across q values for each rpm trials 2
          trial3_rpm_vec(i) = mean(rpm_raw( trial3_index_range(1)+(i-1)*(trial3_num_q):(trial3_index_range(1)+(i)*(trial3_num_q)-1) ));
      end

      for i = 1:trial1_num_rpm %average rpm value across q values for each rpm trials 1
          trial1_q_vec(:, i) = q_raw(trial1_index_range(1)+(i-1)*(trial1_num_q):(trial1_index_range(1)+(i)*(trial1_num_q)-1));
      end
      for i = 1:trial2_num_rpm %average rpm value across q values for each rpm trials 2
          trial2_q_vec(:, i) = q_raw(trial2_index_range(1)+(i-1)*(trial2_num_q):(trial2_index_range(1)+(i)*(trial2_num_q)-1));
      end
      for i = 1:trial3_num_rpm %average rpm value across q values for each rpm trials 2
          trial3_q_vec(:, i) = q_raw( trial3_index_range(1)+(i-1)*(trial3_num_q):(trial3_index_range(1)+(i)*(trial3_num_q)-1) );
      end
      trial1_q_vec = trial1_q_vec/Q_inc;
      trial2_q_vec = trial2_q_vec/Q_inc;
      trial3_q_vec = trial3_q_vec/Q_inc;

      % % % combining data for fit

      trial1_rpm_mat = [trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec; trial1_rpm_vec];
      trial2_rpm_mat = [trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec; trial2_rpm_vec];
      trial3_rpm_mat = [trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec; trial3_rpm_vec];

      trial1_omega_vec = 2*pi/60* trial1_rpm_vec;
      trial2_omega_vec = 2*pi/60* trial2_rpm_vec;
      trial3_omega_vec = 2*pi/60* trial3_rpm_vec;

      rpm_mat_all = [trial1_rpm_mat, trial2_rpm_mat, trial3_rpm_mat];

      q_fit = [trial1_q_vec(:); trial2_q_vec(:); trial3_q_vec(:)];
      rpm_fit = [trial1_rpm_mat(:); trial2_rpm_mat(:); trial3_rpm_mat(:)];
      omega_fit = 2*pi/60*rpm_fit;
      phi_fit = [phi1_mat(:); phi2_mat(:); phi3_mat(:)];

      q_mat_all = [trial1_q_vec, trial2_q_vec, trial3_q_vec];
      phi_mat_all = [phi1_mat, phi2_mat, phi3_mat];

      sf = fit([q_fit, rpm_fit],phi_fit,'poly55');


      figure('Name', 'SUBFig phi vs. q vs. omegai' ,'Renderer', 'painters', 'Position', [0 500 580 325])
      view([135, 45]);
      hold on
      % for i = 1:6
      %     plot3(q_vector, rpm0_vec(i)*ones(1, length(q_vector)), phi0_mat(:, i), 'o - k','LineWidth', 2, 'DisplayName', 'trial 0')
      % end
      % for i = 1:trial1_num_rpm
      %     scatter3(trial1_q_vec(:, i), trial1_omega_vec(i)*ones(1, length(trial1_q_vec)), phi1_mat(:, i), '*','LineWidth', 2, 'DisplayName', 'trial 1')
      % end
      % for i = 1:trial2_num_rpm
      %     scatter3(trial2_q_vec(:, i), trial2_omega_vec(i)*ones(1, length(trial2_q_vec)), phi2_mat(:, i), '*','LineWidth', 2, 'DisplayName', 'trial 1')
      % end
      % for i = 1:trial3_num_rpm
      %     scatter3(trial3_q_vec(:, i), trial3_omega_vec(i)*ones(1, length(trial3_q_vec)), phi3_mat(:, i), '*','LineWidth', 2, 'DisplayName', 'trial 1')
      % end
      % colorbar
      % surf(trial1_q_vec, trial1_rpm_mat, phi1_mat)
      % surf(trial2_q_vec, trial2_rpm_mat, phi2_mat)
      % surf(trial3_q_vec, trial3_rpm_mat, phi3_mat)

      % plot(sf, [q_fit,rpm_fit],phi_fit)
      scatter3(q_fit(:), omega_fit(:), phi_fit(:), ' *','LineWidth', 1, 'CData', phi_fit(:))
      colorbar
      zlabel('$$\phi$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
      ylabel('$$\omega_i [rad/s]$$', 'Interpreter', 'LaTeX','FontSize',12)
      xlabel('$$q = \frac{Q}{Q_{inc}}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)

      rpm_500_all = [rpm_mat_all(:, 10); rpm_mat_all(:, 16)];
      omega_50_all = 2*pi/60*rpm_500_all;
      q_50_all = [q_mat_all(:, 10); q_mat_all(:, 16)];
      phi_50_all = [phi_mat_all(:, 10); phi_mat_all(:, 16)];


      omega50_phi_vs_q_fit = fit([q_50_all(6:11); q_50_all(15:22)],[phi_50_all(6:11); phi_50_all(15:22)],'poly1');
      figure('Name', 'SUBFig phi vs. q vs. omegai' ,'Renderer', 'painters', 'Position', fig_positions(3, :))
      hold on
      plot(q_50_all(1:11), phi_50_all(1:11), ' * r','LineWidth', 1)
      plot(q_50_all(12:22), phi_50_all(12:22), ' * b','LineWidth', 1)
      fplot( @(q) (omega50_phi_vs_q_fit.p1)*(q) + omega50_phi_vs_q_fit.p2, [1.4 2.6],'--', 'Color', [0 0 0],'Linewidth', 2, 'DisplayName', 'Couette flow')
      ylabel('$$\phi$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)
      xlabel('$$q = Q/Q_{inc}$$ [dimensionless]', 'Interpreter', 'LaTeX','FontSize',12)


    end

  end
end

function raw = rheometer_height_change_parser(file, sheet)
  if ~ischar(file)
      error('input must be char array of filename')
  end
  raw = xlsread(file, sheet);
  raw = [NaN(2,5); raw];
  import = 'success'
end

function range_array = range_finder_height_change( raw_array )
% given a standard 'raw' matrix, constructs 2xn array of values consisting
% of the range of each swept rpm measurements

range_array = [7; 0];
index_check = 0;


for i = 7:length(raw_array)
%     checks if current i is NaN
    if index_check == 0
        if isnan(raw_array(i, 5)) == 1 % if beginning of gap or invalid value
            if isnan(raw_array(i+1, 5)) == 1 && isnan(raw_array(i+2,5)) == 1 && (raw_array(i+7) == 1 || raw_array(i+4) == 1) % if gap, begins skip routine
%             if raw_array(i+7, 1) == 1 % if gap, because new data set, begins skip routine
                index_check = 1;
                continue
            else % if invalid value, treats same as if it is not equal to nan
                dimensions = size(range_array);
                range_array(2, dimensions(2)) = i;
                continue
            end

        else % if this is not equal to nan, regular data point
            dimensions = size(range_array);
            range_array(2, dimensions(2)) = i;
            continue
        end
    else
        if isnan(raw_array(i+1, 5)) == 0 %if next value detected to be data point
            index_check = 0;
            new_column = [i + 1; 0];
            range_array = [range_array, new_column];
            continue
        else
            continue
        end
    end
end
