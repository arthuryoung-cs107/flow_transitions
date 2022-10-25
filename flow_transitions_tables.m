classdef flow_transitions_tables
    properties (Constant)
    save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];

    end
    properties

    end
    methods
        function obj = flow_transitions_tables()

        end
    end
    methods (Static)
        function write_FB_fit_results(FB1,FB2,name_)
            FB1_len = length(FB1.exp);
            FB2_len = length(FB2.exp);
            FB_len = length(FB1.exp)+length(FB2.exp);

            % dat_width=7;
            dat_width=5;

            eng_precision = '%.3e';
            percent_precision = '%.1f';

            % raw_collabels = {'\mu_0', '\lambda', 'n', 'k', '\mu_{\infty}', '\tau_y', '\mu_p'};
            raw_collabels = {'\mu_0', '\lambda', 'n','\tau_y', '\mu_p'};
            raw_rowlabels = cell(1,FB_len);

            data_raw = nan(FB_len, dat_width);
            j = 0;
            for i = 1:FB1_len
                j = j+1;
                fb = FB1.exp(i);
                % data_raw(j,:) = [fb.mu0_Carreau fb.lambda_Carreau fb.n_Carreau fb.k_Carreau fb.muinf_Carreau fb.tau_y fb.mu_p];
                data_raw(j,:) = [fb.mu0_Carreau fb.lambda_Carreau fb.n_Carreau fb.tau_y fb.mu_p];
                raw_rowlabels{j} = fb.label;
            end
            for i = 1:FB2_len
                j = j+1;
                fb = FB2.exp(i);
                % data_raw(j,:) = [fb.mu0_Carreau fb.lambda_Carreau fb.n_Carreau fb.k_Carreau fb.muinf_Carreau fb.tau_y fb.mu_p];
                data_raw(j,:) = [fb.mu0_Carreau fb.lambda_Carreau fb.n_Carreau fb.tau_y fb.mu_p];
                raw_rowlabels{j} = fb.label;
            end

            tbl_input.tableColLabels = ltx_array(raw_collabels);
            tbl_input.tableRowLabels = ltx_array(raw_rowlabels);
            tbl_input.data = data_raw;

            tbl_input.dataFormatMode = 'column';
            tbl_input.dataFormat = {eng_precision,dat_width};
            tbl_input.tableColumnAlignment = 'c';
            % tbl_input.tableCaption = fill_mhevent_cap(cap_, fig_name_, solve_ref_, mh.nbeads, mhend_.gen_count, final_frame, mhend_.event_block_count);
            tbl_input.tableLabel = name_;
            tbl_input.makeCompleteLatexDocument = 0;

            tbl_latex = latexTable(tbl_input);

            %% post processing to clean up
            % tbl_latex{1,1} = [tbl_latex{1,1} '[ht]'];
            % data_row_indices=7:2:((2*(mhend_.event_block_count-1))+7);
            % tbl_latex = apply_num_percent_data(tbl_latex, data_row_indices, tbl_input.dataFormat, raw_data, percent_precision);

            file_name = [flow_transitions_tables.save_dir tbl_input.tableLabel '.tex'];
            file = fopen(file_name, 'w');
            for i=1:size(tbl_latex, 1)
                fprintf(file, '%s\n', tbl_latex{i,:});
            end
            fclose(file);
        end
    end
end


function fixed_row = apply_num_data(row_, format_, data_, prefix_, suffix_)
    if (nargin == 3)
        prefix = '';
        suffix = '';
    else
        prefix = prefix_;
        suffix = suffix_;
    end

    fixed_row = row_;
    j=1;
    for itype = 2:2:length(format_)
        type=format_{1,itype-1};
        for i = 1:format_{1,itype}
            num_old = sprintf(type, data_(j));
            fixed_row = strrep(fixed_row, [prefix num_old suffix], [prefix '\num{' num_old '}' suffix]);
            j = j+1;
        end
    end
    fixed_row = strrep(fixed_row, '\num{\num{', '\num{');
    fixed_row = strrep(fixed_row, '}}', '}');
end

function fixed_rows = apply_num_percent_data(rows_, indices_, format_, data_, percent_type_)
    fixed_rows = rows_;
    i = 1;
    for irow = indices_
        fixed_rows{irow,1} = apply_num_data(fixed_rows{irow,1}, format_, data_(i,:), ' ', ' ');
        for j = 1:size(data_, 2)
            num_old = ['\num{' sprintf(percent_type_, data_(i,j)) '}'];
            fixed_rows{irow,1} = strrep(fixed_rows{irow,1}, num_old, [num_old ' \%']);
        end
        fixed_rows{irow,1} = strrep(fixed_rows{irow,1}, '\% \%', '\%');
        i = i + 1;
    end
end


function out = ltx_array(in_)
    out = in_;
    for i = 1:length(in_)
        out{i} = ['$' in_{i} '$'];
    end
end

function out = ltx(in_)
    out = ['$' in_ '$'];
end
