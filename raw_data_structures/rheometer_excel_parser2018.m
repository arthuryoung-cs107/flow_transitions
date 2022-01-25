% anton parr excel spreadsheets parsing script. Keep in same directory as
% excel files
function raw = rheometer_excel_parser2018(file, sheet)
if ~ischar(file)
    error('input must be char array of filename')
end
raw = xlsread(file, sheet);
raw = [NaN(2,8); raw];
import = 'success'
end
