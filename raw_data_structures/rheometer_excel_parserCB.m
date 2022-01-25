% anton parr excel spreadsheets parsing script. Keep in same directory as
% excel files
function raw = rheometer_excel_parserCB(file, sheet)
if ~ischar(file)
    error('input must be char array of filename')
end
raw = xlsread(file, sheet);
raw = [NaN(2,9); raw];
import = 'success'
end
