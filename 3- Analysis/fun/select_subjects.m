function in = select_subjects(T, spec)
%SELECT_SUBJECTS Rows of the subject table that meet every condition of SPEC.
%   in = select_subjects(T, spec) with spec an n x 3 cell {column, operator,
%   value; ...}, e.g. {'group', '==', 'patient'; 'CPZ_equival', '>', 0}.
%   Operators: '==', '~=', '>', '>=', '<', '<=' (text compared without case),
%   'in' / 'notin' (value: cell of values), 'missing', '~missing'.
%   Missing values never meet a comparison. Any group or subgroup of any
%   table can be compared.
if isempty(spec), in = true(height(T), 1); return; end
if size(spec, 2) ~= 3, error('select_subjects:spec', 'Conditions are {column, operator, value} rows'); end
in = true(height(T), 1);
for k = 1:size(spec, 1)
    [col, op, val] = spec{k, :};
    if ~ismember(col, T.Properties.VariableNames)
        error('select_subjects:column', 'Column %s is not in the subject table', col);
    end
    v = T.(col);
    if isnumeric(v) || islogical(v)
        v = double(v);  miss = isnan(v);
        if ischar(val) || isstring(val) || iscell(val)      % '1' is the number 1, not the code of '1'
            val = str2double(string(val));
            if any(isnan(val)) && ~any(strcmp(op, {'missing', '~missing'}))
                error('select_subjects:value', 'Column %s is numeric; condition value is not a number', col);
            end
        end
    else
        v = lower(strtrim(string(v)));  miss = ismissing(v) | v == "";
        if ischar(val) || isstring(val), val = lower(strtrim(string(val))); end
        if iscell(val), val = lower(strtrim(string(val))); end
    end
    switch op
        case '==',  c = v == val;
        case '~=',  c = v ~= val;
        case '>',   c = v > val;
        case '>=',  c = v >= val;
        case '<',   c = v < val;
        case '<=',  c = v <= val;
        case 'in',  c = ismember(v, val);
        case 'notin', c = ~ismember(v, val);
        case 'missing', c = miss;  miss = false(size(miss));
        case '~missing', c = ~miss;
        otherwise, error('select_subjects:op', 'Unknown operator %s', op);
    end
    in = in & c(:) & ~miss(:);
end
end
